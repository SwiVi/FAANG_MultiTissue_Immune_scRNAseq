library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)

## Read in data
seu1 <- readRDS('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/QCed/FAANG1ScrubbedSeurat.rds')
seu2 <- readRDS('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/QCed/FAANG2ScrubbedSeurat.rds')

#See how many cells/genes in Seurat object:
seu1
seu2

#See how many cells in each sample:
table(seu1$SampleID)
table(seu2$SampleID)

# Merge Seurat objects:
seu1$Experiment <- rep('FAANG1', ncol(seu1))
seu2$Experiment <- rep('FAANG2', ncol(seu2))
seu <- merge(seu1, y = seu2)

# Add more meta data info
##Sample+batch info
seu$ID <- paste(seu$SampleID, seu$Experiment, sep = '_')
table(seu$ID)
##10X chemistry info
seu$Chemistry <- seu$ID
Idents(seu) <- seu$Chemistry
levels(seu)
v <- c('v2', 'v2', 'v2', 'v3', 'v2', 'v2', 'v3', 'v3', 'v3', 'v3', 'v3')
names(v) <- levels(seu)
seu <- RenameIdents(seu, v)
seu$Chemistry <- Idents(seu)
##Tissue info
seu$tissue <- gsub('.{2}$', '', seu$SampleID)
##Pig info
seu$PigID <- substr(seu$SampleID, nchar(seu$SampleID) - 2 + 1, nchar(seu$SampleID))

#Split Seurat object into list by sample ID:
seu.list <- SplitObject(seu, split.by = "ID") # split by sample IDs
rm(seu, seu1, seu2)

#See how many cells/genes in each listed Seurat object:
seu.list

## Normalize data using SCT method

#Performed SCTransform on each sample from list:
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

# subset Seurat object list (will need for individual tissue integrations):
seu.list # see list components and which belong to each tissue type
seu.list.thy <- c(seu.list[[7]], seu.list[[8]])
seu.list.bm <- c(seu.list[[1]], seu.list[[2]], seu.list[[9]], seu.list[[10]])
seu.list.sp <- c(seu.list[[5]], seu.list[[6]])
seu.list.ln <- c(seu.list[[3]], seu.list[[4]], seu.list[[11]])

## All samples ----
#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features)
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT")
rm(seu.features, seu.anchors)

#### Calculate principle components
#Calculate 50 PCs:
seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(seu,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
# Tinker with distance and spread parameters to alter how UMAP looks. Example here: https://smorabit.github.io/blog/2020/umap/
seu <- RunUMAP(seu, 
               reduction = "pca",
               dims = 1:30,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

seu <- RunTSNE(seu, 
               dims = 1:30,
               k.seed = 777,
               reduction = "pca") # create UMAP

#Visualize UMAP:
DimPlot(seu,
        reduction = 'umap',
        group.by = 'ID',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        split.by = 'ID',
        shuffle = TRUE,
        group.by = 'ID')

DimPlot(seu,
        reduction = 'umap',
        group.by = 'Chemistry',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'tissue',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'PigID',
        shuffle = TRUE)

FeaturePlot(seu, 
            features = c('nCount_SCT', 'nFeature_SCT', 'prcntMito'),
            reduction = 'umap')

DefaultAssay(seu) <- 'SCT'
FeaturePlot(seu, 
            features = c('AIF1', 'CSF1R', 'CD14',
                         'FLT3',
                         'CD3E', 'CD3G', 'CD4', 'CD8B', 'TRDC', 'CD8A', 'CD2',
                         'CD19', 'CD79B', 'JCHAIN', 'MZB1', 'IRF4',
                         'PCLAF', 'PCNA', 
                         'DNTT', 'RAG1', 'RAG2'),
            reduction = 'umap')

#### Add normalized/scaled data to RNA assay
#dim(seu[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu <- NormalizeData(seu,  # normalize the RNA counts data per cell
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
seu <- ScaleData(seu, # scale the RNA counts data relative to other cells
                 assay = "RNA")
seu <- ScaleData(seu, # scale the SCT counts data relative to other cells
                 assay = "SCT")
#dim(seu[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#### Save the Seurat object
SaveH5Seurat(seu, 
             filename = "/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_AllCells.h5Seurat", overwrite = TRUE)
rm(seu, seu.list)

## Thymus only ----
#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list.thy, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list.thy <- PrepSCTIntegration(seu.list.thy, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list.thy, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features)
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT")
rm(seu.features, seu.anchors)

#### Calculate principle components
#Calculate 50 PCs:
seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(seu,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
# Tinker with distance and spread parameters to alter how UMAP looks. Example here: https://smorabit.github.io/blog/2020/umap/
seu <- RunUMAP(seu, 
               reduction = "pca",
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

seu <- RunTSNE(seu, 
               dims = PCdims,
               k.seed = 777,
               reduction = "pca") # create UMAP

#Visualize UMAP:
DimPlot(seu,
        reduction = 'umap',
        group.by = 'ID',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        split.by = 'ID',
        shuffle = TRUE,
        group.by = 'ID')

DimPlot(seu,
        reduction = 'umap',
        group.by = 'Chemistry',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'tissue',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'PigID',
        shuffle = TRUE)

FeaturePlot(seu, 
            features = c('nCount_SCT', 'nFeature_SCT', 'prcntMito'),
            reduction = 'umap')

DefaultAssay(seu) <- 'SCT'
FeaturePlot(seu, 
            features = c('AIF1', 'CSF1R', 'CD14',
                         'FLT3',
                         'CD3E', 'CD3G', 'CD4', 'CD8B', 'TRDC', 'CD8A', 'CD2',
                         'CD19', 'CD79B', 'JCHAIN', 'MZB1', 'IRF4',
                         'PCLAF', 'PCNA', 
                         'DNTT', 'RAG1', 'RAG2'),
            ncol = 6,
            reduction = 'umap')

#### Add normalized/scaled data to RNA assay
#dim(seu[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu <- NormalizeData(seu,  # normalize the RNA counts data per cell
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
seu <- ScaleData(seu, # scale the RNA counts data relative to other cells
                 assay = "RNA")
seu <- ScaleData(seu, # scale the SCT counts data relative to other cells
                 assay = "SCT")
#dim(seu[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#### Save the Seurat object
SaveH5Seurat(seu, 
             filename = "/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Thymus.h5Seurat", overwrite = TRUE)
rm(seu, seu.list.thy)

## Bone Marrow ----
#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list.bm, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list.bm <- PrepSCTIntegration(seu.list.bm, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list.bm, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features)
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT")
rm(seu.features, seu.anchors)

#### Calculate principle components
#Calculate 50 PCs:
seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(seu,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
# Tinker with distance and spread parameters to alter how UMAP looks. Example here: https://smorabit.github.io/blog/2020/umap/
seu <- RunUMAP(seu, 
               reduction = "pca",
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

seu <- RunTSNE(seu, 
               dims = PCdims,
               k.seed = 777,
               reduction = "pca") # create UMAP

#Visualize UMAP:
DimPlot(seu,
        reduction = 'umap',
        group.by = 'ID',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        split.by = 'ID',
        shuffle = TRUE,
        group.by = 'ID')

DimPlot(seu,
        reduction = 'umap',
        group.by = 'Chemistry',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'tissue',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'PigID',
        shuffle = TRUE)

FeaturePlot(seu, 
            features = c('nCount_SCT', 'nFeature_SCT', 'prcntMito'),
            reduction = 'umap')

DefaultAssay(seu) <- 'SCT'
FeaturePlot(seu, 
            features = c('AIF1', 'CSF1R', 'CD14',
                         'FLT3',
                         'CD3E', 'CD3G', 'CD4', 'CD8B', 'TRDC', 'CD8A', 'CD2',
                         'CD19', 'CD79B', 'JCHAIN', 'MZB1', 'IRF4',
                         'PCLAF', 'PCNA', 
                         'DNTT', 'RAG1', 'RAG2'),
            ncol = 6,
            reduction = 'umap')

#### Add normalized/scaled data to RNA assay
#dim(seu[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu <- NormalizeData(seu,  # normalize the RNA counts data per cell
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
seu <- ScaleData(seu, # scale the RNA counts data relative to other cells
                 assay = "RNA")
seu <- ScaleData(seu, # scale the SCT counts data relative to other cells
                 assay = "SCT")
#dim(seu[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#### Save the Seurat object
SaveH5Seurat(seu, 
             filename = "/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_BoneMarrow.h5Seurat", overwrite = TRUE)
rm(seu, seu.list.bm)

## Spleen ----
#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list.sp, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list.sp <- PrepSCTIntegration(seu.list.sp, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list.sp, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features)
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT")
rm(seu.features, seu.anchors)

#### Calculate principle components
#Calculate 50 PCs:
seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(seu,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
# Tinker with distance and spread parameters to alter how UMAP looks. Example here: https://smorabit.github.io/blog/2020/umap/
seu <- RunUMAP(seu, 
               reduction = "pca",
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

seu <- RunTSNE(seu, 
               dims = PCdims,
               k.seed = 777,
               reduction = "pca") # create UMAP

#Visualize UMAP:
DimPlot(seu,
        reduction = 'umap',
        group.by = 'ID',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        split.by = 'ID',
        shuffle = TRUE,
        group.by = 'ID')

DimPlot(seu,
        reduction = 'umap',
        group.by = 'Chemistry',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'tissue',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'PigID',
        shuffle = TRUE)

FeaturePlot(seu, 
            features = c('nCount_SCT', 'nFeature_SCT', 'prcntMito'),
            reduction = 'umap')

DefaultAssay(seu) <- 'SCT'
FeaturePlot(seu, 
            features = c('AIF1', 'CSF1R', 'CD14',
                         'FLT3',
                         'CD3E', 'CD3G', 'CD4', 'CD8B', 'TRDC', 'CD8A', 'CD2',
                         'CD19', 'CD79B', 'JCHAIN', 'MZB1', 'IRF4',
                         'PCLAF', 'PCNA', 
                         'DNTT', 'RAG1', 'RAG2'),
            ncol = 6,
            reduction = 'umap')

#### Add normalized/scaled data to RNA assay
#dim(seu[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu <- NormalizeData(seu,  # normalize the RNA counts data per cell
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
seu <- ScaleData(seu, # scale the RNA counts data relative to other cells
                 assay = "RNA")
seu <- ScaleData(seu, # scale the SCT counts data relative to other cells
                 assay = "SCT")
#dim(seu[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#### Save the Seurat object
SaveH5Seurat(seu, 
             filename = "/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen.h5Seurat", overwrite = TRUE)
rm(seu, seu.list.sp)

## Lymph Node ----
#Find variable features and integrate data:
seu.features <- SelectIntegrationFeatures(seu.list.ln, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list.ln <- PrepSCTIntegration(seu.list.ln, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list.ln, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features)
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT")
rm(seu.features, seu.anchors)

#### Calculate principle components
#Calculate 50 PCs:
seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(seu,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
rm(pct, cumu, co1, co2, pcs)

#### Perform multidimensional visualization of data:
#Create UMAP dimensions:
# Tinker with distance and spread parameters to alter how UMAP looks. Example here: https://smorabit.github.io/blog/2020/umap/
seu <- RunUMAP(seu, 
               reduction = "pca",
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

seu <- RunTSNE(seu, 
               dims = PCdims,
               k.seed = 777,
               reduction = "pca") # create UMAP

#Visualize UMAP:
DimPlot(seu,
        reduction = 'umap',
        group.by = 'ID',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        split.by = 'ID',
        shuffle = TRUE,
        group.by = 'ID')

DimPlot(seu,
        reduction = 'umap',
        group.by = 'Chemistry',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'tissue',
        shuffle = TRUE)

DimPlot(seu,
        reduction = 'umap',
        group.by = 'PigID',
        shuffle = TRUE)

FeaturePlot(seu, 
            features = c('nCount_SCT', 'nFeature_SCT', 'prcntMito'),
            reduction = 'umap')

DefaultAssay(seu) <- 'SCT'
FeaturePlot(seu, 
            features = c('AIF1', 'CSF1R', 'CD14',
                         'FLT3',
                         'CD3E', 'CD3G', 'CD4', 'CD8B', 'TRDC', 'CD8A', 'CD2',
                         'CD19', 'CD79B', 'JCHAIN', 'MZB1', 'IRF4',
                         'PCLAF', 'PCNA', 
                         'DNTT', 'RAG1', 'RAG2'),
            ncol = 6,
            reduction = 'umap')

#### Add normalized/scaled data to RNA assay
#dim(seu[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu <- NormalizeData(seu,  # normalize the RNA counts data per cell
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
seu <- ScaleData(seu, # scale the RNA counts data relative to other cells
                 assay = "RNA")
seu <- ScaleData(seu, # scale the SCT counts data relative to other cells
                 assay = "SCT")
#dim(seu[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(seu[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

#### Save the Seurat object
SaveH5Seurat(seu, 
             filename = "/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_LymphNode.h5Seurat", overwrite = TRUE)
rm(seu, seu.list.ln)

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scales_1.2.1                writexl_1.4.2               readxl_1.4.3                dplyr_1.1.3                 SeuratDisk_0.0.0.9020      
#[6] ggplot2_3.4.3               SeuratObject_4.1.3          Seurat_4.3.0.1              DropletUtils_1.16.0         SingleCellExperiment_1.18.1
#[11] SummarizedExperiment_1.26.1 Biobase_2.56.0              GenomicRanges_1.48.0        GenomeInfoDb_1.32.4         IRanges_2.30.1             
#[16] S4Vectors_0.34.0            BiocGenerics_0.42.0         MatrixGenerics_1.8.1        matrixStats_1.0.0          

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.3                spatstat.explore_3.2-3    reticulate_1.32.0         R.utils_2.12.2            tidyselect_1.2.0          htmlwidgets_1.6.2        
#[7] grid_4.2.2                BiocParallel_1.30.4       Rtsne_0.16                munsell_0.5.0             ScaledMatrix_1.4.1        codetools_0.2-19         
#[13] ica_1.0-3                 statmod_1.5.0             scran_1.24.1              xgboost_1.7.5.1           future_1.33.0             miniUI_0.1.1.1           
#[19] withr_2.5.0               spatstat.random_3.1-6     colorspace_2.1-0          progressr_0.14.0          rstudioapi_0.15.0         ROCR_1.0-11              
#[25] tensor_1.5                listenv_0.9.0             labeling_0.4.3            GenomeInfoDbData_1.2.8    polyclip_1.10-4           farver_2.1.1             
#[31] bit64_4.0.5               rhdf5_2.40.0              parallelly_1.36.0         vctrs_0.6.3               generics_0.1.3            R6_2.5.1                 
#[37] ggbeeswarm_0.7.2          rsvd_1.0.5                locfit_1.5-9.8            hdf5r_1.3.8               bitops_1.0-7              rhdf5filters_1.8.0       
#[43] spatstat.utils_3.0-3      DelayedArray_0.22.0       promises_1.2.1            BiocIO_1.6.0              beeswarm_0.4.0            gtable_0.3.4             
#[49] beachmat_2.12.0           Cairo_1.6-1               globals_0.16.2            goftest_1.2-3             rlang_1.1.1               splines_4.2.2            
#[55] rtracklayer_1.56.1        lazyeval_0.2.2            spatstat.geom_3.2-5       yaml_2.3.7                reshape2_1.4.4            abind_1.4-5              
#[61] httpuv_1.6.11             tools_4.2.2               ellipsis_0.3.2            RColorBrewer_1.1-3        ggridges_0.5.4            Rcpp_1.0.11              
#[67] plyr_1.8.8                sparseMatrixStats_1.8.0   zlibbioc_1.42.0           purrr_1.0.2               RCurl_1.98-1.12           deldir_1.0-9             
#[73] pbapply_1.7-2             viridis_0.6.4             cowplot_1.1.1             zoo_1.8-12                ggrepel_0.9.3             cluster_2.1.4            
#[79] magrittr_2.0.3            data.table_1.14.8         scattermore_1.2           lmtest_0.9-40             RANN_2.6.1                fitdistrplus_1.1-11      
#[85] patchwork_1.1.3           mime_0.12                 xtable_1.8-4              XML_3.99-0.14             gridExtra_2.3             compiler_4.2.2           
#[91] scater_1.24.0             tibble_3.2.1              KernSmooth_2.23-22        crayon_1.5.2              R.oo_1.25.0               htmltools_0.5.6          
#[97] later_1.3.1               tidyr_1.3.0               MASS_7.3-60               Matrix_1.6-1              cli_3.6.1                 R.methodsS3_1.8.2        
#[103] parallel_4.2.2            metapod_1.4.0             igraph_1.5.1              pkgconfig_2.0.3           GenomicAlignments_1.32.1  sp_2.0-0                 
#[109] plotly_4.10.2             scuttle_1.6.3             spatstat.sparse_3.0-2     vipor_0.4.5               dqrng_0.3.1               XVector_0.36.0           
#[115] stringr_1.5.0             digest_0.6.33             sctransform_0.3.5         RcppAnnoy_0.0.21          spatstat.data_3.0-1       Biostrings_2.64.1        
#[121] cellranger_1.1.0          leiden_0.4.3              uwot_0.1.16               edgeR_3.38.4              DelayedMatrixStats_1.18.2 restfulr_0.0.15          
#[127] shiny_1.7.5               Rsamtools_2.12.0          rjson_0.2.21              lifecycle_1.0.3           nlme_3.1-163              jsonlite_1.8.7           
#[133] Rhdf5lib_1.18.2           BiocNeighbors_1.14.0      viridisLite_0.4.2         limma_3.52.4              fansi_1.0.4               pillar_1.9.0             
#[139] lattice_0.21-8            ggrastr_1.0.2             fastmap_1.1.1             httr_1.4.7                survival_3.5-7            glue_1.6.2               
#[145] png_0.1-8                 bluster_1.6.0             bit_4.0.5                 stringi_1.7.12            HDF5Array_1.24.2          BiocSingular_1.12.0      
#[151] irlba_2.3.5.1             future.apply_1.11.0   
