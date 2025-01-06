library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(scales)

# Spleen (Madissoon et al.) ----

## Process reference data ----
# .rds downloaded from https://www.tissuestabilitycellatlas.org/; data in article at https://doi.org/10.1186/s13059-019-1906-x
ref <- readRDS('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceData/spleen_ts.rds')
Idents(ref) <- ref$Celltypes_GenomeBiol_2019
levels(ref) <- c('B_T_doublet', 'B_follicular', 'B_Hypermutation', 'B_mantle', 'CD34_progenitor',
                 'DC_1', 'DC_2', 'DC_activated', 'DC_plasmacytoid', 'ILC', 'Macrophage', 'Monocyte',
                 'NK_CD160pos', 'NK_dividing', 'NK_FCGR3Apos', 'Plasma_IgG', 'Plasma_IgM', 'Plasmablast', 'Platelet',
                 'T_CD4_conv', 'T_CD4_fh', 'T_CD4_naive', 'T_CD4_reg',
                 'T_CD8_activated', 'T_CD8_CTL', 'T_CD8_gd', 'T_CD8_MAIT', 'T_cell_dividing', 'Unknown') # Reorder the clusters based on phylo order from hierarchical tree
ref$Celltypes_GenomeBiol_2019 <- Idents(ref)
DimPlot(ref, cols = c('deepskyblue4', 'chocolate3', 'darkorange', 'chocolate2',
                      'goldenrod3', 'mediumpurple1', 'mediumpurple2', 'mediumpurple4', 'mediumorchid4',
                      'forestgreen', 'darkorchid3', 'darkorchid1',
                      'darkolivegreen3', 'darkolivegreen4', 'darkolivegreen1',
                      'magenta', 'maroon3', 'magenta3', 'gold',
                      'dodgerblue', 'dodgerblue1', 'dodgerblue2', 'dodgerblue3', 
                      'deepskyblue', 'deepskyblue1', 'deepskyblue2', 'deepskyblue3', 'steelblue1',
                      'grey60'))

# Convert to one-to-one ortholog genes only
# Load in ortho genes
orthoGenes <- read.delim("/home/Jayne.Wiarda/scRNAseqEpSI/PigToHuman_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Human.homology.type == 'ortholog_one2one') # subset to only one to one orthologs

# Filter to only one-to-one gene orthologs:
refgenes <- rownames(ref[['RNA']]@counts) # extract gene names from reference dataset
reforthos <- intersect(refgenes, orthoGenes$Human.gene.name) # find which gene names from reference are also one-to-one orthologs
#length(reforthos) # how many genes are orthologs?

# Create new Seurat object:
refcounts <- ref[['RNA']]@counts[rownames(ref[['RNA']]@counts) %in% reforthos,] # make count matrix from referemce, only taking counts from one-to-one ortholog genes
refMeta <- ref@meta.data # extract all the meta data from reference
ref <- CreateSeuratObject( # now create new Seurat object with only the filtered cells, orthologous genes, and meta data for filtered cells
  counts = refcounts, 
  meta.data = refMeta)

#Split object into list by sample ID:
ref.list <- SplitObject(ref, split.by = "sample") # split by sample IDs

#See how many cells/genes in each listed refrat object:
ref.list

## Normalize data using SCT method

#Performed SCTransform on each sample from list:
for (i in 1:length(ref.list)) { # normalize data using SCTransform method
  ref.list[[i]] <- SCTransform(ref.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

#Find variable features and integrate data:
ref.features <- SelectIntegrationFeatures(ref.list, # select the genes to use for integration
                                          verbose = TRUE) 
ref.list <- PrepSCTIntegration(ref.list, 
                               anchor.features = ref.features,
                               verbose = TRUE)
ref.anchors <- FindIntegrationAnchors(ref.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = ref.features)
ref <- IntegrateData(ref.anchors, # integrate data
                     normalization.method = "SCT")
rm(ref.features, ref.anchors)

#### Calculate principle components
#Calculate 50 PCs:
ref <- RunPCA(ref, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(ref,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- ref[["pca"]]@stdev / sum(ref[["pca"]]@stdev) * 100 # find standard deviation for each PC
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
ref <- RunUMAP(ref, 
               reduction = "pca",
               dims = 1:30,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

#Visualize UMAP:
DimPlot(ref,
        reduction = 'umap',
        group.by = 'sample',
        shuffle = TRUE)

DimPlot(ref,
        reduction = 'umap',
        group.by = 'Donor',
        shuffle = TRUE)

DimPlot(ref,
        reduction = 'umap',
        group.by = 'Time',
        shuffle = TRUE)

DimPlot(ref,
        reduction = 'umap',
        group.by = 'Celltypes_GenomeBiol_2019',
        shuffle = TRUE)

#### Add normalized/scaled data to RNA assay
#dim(ref[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
ref <- NormalizeData(ref,  # normalize the RNA counts data per cell
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
ref <- ScaleData(ref, # scale the RNA counts data relative to other cells
                 assay = "RNA")
ref <- ScaleData(ref, # scale the SCT counts data relative to other cells
                 assay = "SCT")
#dim(ref[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now
#dim(ref[["SCT"]]@scale.data) # see that all genes are scaled in SCT assay now

SaveH5Seurat(ref, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceData/HsSpleenReference_Madissoon.h5Seurat')

## Process query data ----

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen_CellTypeAnnotated_HBBhiRemoved.h5Seurat')

querygenes <- as.data.frame(rownames(seu[['RNA']]@counts)) # extract pig gene names from dataset
colnames(querygenes) <- 'gene'
pigGenes <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$Symbol97_Original <-gsub("_", "-", pigGenes$Symbol97_Original) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[pigGenes$Symbol97_Original %in% querygenes$gene, ] # slim down to only genes in our dataset
queryorthos <- intersect(pigGenes$ENSID, orthoGenes$Gene.stable.ID) # find which genes from reference are also one-to-one orthologs
length(queryorthos) # how many genes are orthologs?
pigGenes <- pigGenes[pigGenes$ENSID %in% queryorthos, ]
#dim(pigGenes)
orthoGenes <- orthoGenes[orthoGenes$Gene.stable.ID %in% pigGenes$ENSID, ] # slim down to only ortho genes in our dataset
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Human.gene.stable.ID, .keep_all = TRUE)  # retain only unique combinations of pig & human Ensembl IDs, ignoring transcript IDs
#dim(orthoGenes) # should not have same number of rows as in pigGenes
querycounts <- seu[['RNA']]@counts[rownames(seu[['RNA']]@counts) %in% pigGenes$Symbol97_Original,]
pigGenes <- pigGenes %>% arrange(factor(Symbol97_Original, levels = rownames(querycounts))) # arrange pigGenes in same order as querycounts
orthoGenes <- orthoGenes %>% arrange(factor(Gene.stable.ID, levels = pigGenes$ENSID)) # arrange orthoGenes in same order as pigGenes (and consequently querycounts)
rownames(querycounts) <- orthoGenes$Human.gene.name # change pig genes to human gene names
queryMeta <- seu@meta.data
seu <- CreateSeuratObject(counts = querycounts, # create new Seurat object for our query dataset
                          meta.data = queryMeta)
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceData/HumanizedSpleenQuery.h5Seurat')

#Split object into list by sample ID:
seu.list <- SplitObject(seu, split.by = "SampleID") # split by sample IDs

#See how many cells/genes in each listed seurat object:
seu.list

## Normalize data using SCT method

#Performed SCTransform on each sample from list:
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

for(i in 1:length(seu.list)) {
  seu.list[[i]] <- RunPCA(seu.list[[i]], npcs = 50) # have to calculate and store PCA to run MappingScore() function later
}

## Perform mapping/prediction ----

CellTypePredictions <- list()
MappingScores <- list()
for(i in 1:length(seu.list)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = seu.list[[i]],
    reduction = "cca", 
    dims = 1:30, 
    normalization.method = "SCT",
    recompute.residuals = FALSE) # have to set to FALSE to use SCT in Seurat v4 
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = ref$Celltypes_GenomeBiol_2019), 
                              dims = 1:30,
                              weight.reduction = "cca")
  MapScores <- MappingScore(
    anchors = anchors@anchors,
    combined.object = anchors@object.list[[1]],
    query.neighbors =  slot(object = seu.list[[i]], name = "neighbors")[["seu.list_ref.nn"]],
    query.weights = Tool(object = seu.list[[i]], slot = "TransferData")$weights.matrix,
    query.embeddings = Embeddings(object = seu.list[[i]]),
    ref.embeddings = Embeddings(object = ref),
    nn.method = "annoy",
    # n.trees = n.trees
  )
  CellTypePredictions[[i]] <- predictions
  MappingScores[[i]] <- MapScores
}


CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)
MappingScores <- Reduce(c,MappingScores)
MappingScores <- as.data.frame(MappingScores)

#Save the mapping & prediction results:

CellTypePredictions$CellBarcodes <- rownames(CellTypePredictions)
MappingScores$CellBarcodes <- rownames(MappingScores)
write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceMappingPrediction/HumanSpleenToPigSpleen_Madissoon_CellTypePredictions.xlsx')
write_xlsx(MappingScores, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceMappingPrediction/HumanSpleenToPigSpleen_Madissoon_MappingScores.xlsx')

#Incorporate cell prediction & mapping scores into original Seurat object of query data:
#CellTypePredictions <- read_xlsx('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceMappingPrediction/HumanSpleenToPigSpleen_Madissoon_CellTypePredictions.xlsx')
#MappingScores <- read_xlsx('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceMappingPrediction/HumanSpleenToPigSpleen_Madissoon_MappingScores.xlsx')
seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen_CellTypeAnnotated_HBBhiRemoved.h5Seurat')
DefaultAssay(seu) <- 'RNA'
seu <- AddMetaData(object = seu, 
                   metadata = c(CellTypePredictions, MappingScores))

#Visualize prediction & mapping results:

FeaturePlot(seu, reduction = 'tsne', features = 'MappingScores', min.cutoff = 0.4) & 
  scale_color_gradientn( colours = c('white', 'yellow', 'gold','orange', 'red')) & DarkTheme() # do we see high mapping scores, indicating good representation of query data by reference dataset? (Answer: yes)

VlnPlot(seu, split.by = 'CellType', features = 'MappingScores') # any major differences in mapping scores between ileum and jejunum data?
# Since jejunum and ileum samples have similar high mapping scores, suggests that the ileum atlas reference dataset annotations are cross-applicable to the jejunal samples used here

# Show example of a few prediction scores:
FeaturePlot(seu, # plot all prediction scores
            features = c('prediction.score.B_T_doublet', 'prediction.score.B_follicular', 'prediction.score.B_Hypermutation', 'prediction.score.B_mantle', 'prediction.score.CD34_progenitor',
                         'prediction.score.DC_1', 'prediction.score.DC_2', 'prediction.score.DC_activated', 'prediction.score.DC_plasmacytoid', 'prediction.score.ILC', 'prediction.score.Macrophage', 'prediction.score.Monocyte',
                         'prediction.score.NK_CD160pos', 'prediction.score.NK_dividing', 'prediction.score.NK_FCGR3Apos', 'prediction.score.Plasma_IgG', 'prediction.score.Plasma_IgM', 'prediction.score.Plasmablast', 'prediction.score.Platelet',
                         'prediction.score.T_CD4_conv', 'prediction.score.T_CD4_fh', 'prediction.score.T_CD4_naive', 'prediction.score.T_CD4_reg',
                         'prediction.score.T_CD8_activated', 'prediction.score.T_CD8_CTL', 'prediction.score.T_CD8_gd', 'prediction.score.T_CD8_MAIT', 'prediction.score.T_cell_dividing', 'prediction.score.Unknown'),
            reduction = 'tsne', # change to 'tsne' to overlay onto t-SNE plot
            ncol = 7, cols = c('grey75', 'red3')) & NoAxes()

VlnPlot(seu, group.by = 'predicted.id', features = 'prediction.score.max') # any major differences in prediction scores between ileum and jejunum data? 
# See mostly high prediction scores, giving confidence to cell type identities assigned via reference annotation-based prediction
DimPlot(seu, group.by = 'predicted.id', reduction = 'tsne', label = FALSE, cols = c('deepskyblue4', 'chocolate3', 'darkorange', 'chocolate2',
                                                                                    'goldenrod3', 'mediumpurple2', 'mediumorchid4',
                                                                                    'forestgreen', 'darkorchid3', 'darkorchid1',
                                                                                    'darkolivegreen3', 'darkolivegreen4', 'darkolivegreen1',
                                                                                    'magenta', 'maroon3', 'magenta3',
                                                                                    'dodgerblue', 'dodgerblue1', 'dodgerblue2', 'dodgerblue3', 
                                                                                    'deepskyblue', 'deepskyblue1', 'deepskyblue2', 'deepskyblue3', 'steelblue1')) # plot predicted cell IDs (based on highest prediction scores)

# Dot Plot annotation comparison ----

df <- data.frame(seu$CellType, seu$predicted.id, seu$MappingScores)
colnames(df) <- c('Annotation', 'PredictedID', 'MappingScore')

count <- data.frame(table(df$Annotation, df$PredictedID))
colnames(count) <- c('Annotation', 'PredictedID', 'Count')

pct_summary <- count %>%
  group_by(Annotation, PredictedID) %>%
  summarise(Count = sum(Count)) %>%
  mutate(Percent = prop.table(Count) * 100) %>%
  ungroup


map_summary <- df %>%
  group_by(Annotation, PredictedID) %>%
  summarise(MappingScore = mean(MappingScore))

pct_summary$combo <- paste(pct_summary$Annotation, pct_summary$PredictedID, sep = '_')
map_summary$combo <- paste(map_summary$Annotation, map_summary$PredictedID, sep = '_')
newdf <- merge(pct_summary, map_summary, by = "combo", all = T)

newdf$Annotation.x <- factor(newdf$Annotation.x, levels=c('Bcell', 'cDC', 'Mac/Mono',
                                                          'ILC_NCR1posEOMESpos', 'ILC_Cytotoxic',
                                                          'ASC', 'gdT_CD2neg', 'abT_CCL5neg', 
                                                          'abT_Cytotoxic', 'abT_CCL5pos', 
                                                          'gdT_CD2pos', 'T/ILC_Cycling'))
newdf$PredictedID.x <- factor(newdf$PredictedID.x, levels=c('B_T_doublet', 'B_follicular', 'B_Hypermutation', 'B_mantle', 'CD34_progenitor',
                 'DC_1', 'DC_2', 'DC_activated', 'DC_plasmacytoid', 'ILC', 'Macrophage', 'Monocyte',
                 'NK_CD160pos', 'NK_dividing', 'NK_FCGR3Apos', 'Plasma_IgG', 'Plasma_IgM', 'Plasmablast', 'Platelet',
                 'T_CD4_conv', 'T_CD4_fh', 'T_CD4_naive', 'T_CD4_reg',
                 'T_CD8_activated', 'T_CD8_CTL', 'T_CD8_gd', 'T_CD8_MAIT', 'T_cell_dividing', 'Unknown'))

ggplot(newdf, aes(y=Annotation.x, x=PredictedID.x)) +
  geom_point(aes(size = Percent, fill = MappingScore), color="black", shape=21) +
  scale_size("Percentage of Annotated Cells", range = c(0,6)) +
  scale_fill_gradient(high = 'navy', low = 'yellow',
                      name = "Average Mapping Score",
                      limits=c(0.5,1), oob = squish) +
  ylab("Annotation") + xlab("Predicted ID") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14))

# Merge pig & human ILCs ----
que <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceData/HumanizedSpleenQuery.h5Seurat')
Idents(que) <- seu$CellType
subpor <- subset(que, idents = c('ILC_NCR1posEOMESpos', 'ILC_Cytotoxic'))
subpor$sample <- subpor$SampleID
Idents(ref) <- ref$Celltypes_GenomeBiol_2019
subhum <- subset(ref, idents = c('NK_FCGR3Apos', 'NK_CD160pos'))
subhum$CellType <- subhum$Celltypes_GenomeBiol_2019
ILC <- merge(subpor, subhum)
DefaultAssay(ILC) <- 'SCT'

#Find variable features and integrate data:
ILC.list <- SplitObject(ILC, split.by = "sample") # split by sample IDs
for (i in 1:length(ILC.list)) { # normalize data using SCTransform method
  ILC.list[[i]] <- SCTransform(ILC.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
ILC.features <- SelectIntegrationFeatures(ILC.list, # select the genes to use for integration
                                          verbose = TRUE) 
ILC.list <- PrepSCTIntegration(ILC.list, 
                               anchor.features = ILC.features,
                               verbose = TRUE)
ILC.anchors <- FindIntegrationAnchors(ILC.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = ILC.features)
ILC <- IntegrateData(ILC.anchors, # integrate data
                     normalization.method = "SCT")
rm(ILC.features, ILC.anchors)

#### Calculate principle components
#Calculate 50 PCs:
ILC <- RunPCA(ILC, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(ILC,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- ILC[["pca"]]@stdev / sum(ILC[["pca"]]@stdev) * 100 # find standard deviation for each PC
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
ILC <- RunUMAP(ILC, 
               reduction = "pca",
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

#Visualize UMAP:
Idents(ILC) <- ILC$CellType
levels(ILC) <- c('ILC_NCR1posEOMESpos', 'NK_CD160pos',
                 'ILC_Cytotoxic', 'NK_FCGR3Apos')
ILC$CellType <- Idents(ILC)

DimPlot(ILC,
        reduction = 'umap',
        group.by = 'sample',
        shuffle = TRUE, cols = c(rep('cornflowerblue', 19), rep('red3', 2)))

DimPlot(ILC,
        reduction = 'umap',
        group.by = 'CellType', split.by = 'CellType',
        shuffle = TRUE, ncol = 4, cols = c(rep('black', 4)))

ILC <- BuildClusterTree(ILC, 
                        dims = PCdims, 
                        assay = "PCA")
PlotClusterTree(ILC, 
                edge.width = 3) # plot tree with node labels
data.tree <- Tool(object = ILC, 
                  slot = "BuildClusterTree") 
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)

SaveH5Seurat(ILC, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceData/MergedHumanPigILCs.h5Seurat')

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.5 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scales_1.3.0          writexl_1.4.2         readxl_1.4.3          dplyr_1.1.4           SeuratDisk_0.0.0.9020 ggplot2_3.5.0         SeuratObject_4.1.3    Seurat_4.3.0.1       

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.16             colorspace_2.1-0       rjson_0.2.21           deldir_1.0-9           ellipsis_0.3.2         ggridges_0.5.4         circlize_0.4.16        GlobalOptions_0.1.2    clue_0.3-64            rstudioapi_0.15.0      spatstat.data_3.0-1    leiden_0.4.3           listenv_0.9.1          ggrepel_0.9.5          bit64_4.0.5           
#[16] fansi_1.0.6            codetools_0.2-19       splines_4.2.2          doParallel_1.0.17      polyclip_1.10-4        jsonlite_1.8.8         ica_1.0-3              cluster_2.1.4          png_0.1-8              uwot_0.1.16            shiny_1.8.0            sctransform_0.3.5      spatstat.sparse_3.0-2  compiler_4.2.2         httr_1.4.7            
#[31] Matrix_1.6-1           fastmap_1.1.1          lazyeval_0.2.2         cli_3.6.2              later_1.3.2            formatR_1.14           htmltools_0.5.7        tools_4.2.2            igraph_1.5.1           gtable_0.3.4           glue_1.7.0             RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.12            scattermore_1.2       
#[46] cellranger_1.1.0       vctrs_0.6.5            ape_5.7-1              spatstat.explore_3.2-3 nlme_3.1-163           progressr_0.14.0       iterators_1.0.14       lmtest_0.9-40          spatstat.random_3.1-6  stringr_1.5.1          globals_0.16.2         mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1         
#[61] goftest_1.2-3          future_1.33.1          MASS_7.3-60            zoo_1.8-12             promises_1.2.1         spatstat.utils_3.0-3   parallel_4.2.2         lambda.r_1.2.4         RColorBrewer_1.1-3     ComplexHeatmap_2.12.1  reticulate_1.35.0      pbapply_1.7-2          gridExtra_2.3          stringi_1.8.3          S4Vectors_0.34.0      
#[76] foreach_1.5.2          BiocGenerics_0.42.0    shape_1.4.6.1          rlang_1.1.3            pkgconfig_2.0.3        matrixStats_1.0.0      lattice_0.21-8         ROCR_1.0-11            purrr_1.0.2            tensor_1.5             patchwork_1.2.0        htmlwidgets_1.6.4      cowplot_1.1.3          bit_4.0.5              tidyselect_1.2.0      
#[91] parallelly_1.37.1      RcppAnnoy_0.0.21       plyr_1.8.9             magrittr_2.0.3         R6_2.5.1               IRanges_2.30.1         generics_0.1.3         pillar_1.9.0           withr_3.0.0            fitdistrplus_1.1-11    survival_3.5-7         abind_1.4-5            sp_2.0-0               tibble_3.2.1           future.apply_1.11.1   
#[106] crayon_1.5.2           hdf5r_1.3.8            futile.options_1.0.1   KernSmooth_2.23-22     utf8_1.2.4             spatstat.geom_3.2-5    plotly_4.10.4          GetoptLong_1.0.5       data.table_1.15.2      digest_0.6.34          xtable_1.8-4           tidyr_1.3.1            httpuv_1.6.14          stats4_4.2.2           munsell_0.5.0         
#[121] viridisLite_0.4.2  
