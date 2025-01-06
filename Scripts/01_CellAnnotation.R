library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)
library(cowplot)

cols = c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown',
         'cornflowerblue', 'navy', 'lightpink', 'salmon', 
         'deepskyblue2', 'tan4', 'mediumpurple1', 
         'darkgreen', 'gray50', 'darkmagenta', 'red', 'hotpink', 'khaki', 
         'orange3', 'limegreen', 'cadetblue3', 'firebrick', 'deepskyblue4', 'darkseagreen', 'burlywood3', 'black',
         'goldenrod3', 'blue', 'deeppink', 'lightskyblue3', 'mistyrose3', 'purple')

# Spleen ----

## Load data ----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen.h5Seurat')

## Define clusters ----
DefaultAssay(seu) <- 'integrated'
seu <- FindClusters(seu, 
                    resolution = 1.5) # toggled with this parameter in increments of 0.5 to find lowest resolution still suited to defining general cell types
DimPlot(seu, label = TRUE, reduction = 'tsne')

## Assess canonical gene expression ----
DefaultAssay(seu) <- 'RNA'
DotPlot(seu, features = unique(c('CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRDC', 'CD2',
                                 'NCR1', 'ID3', 'ID2', 'IL7R', 'CTLA4', 'ICOS', 'CD40LG',
                                 'PCLAF', 'LTB', 'CD69', 'PDCD1', 'KLF2', 'LEF1', 'S1PR1', 'CCR7', 'CXCR4', 'CXCR5', 'BCL6',
                                 'KLRK1', 'GNLY', 'GZMB', 'ITGAE', 'CCL5', 'RORC', 'KLRB1', 'ITGB7', 'CCR9',
                                 'IFNG', 'TGFB1', 'TNF', 'IL17F', 'IL10', 'IL2', 'CD79A', 'CD79B', 'CD19', 'PAX5', 'MZB1',
                                 'PRDM1', 'IRF4', 'JCHAIN', 'AICDA','XBP1',  'CD83', 'SLA-DQB1', 'HLA-DRA',
                                 'GPR183','KLF2', 'CD40', 'SELL', 'CD86',  'AIF1', 'CST3', 'CD163', 'CSF1R', 'ENSSSCG00000036618',
                                 'CD14', 'FLT3',  'AIRE', 'C1QB', 'MS4A2', 'NLRP3', 'CD68', 'TLR4',
                                 'CXCL14', 'CCR2', 'FCER1A',
                                 'TLR4', 'NLRP3',  'CD80', 'LMNA',
                                 'KIT', 'DNTT', 'RAG1', 'RAG2', 'HBB'))) + RotatedAxis()

## Identify DEGs in clusters ----
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_Clusters_Spleen.xlsx')

## Assign cell annotations ----
Idents(seu) <- seu$seurat_clusters
seu$CellType <- seu$seurat_clusters
Idents(seu) <- seu$CellType
types <- c('gdT_CD2neg', 'ILC_Cytotoxic', 'HBBhi', 'gdT_CD2neg', 'abT_CCL5neg', 'gdT_CD2neg', 
           'Bcell', 'ILC_NCR1posEOMESpos', 'abT_Cytotoxic', 'ILC_Cytotoxic', 'abT_CCL5pos',
           'gdT_CD2neg', 'ILC_Cytotoxic', 'T/ILC_Cycling', 'Mac/Mono', 'abT_CCL5neg',
           'ASC', 'abT_Cytotoxic', 'gdT_CD2pos', 'cDC') 
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$CellType <- Idents(seu)
Idents(seu) <- seu$CellType

## Visualize data ----
#DimPlot(seu, reduction = 'tsne', label = TRUE, cols = cols[1:(length(unique(seu$CellType)))])
DimPlot(seu, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(seu$CellType)))], group.by = 'CellType')
DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'seurat_clusters')

## Perform hierarchical clustering ----
#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our seusequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for seusequent steps
rm(pct, cumu, co1, co2, pcs)

Idents(seu) <- seu$CellType
seu <- BuildClusterTree(seu, 
                        dims = PCdims)
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('ASC', 'Mac/Mono', 'cDC', 'Bcell',
                                            'T/ILC_Cycling', 'ILC_NCR1posEOMESpos', 'ILC_Cytotoxic',
                                            'gdT_CD2pos', 'gdT_CD2neg', 'abT_CCL5neg', 'abT_Cytotoxic', 'abT_CCL5pos',
                                            'HBBhi'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)
Idents(seu) <- seu$CellType
levels(seu) <- rev(c('ASC', 'Mac/Mono', 'cDC', 'Bcell',
                     'T/ILC_Cycling', 'ILC_NCR1posEOMESpos', 'ILC_Cytotoxic',
                     'gdT_CD2pos', 'gdT_CD2neg', 'abT_CCL5neg', 'abT_Cytotoxic', 'abT_CCL5pos',
                     'HBBhi')) # Reorder the clusters based on phylo order from hierarchical tree

## Remove HBBhi cluster ----
seu$CellType <- Idents(seu)
DotPlot(seu,
        features = unique(c('CD3E', 'CD3G', 'CD2', 'TRDC', 'CD4', 'CD8A', 'CD8B', 'KLRB1', 'NCR1', 'EOMES', 
                            'GZMB', 'PRF1', 'GNLY', 
                            'CCL5', 'PCLAF', 'ENSSSCG00000026302', #MKI67
                            'CD79B', 'PAX5', 'CD19', 'JCHAIN', 'XBP1', 'PRDM1',
                            'AIF1', 'CST3', 'CSF1R', 'FLT3', 'HBB', 'RPS6')),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()

SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen_CellTypeAnnotated_HBBhiIncluded.h5Seurat')

Idents(seu) <- seu$CellType
seu <- subset(seu, idents = 'HBBhi', invert = TRUE)
DefaultAssay(seu) <- 'RNA'

## Identify DEGs in cell types ----
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_CellTypes_Spleen.xlsx')

## Redo processing after HBBhi cell removal ----
#Find variable features and integrate data:
DefaultAssay(seu) <- 'RNA'
seu.list <- SplitObject(seu, split.by = "ID") # split by sample IDs
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
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
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

seu <- RunTSNE(seu, 
               dims = PCdims,
               k.seed = 777,
               reduction = "pca") # create UMAP

Idents(seu) <- seu$CellType
seu <- BuildClusterTree(seu, 
                        dims = PCdims)
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('ASC', 'Mac/Mono', 'cDC', 'Bcell',
                                            'T/ILC_Cycling', 'ILC_NCR1posEOMESpos', 'ILC_Cytotoxic',
                                            'gdT_CD2pos', 'gdT_CD2neg', 'abT_CCL5neg', 'abT_Cytotoxic', 'abT_CCL5pos'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)
Idents(seu) <- seu$CellType
levels(seu) <- rev(c('ASC', 'Mac/Mono', 'cDC', 'Bcell',
                     'T/ILC_Cycling', 'ILC_NCR1posEOMESpos', 'ILC_Cytotoxic',
                     'gdT_CD2pos', 'gdT_CD2neg', 'abT_CCL5neg', 'abT_Cytotoxic', 'abT_CCL5pos')) # Reorder the clusters based on phylo order from hierarchical tree
seu$CellType <- Idents(seu)
DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = unique(c('CD3E', 'CD3G', 'CD2', 'TRDC', 'CD4', 'CD8A', 'CD8B', 'KLRB1', 'NCR1', 'EOMES', 
                            'GZMB', 'PRF1', 'GNLY', 
                            'CCL5', 'PCLAF', 'ENSSSCG00000026302', #MKI67
                            'CD79B', 'PAX5', 'CD19', 'JCHAIN', 'XBP1', 'PRDM1',
                            'AIF1', 'CST3', 'CSF1R', 'FLT3', 'HBB', 'RPS6')),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()

## Visualize data ----
#DimPlot(seu, reduction = 'tsne', label = TRUE, cols = cols[1:(length(unique(seu$CellType)))])
DimPlot(seu, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(seu$CellType)))], group.by = 'CellType')
DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'seurat_clusters')

## Save final dataset ----
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen_CellTypeAnnotated_HBBhiRemoved.h5Seurat')

# Lymph Node ----

## Load data ----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_LymphNode.h5Seurat')

## Define clusters ----
DefaultAssay(seu) <- 'integrated'
seu <- FindClusters(seu, 
                    resolution = 4.5) # toggled with this parameter in increments of 0.5 to find lowest resolution still suited to defining general cell types
DimPlot(seu, label = TRUE, reduction = 'tsne')

## Assess canonical gene expression ----
DefaultAssay(seu) <- 'RNA'
DotPlot(seu, features = unique(c('CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRDC', 'CD2',
                                 'NCR1', 'ID3', 'ID2', 'IL7R', 'CTLA4', 'ICOS', 'CD40LG',
                                 'PCLAF', 'LTB', 'CD69', 'PDCD1', 'KLF2', 'LEF1', 'S1PR1', 'CCR7', 'CXCR4', 'CXCR5', 'BCL6',
                                 'KLRK1', 'GNLY', 'GZMB', 'ITGAE', 'CCL5', 'RORC', 'KLRB1', 'ITGB7', 'CCR9',
                                 'IFNG', 'TGFB1', 'TNF', 'IL17F', 'IL10', 'IL2', 'CD79A', 'CD79B', 'CD19', 'PAX5', 'MZB1',
                                 'PRDM1', 'IRF4', 'JCHAIN', 'AICDA','XBP1',  'CD83', 'SLA-DQB1', 'HLA-DRA',
                                 'GPR183','KLF2', 'CD40', 'SELL', 'CD86',  'AIF1', 'CST3', 'CD163', 'CSF1R', 'ENSSSCG00000036618',
                                 'CD14', 'FLT3',  'CD3E', 'CD3G', 'CD2', 'TRDC', 'CD8A', 'CD8B', 'CD4', 'KLRB1', 
                                 'CD79B', 'JCHAIN', 'AICDA', 'PCLAF', 'BCL6',
                                 'CD40LG', 'CD40', 'IL7R', 'CTLA4', 'ICOS', 'CD40LG', 'LTB', 'CD69', 'PDCD1', 'KLF2', 'LEF1', 'S1PR1', 'CXCR4', 'CXCR5', 
                                 'KLRK1', 'GNLY', 'GZMB', 'ITGAE', 'CCL5', 'CCR9', 'CCR7',
                                 'AIF1', 'CSF1R', 'FLT3'))) + RotatedAxis()

## Identify DEGs in clusters ----
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_Clusters_LymphNode.xlsx')

## Assign cell annotations ----
Idents(seu) <- seu$seurat_clusters
seu$CellType <- seu$seurat_clusters
Idents(seu) <- seu$CellType
types <- c('B_AICDAneg', 'CD8T_KLF2pos', 'CD4T_KLF2pos', 'B_AICDAneg', 'CD4T_KLF2neg', 'B_AICDApos',
           'B_AICDAneg', 'CD8T_KLF2pos', 'B_AICDAneg', 'B_AICDAneg', 'CD4T_KLF2pos',
           'B_AICDAneg', 'CD4T_KLF2pos', 'B_AICDAneg', 'CD4T_KLF2pos', 'B_AICDApos',
           'B_AICDAneg', 'B_AICDApos_Cycling', 'CD4T_KLF2pos', 'B_AICDApos', 'B_AICDAneg',
           'CD8T_KLF2neg', 'CD4T_KLF2pos', 'B_AICDAneg', 'B_AICDAneg', 'B_AICDAneg',
           'CD4Tfh', 'CD8T_KLF2pos', 'CD4T_KLF2pos', 'B_AICDAneg', 'B_AICDAneg',
           'B_AICDAneg', 'CD4T_KLF2pos', 'CD8T_KLF2pos', 'CD4T_KLF2pos', 'B_AICDAneg',
           'B_AICDApos_Cycling', 'B_AICDAneg', 'B_AICDAneg', 'B_AICDApos_Cycling', 'B_AICDAneg',
           'B_AICDAneg', 'CD8T_KLF2pos', 'B_AICDAneg', 'B_AICDAneg_Cycling', 'B_AICDApos_Cycling',
           'B_AICDApos_Cycling', 'CD8T_KLF2pos', 'gdT_CD2neg', 'CD4T_KLF2pos', 'B_AICDAneg_Cycling',
           'B_AICDAneg', 'B_AICDAneg_Cycling', 'TILC_CCL5pos', 'ILC_KITpos', 'gdT_CD2pos',
           'Mac/Mono/cDC', 'ASC') 
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$CellType <- Idents(seu)
Idents(seu) <- seu$CellType

## Visualize data ----
#DimPlot(seu, reduction = 'tsne', label = TRUE, cols = cols[1:(length(unique(seu$CellType)))])
DimPlot(seu, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(seu$CellType)))], group.by = 'CellType')
DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'seurat_clusters')

## Perform hierarchical clustering ----
#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our seusequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for seusequent steps
rm(pct, cumu, co1, co2, pcs)

Idents(seu) <- seu$CellType
seu <- BuildClusterTree(seu, 
                         dims = PCdims)
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                                            'Mac/Mono/cDC',
                                            'TILC_CCL5pos', 'ILC_KITpos',
                                            'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                                            'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)
Idents(seu) <- seu$CellType
levels(seu) <- rev(c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                 'Mac/Mono/cDC',
                 'TILC_CCL5pos', 'ILC_KITpos',
                 'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                 'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg')) # Reorder the clusters based on phylo order from hierarchical tree
seu$CellType <- Idents(seu)

DotPlot(seu,
        features = unique(c('CD3E', 'CD3G', 'CD2', 'TRDC', 'CD4', 'CD8A', 'CD8B', 'KLRB1', 'RORC', 
                     'CCR7', 'S1PR1', 'KLF2', 'CD69', 'GZMB', 'CCR9',
                     'ICOS', 'CTLA4', 'PDCD1',
                     'CCL5', 'KLRK1', 'IL7R',
                     'CD79B', 'PAX5', 'CD19', 'JCHAIN', 'XBP1', 'PRDM1',
                     'CXCR4', 'CCR7', 'CD69', 'CD40', 'S1PR1', 'KLF2',
                     'AICDA', 'BCL6', 'PCLAF', 'ENSSSCG00000026302', #MKI67
                     'AIF1', 'CST3', 'CSF1R', 'FLT3', 'RPS6')),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()

## Identify DEGs in cell types ----
Idents(seu) <- seu$CellType
DefaultAssay(seu) <- 'RNA'
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_CellTypes_LymphNode.xlsx')

## Save final dataset ----
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_LymphNode_CellTypeAnnotated.h5Seurat')

# Thymus ----

## Load data ----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Thymus.h5Seurat')

## Define clusters ----
DefaultAssay(seu) <- 'integrated'
seu <- FindClusters(seu, 
                    resolution = 4.5) # toggled with this parameter in increments of 0.5 to find lowest resolution still suited to defining general cell types
DimPlot(seu, label = TRUE, reduction = 'tsne')

## Assess canonical gene expression ----
DefaultAssay(seu) <- 'RNA'
DotPlot(seu, features = unique(c('CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRDC', 'CD2',
                                 'NCR1', 'ID3', 'ID2', 'IL7R', 'CTLA4', 'ICOS', 'CD40LG',
                                 'PCLAF', 'LTB', 'CD69', 'PDCD1', 'KLF2', 'LEF1', 'S1PR1', 'CCR7', 'CXCR4', 'CXCR5', 'BCL6',
                                 'KLRK1', 'GNLY', 'GZMB', 'ITGAE', 'CCL5', 'FOXP3', 'RORC', 'KLRB1', 'ITGB7', 'CCR9',
                                 'IFNG', 'TGFB1', 'TNF', 'IL17F', 'IL10', 'IL2', 'CD79A', 'CD79B', 'CD19', 'PAX5', 'MZB1',
                                 'PRDM1', 'IRF4', 'JCHAIN', 'AICDA','XBP1',  'CD83', 'SLA-DQB1', 'HLA-DRA',
                                 'GPR183','KLF2', 'CD40', 'SELL', 'CD86',  'AIF1', 'CST3', 'CD163', 'CSF1R', 'ENSSSCG00000036618',
                                 'CD14', 'FLT3',  'AIRE', 'C1QB', 'MS4A2', 'NLRP3', 'CD68', 'TLR4',
                                 'CXCL14', 'CCR2', 'FCER1A',
                                 'TLR4', 'NLRP3',  'CD80', 'LMNA',
                                 'KIT', 'DNTT', 'RAG1', 'RAG2', 'HBB'))) + RotatedAxis()

## Identify DEGs in clusters ----
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_Clusters_Thymus.xlsx')

## Assign cell annotations ----
Idents(seu) <- seu$seurat_clusters
seu$CellType <- seu$seurat_clusters
Idents(seu) <- seu$CellType
types <- c('ThymocyteDP', 'CD4T_KLF2pos', 'ThymocyteDN_Cycling', 'ThymocyteDN_Cycling', 'gdT_CD2neg', 'gdT_CD2posCD8Aneg', 
           'CD4T_KLF2neg_CCR9neg', 'CD8T_KLF2neg', 'Bcell', 'CD4T_KLF2neg_CCR9neg', 'CD4T_KLF2neg_CCR9pos',
           'CD8T_KLF2pos', 'ThymocyteDN', 'ThymocyteDP_Cycling', 'ThymocyteDP_Cycling', 'ThymocyteDP_Cycling',
           'CD4T_KLF2neg_CCR9pos', 'gdT_CD2pos_Cycling', 'ThymocyteDP_Cycling', 'gdT_CD2neg', 'ThymocyteDN_Cycling',
           'gdT_CD2neg', 'ThymocyteDP_Cycling', 'CD8T_KLF2pos', 'Bcell', 'ThymocyteDP_Cycling',
           'ThymocyteDP_Cycling', 'CD4T_KLF2neg_CCR9pos', 'ThymocyteDP', 'abT_CD1DposCD1Epos', 'CD8T_KLF2pos',
           'ThymocyteDN_Cycling', 'CD8T_KLF2neg', 'gdT_CD2posCD8Apos', 'gdT_CD2neg', 'ThymocyteDP_Cycling',
           'ThymocyteDN', 'TILC_CCL5pos', 'Bcell', 'ThymocyteDN_Cycling', 'ThymocyteDP_Cycling', 
           'TILC_CCL5pos', 'gdT_CD2pos_Cycling', 'ThymocyteDN_Cycling', 'Progenitor', 'gdT_CD2neg_Cycling',
           'CD4T_KLF2neg_CCR9pos', 'Mac/Mono/cDC', 'Bcell', 'abT_Cycling', 'ASC', 'pDC') 
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$CellType <- Idents(seu)
Idents(seu) <- seu$CellType

## Visualize data ----
#DimPlot(seu, reduction = 'tsne', label = TRUE, cols = cols[1:(length(unique(seu$CellType)))])
DimPlot(seu, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(seu$CellType)))], group.by = 'CellType')
DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'seurat_clusters')

## Perform hierarchical clustering ----
#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our seusequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for seusequent steps
rm(pct, cumu, co1, co2, pcs)

Idents(seu) <- seu$CellType
seu <- BuildClusterTree(seu, 
                        dims = PCdims)
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('pDC', 'Mac/Mono/cDC', 'ASC', 'Bcell', 
                                            'abT_Cycling', 'gdT_CD2neg_Cycling', 
                                            'gdT_CD2neg', 'gdT_CD2posCD8Aneg', 'CD4T_KLF2pos', 'CD8T_KLF2pos', 'CD4T_KLF2neg_CCR9neg',
                                            'CD4T_KLF2neg_CCR9pos', 'CD8T_KLF2neg', 'abT_CD1DposCD1Epos',
                                            'gdT_CD2pos_Cycling', 'ThymocyteDN_Cycling', 'ThymocyteDP_Cycling', 'ThymocyteDN', 'ThymocyteDP',
                                            'Progenitor', 'gdT_CD2posCD8Apos', 'TILC_CCL5pos'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)
Idents(seu) <- seu$CellType
levels(seu) <- rev(c('pDC', 'Mac/Mono/cDC', 'ASC', 'Bcell', 
                     'abT_Cycling', 'gdT_CD2neg_Cycling', 
                     'gdT_CD2neg', 'gdT_CD2posCD8Aneg', 'CD4T_KLF2pos', 'CD8T_KLF2pos', 'CD4T_KLF2neg_CCR9neg',
                     'CD4T_KLF2neg_CCR9pos', 'CD8T_KLF2neg', 'abT_CD1DposCD1Epos',
                     'gdT_CD2pos_Cycling', 'ThymocyteDN_Cycling', 'ThymocyteDP_Cycling', 'ThymocyteDN', 'ThymocyteDP',
                     'Progenitor', 'gdT_CD2posCD8Apos', 'TILC_CCL5pos')) # Reorder the clusters based on phylo order from hierarchical tree
seu$CellType <- Idents(seu)
DotPlot(seu,
        features = unique(c('CD3E', 'CD3G', 'CD2', 'TRDC', 'CD4', 'CD8A', 'CD8B', 'CD1D', 'CD1E',
                            'KLRB1', 'KLRK1', 
                            'CCL5', 'KLF2', 'S1PR1', 'IL7R', 'IL2RA',
                            'PCLAF', 'ENSSSCG00000026302', #MKI67
                            'DNTT', 'RAG1', 'RAG2', 'KIT', 'RORC', 
                            'CD79B', 'PAX5', 'CD19', 'JCHAIN', 'XBP1', 'PRDM1',
                            'AIF1', 'CST3', 'CSF1R', 'FLT3', 
                            'CD93', 'TCF4', 'CLEC12A', 'IRF8', 'HBB', 'RPS6')),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()

## Identify DEGs in cell types ----
Idents(seu) <- seu$CellType
DefaultAssay(seu) <- 'RNA'
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_CellTypes_Thymus.xlsx')

## Save final dataset ----
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Thymus_CellTypeAnnotated.h5Seurat')

# Bone Marrow ----

## Load data ----
seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_BoneMarrow.h5Seurat')

## Define clusters ----
DefaultAssay(seu) <- 'integrated'
seu <- FindClusters(seu, 
                    resolution = 4.5) # toggled with this parameter in increments of 0.5 to find lowest resolution still suited to defining general cell types
DimPlot(seu, label = TRUE, reduction = 'tsne')

## Assess canonical gene expression ----
DefaultAssay(seu) <- 'RNA'
DotPlot(seu, features = unique(c('CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRDC', 'CD2',
                                 'NCR1', 'ID3', 'ID2', 'IL7R', 'CTLA4', 'ICOS', 'CD40LG',
                                 'PCLAF', 'LTB', 'CD69', 'PDCD1', 'KLF2', 'LEF1', 'S1PR1', 'CCR7', 'CXCR4', 'CXCR5', 'BCL6',
                                 'KLRK1', 'GNLY', 'GZMB', 'ITGAE', 'CCL5', 'RORC', 'KLRB1', 'ITGB7', 'CCR9',
                                 'IFNG', 'TGFB1', 'TNF', 'IL17F', 'IL10', 'IL2', 'CD79A', 'CD79B', 'CD19', 'PAX5', 'MZB1',
                                 'PRDM1', 'IRF4', 'JCHAIN', 'AICDA','XBP1',  'CD83', 'SLA-DQB1', 'HLA-DRA',
                                 'GPR183','KLF2', 'CD40', 'SELL', 'CD86',  'AIF1', 'CST3', 'CD163', 'CSF1R', 'ENSSSCG00000036618',
                                 'CD14', 'FLT3',  'AIRE', 'C1QB', 'MS4A2', 'NLRP3', 'CD68', 'TLR4',
                                 'CXCL14', 'CCR2', 'FCER1A',
                                 'TLR4', 'NLRP3',  'CD80', 'LMNA',
                                 'KIT', 'DNTT', 'RAG1', 'RAG2', 'HBB',
                                 'CKB', 'APOE', 'ATP6V0D2', 'NFATC1', 'MMP9', 'CALCR', 'SIGLEC15', 'ACP5', 'DCSTAMP', 'OCSTAMP', 'TNFRSF11A'))) + RotatedAxis()

## Identify DEGs in clusters ----
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_Clusters_BoneMarrow.xlsx')

## Subset a cluster to further differentiate ----
# Can't seem to get some of the myeloid, early progenitor, and B cells to split via clustering, so taking the problem cluster as a subset and re-clustering only those cells to divide properly
Idents(seu) <- seu$seurat_clusters
sub <- subset(seu, idents = c('16'))
DefaultAssay(sub) <- 'RNA'
sub <- DietSeurat(sub,
                  assays = 'RNA')
sub.list <- SplitObject(sub, split.by = "ID") # split by sample IDs
for (i in 1:length(sub.list)) { # normalize data using SCTransform method
  sub.list[[i]] <- SCTransform(sub.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

seu.features <- SelectIntegrationFeatures(sub.list, # select the genes to use for integration
                                          verbose = TRUE) 
sub.list <- PrepSCTIntegration(sub.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(sub.list, # identify anchors for integration from top 30 data dimensions
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features,
                                      k.score = min(table(sub$ID))-1, # fixes error: Error in idx[i, ] <- res[[i]][[1]] : number of items to replace is not a multiple of replacement length
                                      dims = 1:(min(table(sub$ID))-1)) #fixes error: Error in FindIntegrationAnchors(sub.list, normalization.method = "SCT",  : Max dimension too large: objects 4 contain fewer than 30 cells. Please specify a maximum dimensions that is less than the number of cells in any object (29).
sub <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT",
                     k.weight = min(table(sub$ID))) # fixes error: Error in idx[i, ] <- res[[i]][[1]] : number of items to replace is not a multiple of replacement length

rm(seu.features, seu.anchors, sub.list)

#### Calculate principle components
#Calculate 50 PCs:
sub <- RunPCA(sub, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

#Visualize PCs:
ElbowPlot(sub,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... we can use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering or take a numeric approach to identifying PCs as below

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- sub[["pca"]]@stdev / sum(sub[["pca"]]@stdev) * 100 # find standard deviation for each PC
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
sub <- RunUMAP(sub, 
               reduction = "pca",
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

sub <- RunTSNE(sub, 
               dims = PCdims,
               k.seed = 777,
               reduction = "pca") # create UMAP

DefaultAssay(sub) <- 'integrated'
sub <- FindNeighbors(sub, 
                     verbose = TRUE,
                     dims = PCdims,
                     assay = 'integrated') 
sub <- FindClusters(sub, 
                    resolution = 1)
DimPlot(sub, label = TRUE)

DefaultAssay(sub) <- 'RNA'
FeaturePlot(sub, features = c('AIF1', 'CSF1R', 'FLT3', 'CD79B', 'CD19', 'PAX5', 'KIT', 'LYZ', 'CST3', 'ENSSSCG00000036618'))
DotPlot(sub,
        features = c('CD19', 'CD79A', 'CD79B', 'PAX5', 'JCHAIN',
                     'CD3E', 'CD3G', 'CD247', 'CD2', 'TRDC', 'CD4', 'CD8B', 'CD8A', 'KLRB1',
                     'AIF1', 'CSF1R', 'CST3', 'CD14', 'ENSSSCG00000036618', 'LYZ', 'KIT', 'FLT3')) & RotatedAxis()

# After some gene query, we came up with annotations for our clusters:
Idents(sub) <- sub$seurat_clusters
sub$CellLineage <- sub$seurat_clusters
Idents(sub) <- sub$CellLineage
types <- c('B', 'Progenitor', 'Progenitor', 'Myeloid') 
names(types) <- levels(sub) # assign to cluster numbers
sub <- RenameIdents(sub, types) # change dataset identity to cell types in new subrat object
sub$CellLineage <- Idents(sub)
Idents(sub) <- sub$CellLineage
DimPlot(sub)

# Apply sub idents to seu data object:
## Merge annotations:
#Identify cell barcodes corresponding to different epithelial cell type predictions:
my <- rownames(sub@meta.data %>% filter(sub$CellLineage == 'Myeloid'))
b <- rownames(sub@meta.data %>% filter(sub$CellLineage == 'B'))
pre <- rownames(sub@meta.data %>% filter(sub$CellLineage == 'Progenitor'))

## Assign cell annotations ----
seu$seurat_clusters2 <- seu$seurat_clusters
Idents(seu) <- seu$seurat_clusters
seu$CellType <- seu$seurat_clusters
Idents(seu) <- seu$CellType
types <- c('cDC', 'gdT_CD2neg', 'TILC_CCL5pos', 'Mac/Mono_CD163pos', 'HBBhi', 'Mac/Mono_CD163neg_CD14neg',
           'Mac/Mono_CD163pos_Cycling', 'Mac/Mono_CD163neg_CD14pos', 'pDC', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos_Cycling', 
           'T/ILC_CCL5neg', 'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163pos_Cycling', 'Progenitor_Cycling', 'Mac/Mono_CD163neg_CD14pos',
           'Split', 'B', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163neg_CD14neg', 
           'Mac/Mono_CD163neg_CD14neg', 'cDC', 'Mac/Mono_CD163pos_Cycling', 'Progenitor', 'Mac/Mono_CD163pos_Cycling',
           'Mac/Mono_CD163pos', 'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163neg_CD14pos', 'Progenitor', 'Mac/Mono_CD163pos',
           'Mac/Mono_CD163pos', 'Early_myeloid_cycling', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos', 'pDC',
           'Early_myeloid_cycling', 'cDC_cycling', 'Mac/Mono_CD163neg_CD14neg', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos', 
           'cDC', 'ASC', 'pDC', 'B_DNTTpos', 'TILC_Cycling', 
           'ILC_KITpos', 'Osteoclast') 
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$CellType <- Idents(seu)
Idents(seu) <- seu$CellType
seu <- SetIdent(seu, cells = my, value = "cDC_cycling")
seu <- SetIdent(seu, cells = b, value = "B_DNTTpos")
seu <- SetIdent(seu, cells = pre, value = "Progenitor_Cycling")
seu$CellType <- Idents(seu)

## Visualize data ----
#DimPlot(seu, reduction = 'tsne', label = TRUE, cols = cols[1:(length(unique(seu$CellType)))])
DimPlot(seu, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(seu$CellType)))], group.by = 'CellType')
DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'seurat_clusters')

## Perform hierarchical clustering ----
#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our seusequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for seusequent steps
rm(pct, cumu, co1, co2, pcs)

Idents(seu) <- seu$CellType
seu <- BuildClusterTree(seu, 
                        dims = PCdims)
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('pDC', 'ASC', 'Progenitor_Cycling', 
                                            'Progenitor', 'B_DNTTpos', 'Early_myeloid_cycling',
                                            'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163neg_CD14neg',
                                            'Mac/Mono_CD163neg_CD14pos', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos_Cycling', 'HBBhi',
                                            'cDC_cycling', 'Osteoclast', 'cDC', 'B', 
                                            'gdT_CD2neg', 'TILC_CCL5pos', 'T/ILC_CCL5neg', 'ILC_KITpos',
                                            'TILC_Cycling'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)
Idents(seu) <- seu$CellType
levels(seu) <- rev(c('pDC', 'ASC', 'Progenitor_Cycling', 
                     'Progenitor', 'B_DNTTpos', 'Early_myeloid_cycling',
                     'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163neg_CD14neg',
                     'Mac/Mono_CD163neg_CD14pos', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos_Cycling', 'HBBhi',
                     'cDC_cycling', 'Osteoclast', 'cDC', 'B', 
                     'gdT_CD2neg', 'TILC_CCL5pos', 'T/ILC_CCL5neg', 'ILC_KITpos',
                     'TILC_Cycling')) # Reorder the clusters based on phylo order from hierarchical tree
seu$CellType <- Idents(seu)
DotPlot(seu,
        features = unique(c('CD3E', 'CD3G', 'CD2', 'TRDC', 'CD4', 'CD8A', 'CD8B', 
                            'CCL5', 'KLRK1', 'IL7R', 'KLF2', 'S1PR1', 'PCLAF', 'ENSSSCG00000026302', #MKI67
                            'CD79B', 'PAX5', 'CD19', 'JCHAIN', 'XBP1', 'PRDM1',
                            'AIF1', 'CST3', 'CSF1R', 'FLT3', 'CD14', 'ENSSSCG00000036618', 'CD163', 'NLRP3', 'CCR2',
                            'KIT', 'DNTT', 'HBB', 'CDH5', 'RPS6')),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()


SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_BoneMarrow_CellTypeAnnotated_HBBhiIncluded.h5Seurat')
SaveH5Seurat(sub, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_BoneMarrow_SubClustersAnnotated.h5Seurat')

## Remove HBBhi cluster ----
seu$CellType <- Idents(seu)
seu <- subset(seu, idents = 'HBBhi', invert = TRUE)
DimPlot(seu, reduction = 'tsne')

## Identify DEGs in cell types ----
DefaultAssay(seu) <- 'RNA'
de <- FindAllMarkers(seu)
write_xlsx(de, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/DEGLists/DEGs_CellTypes_BoneMarrow.xlsx')

## Redo processing after HBBhi removal ----
#Find variable features and integrate data:
DefaultAssay(seu) <- 'RNA'
Idents(seu) <- seu$CellType
seu.list <- SplitObject(seu, split.by = "ID") # split by sample IDs
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
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
               dims = PCdims,
               seed.use = 777,
               min.dist = 0.3,
               spread = 0.3) # create UMAP

seu <- RunTSNE(seu, 
               dims = PCdims,
               k.seed = 777,
               reduction = "pca") # create UMAP

Idents(seu) <- seu$CellType
seu <- BuildClusterTree(seu, 
                        dims = PCdims)
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('pDC', 'ASC', 'Progenitor_Cycling', 
                                            'Progenitor', 'B_DNTTpos', 'Early_myeloid_cycling',
                                            'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163neg_CD14neg',
                                            'Mac/Mono_CD163neg_CD14pos', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos_Cycling',
                                            'cDC_cycling', 'Osteoclast', 'cDC', 'B', 
                                            'gdT_CD2neg', 'TILC_CCL5pos', 'T/ILC_CCL5neg', 'ILC_KITpos',
                                            'TILC_Cycling'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)
Idents(seu) <- seu$CellType
levels(seu) <- rev(c('pDC', 'ASC', 'Progenitor_Cycling', 
                     'Progenitor', 'B_DNTTpos', 'Early_myeloid_cycling',
                     'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163neg_CD14neg',
                     'Mac/Mono_CD163neg_CD14pos', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos_Cycling', 
                     'cDC_cycling', 'Osteoclast', 'cDC', 'B', 
                     'gdT_CD2neg', 'TILC_CCL5pos', 'T/ILC_CCL5neg', 'ILC_KITpos',
                     'TILC_Cycling')) # Reorder the clusters based on phylo order from hierarchical tree
seu$CellType <- Idents(seu)
DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = unique(c('CD3E', 'CD3G', 'CD2', 'TRDC', 'CD4', 'CD8A', 'CD8B', 
                            'CCL5', 'KLRK1', 'IL7R', 'KLF2', 'S1PR1', 'PCLAF', 'ENSSSCG00000026302', #MKI67
                            'CD79B', 'PAX5', 'CD19', 'JCHAIN', 'XBP1', 'PRDM1',
                            'AIF1', 'CST3', 'CSF1R', 'FLT3', 'CD14', 'ENSSSCG00000036618', 'CD163', 'NLRP3', 'CCR2',
                            'KIT', 'DNTT', 'HBB', 'CDH5', 'RPS6')),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()

## Visualize data ----
#DimPlot(seu, reduction = 'tsne', label = TRUE, cols = cols[1:(length(unique(seu$CellType)))])
DimPlot(seu, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(seu$CellType)))], group.by = 'CellType')
DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'seurat_clusters')

## Save final dataset ----
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_BoneMarrow_CellTypeAnnotated_HBBhiRemoved.h5Seurat')

# Merge plots ----

## Load data ----
sp <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen_CellTypeAnnotated_HBBhiRemoved.h5Seurat')
ln <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_LymphNode_CellTypeAnnotated.h5Seurat')
th <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Thymus_CellTypeAnnotated.h5Seurat')
bm <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_BoneMarrow_CellTypeAnnotated_HBBhiRemoved.h5Seurat')

## Visualize tSNE plots ----
a <- DimPlot(bm, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(bm$CellType)))]) 
b <- DimPlot(th, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(th$CellType)))]) 
c <- DimPlot(ln, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(ln$CellType)))]) 
d <- DimPlot(sp, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(sp$CellType)))]) 
plot_grid(a, b, c, d,
          ncol = 2)

a <- DimPlot(bm, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(bm$CellType)))]) & NoLegend()
b <- DimPlot(th, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(th$CellType)))]) & NoLegend()
c <- DimPlot(ln, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(ln$CellType)))]) & NoLegend()
d <- DimPlot(sp, reduction = 'tsne', label = FALSE, cols = cols[1:(length(unique(sp$CellType)))]) & NoLegend()
plot_grid(a, b, c, d,
          ncol = 2)

## Visualize heirarchical clustering trees ----
seu <- bm
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('pDC', 'ASC', 'Progenitor_Cycling', 
                                            'Progenitor', 'B_DNTTpos', 'Early_myeloid_cycling',
                                            'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163neg_CD14neg',
                                            'Mac/Mono_CD163neg_CD14pos', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163pos_Cycling',
                                            'cDC_cycling', 'Osteoclast', 'cDC', 'B', 
                                            'gdT_CD2neg', 'TILC_CCL5pos', 'T/ILC_CCL5neg', 'ILC_KITpos',
                                            'TILC_Cycling'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)

seu <- th
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('pDC', 'Mac/Mono/cDC', 'ASC', 'Bcell', 
                                            'abT_Cycling', 'gdT_CD2neg_Cycling', 
                                            'gdT_CD2neg', 'gdT_CD2posCD8Aneg', 'CD4T_KLF2pos', 'CD8T_KLF2pos', 'CD4T_KLF2neg_CCR9neg',
                                            'CD4T_KLF2neg_CCR9pos', 'CD8T_KLF2neg', 'abT_CD1DposCD1Epos',
                                            'gdT_CD2pos_Cycling', 'ThymocyteDN_Cycling', 'ThymocyteDP_Cycling', 'ThymocyteDN', 'ThymocyteDP',
                                            'Progenitor', 'gdT_CD2posCD8Apos', 'TILC_CCL5pos'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)

seu <- ln
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                                            'Mac/Mono/cDC',
                                            'TILC_CCL5pos', 'ILC_KITpos',
                                            'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                                            'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)

seu <- sp
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('ASC', 'Mac/Mono', 'cDC', 'Bcell',
                                            'T/ILC_Cycling', 'ILC_NCR1posEOMESpos', 'ILC_Cytotoxic',
                                            'gdT_CD2pos', 'gdT_CD2neg', 'abT_CCL5neg', 'abT_Cytotoxic', 'abT_CCL5pos'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)

## Visualize feature plots w/ canonical genes ----
DefaultAssay(bm) <- 'RNA'
DefaultAssay(th) <- 'RNA'
DefaultAssay(ln) <- 'RNA'
DefaultAssay(sp) <- 'RNA'
FeaturePlot(bm, features = c('CD3G', 'CD79B', 'AIF1', 'KIT'), reduction = 'tsne', ncol = 2, cols = c('grey75', 'red3')) & NoAxes()
FeaturePlot(th, features = c('CD3G', 'CD79B', 'AIF1', 'KIT'), reduction = 'tsne', ncol = 2, cols = c('grey75', 'red3')) & NoAxes()
FeaturePlot(ln, features = c('CD3G', 'CD79B', 'AIF1', 'KIT'), reduction = 'tsne', ncol = 2, cols = c('grey75', 'red3')) & NoAxes()
FeaturePlot(sp, features = c('CD3G', 'CD79B', 'AIF1', 'KIT'), reduction = 'tsne', ncol = 2, cols = c('grey75', 'red3')) & NoAxes()

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] cowplot_1.1.3         scales_1.3.0          writexl_1.4.2         readxl_1.4.3          dplyr_1.1.4           SeuratDisk_0.0.0.9020 ggplot2_3.5.0        
#[8] SeuratObject_4.1.3    Seurat_4.3.0.1       

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.16             colorspace_2.1-0       rjson_0.2.21           deldir_1.0-9           ellipsis_0.3.2         ggridges_0.5.4         circlize_0.4.16       
#[8] GlobalOptions_0.1.2    clue_0.3-64            rstudioapi_0.15.0      spatstat.data_3.0-1    farver_2.1.1           leiden_0.4.3           listenv_0.9.1         
#[15] ggrepel_0.9.5          bit64_4.0.5            fansi_1.0.6            codetools_0.2-19       splines_4.2.2          doParallel_1.0.17      polyclip_1.10-4       
#[22] jsonlite_1.8.8         ica_1.0-3              cluster_2.1.4          png_0.1-8              uwot_0.1.16            shiny_1.8.0            sctransform_0.3.5     
#[29] spatstat.sparse_3.0-2  compiler_4.2.2         httr_1.4.7             Matrix_1.6-1           fastmap_1.1.1          lazyeval_0.2.2         cli_3.6.2             
#[36] later_1.3.2            htmltools_0.5.7        tools_4.2.2            igraph_1.5.1           gtable_0.3.4           glue_1.7.0             RANN_2.6.1            
#[43] reshape2_1.4.4         Rcpp_1.0.12            scattermore_1.2        cellranger_1.1.0       vctrs_0.6.5            ape_5.7-1              spatstat.explore_3.2-3
#[50] nlme_3.1-163           progressr_0.14.0       iterators_1.0.14       lmtest_0.9-40          spatstat.random_3.1-6  stringr_1.5.1          globals_0.16.2        
#[57] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3          future_1.33.1          MASS_7.3-60           
#[64] zoo_1.8-12             promises_1.2.1         spatstat.utils_3.0-3   parallel_4.2.2         RColorBrewer_1.1-3     ComplexHeatmap_2.12.1  reticulate_1.35.0     
#[71] pbapply_1.7-2          gridExtra_2.3          stringi_1.8.3          S4Vectors_0.34.0       foreach_1.5.2          BiocGenerics_0.42.0    shape_1.4.6.1         
#[78] rlang_1.1.3            pkgconfig_2.0.3        matrixStats_1.0.0      lattice_0.21-8         ROCR_1.0-11            purrr_1.0.2            tensor_1.5            
#[85] labeling_0.4.3         patchwork_1.2.0        htmlwidgets_1.6.4      bit_4.0.5              tidyselect_1.2.0       parallelly_1.37.1      RcppAnnoy_0.0.21      
#[92] plyr_1.8.9             magrittr_2.0.3         R6_2.5.1               IRanges_2.30.1         generics_0.1.3         pillar_1.9.0           withr_3.0.0           
#[99] fitdistrplus_1.1-11    survival_3.5-7         abind_1.4-5            sp_2.0-0               tibble_3.2.1           future.apply_1.11.1    crayon_1.5.2          
#[106] hdf5r_1.3.8            KernSmooth_2.23-22     utf8_1.2.4             spatstat.geom_3.2-5    plotly_4.10.4          GetoptLong_1.0.5       locfit_1.5-9.8        
#[113] grid_4.2.2             data.table_1.15.2      digest_0.6.34          xtable_1.8-4           tidyr_1.3.1            httpuv_1.6.14          stats4_4.2.2          
#[120] munsell_0.5.0          viridisLite_0.4.2  