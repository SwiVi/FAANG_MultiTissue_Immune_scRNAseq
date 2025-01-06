library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)
library(cowplot)
library(edgeR)

# Merge cell lineages ----

## Load data ----
sp <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Spleen_CellTypeAnnotated_HBBhiRemoved.h5Seurat')
ln <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_LymphNode_CellTypeAnnotated.h5Seurat')
th <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Thymus_CellTypeAnnotated.h5Seurat')
bm <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_BoneMarrow_CellTypeAnnotated_HBBhiRemoved.h5Seurat')

## T/ILC ----
unique(Idents(sp))
unique(Idents(ln))
unique(Idents(th))
unique(Idents(bm))

spT <- subset(sp, idents = c('abT_CCL5neg', 'gdT_CD2neg', 'ILC_Cytotoxic', 'abT_CCL5pos', 
                             'T/ILC_Cycling', 'ILC_NCR1posEOMESpos', 'gdT_CD2pos', 'abT_Cytotoxic'))
lnT <- subset(ln, idents = c('CD4T_KLF2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'CD4Tfh', 
                             'CD4T_KLF2neg', 'TILC_CCL5pos', 'gdT_CD2neg', 'ILC_KITpos', 'gdT_CD2pos'))
thT <- subset(th, idents = c('ThymocyteDP_Cycling', 'ThymocyteDP', 'gdT_CD2posCD8Aneg', 'gdT_CD2neg', 'CD8T_KLF2neg', 
                             'CD4T_KLF2pos', 'ThymocyteDN_Cycling', 'gdT_CD2posCD8Apos', 'gdT_CD2pos_Cycling',
                             'ThymocyteDN', 'CD4T_KLF2neg_CCR9neg', 'CD8T_KLF2pos', 'CD4T_KLF2neg_CCR9pos', 'TILC_CCL5pos',
                             'abT_CD1DposCD1Epos', 'gdT_CD2neg_Cycling', 'abT_Cycling'))
bmT <- subset(bm, idents = c('gdT_CD2neg', 'TILC_CCL5pos', 'ILC_KITpos', 'T/ILC_CCL5neg', 'TILC_Cycling'))

spT$CellType <- paste(spT$CellType, 'Spleen', sep = '_')
lnT$CellType <- paste(lnT$CellType, 'LymphNode', sep = '_')
thT$CellType <- paste(thT$CellType, 'Thymus', sep = '_')
bmT$CellType <- paste(bmT$CellType, 'BoneMarrow', sep = '_')

seu <- merge(x = spT, y = c(lnT, thT, bmT))

sub.list <- SplitObject(seu, split.by = "ID")
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
                                      anchor.features = seu.features)
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT")

seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

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

seu$CellType <- factor(seu$CellType)
Idents(seu) <- seu$CellType
seu$CellType2 <- seu$CellType
Idents(seu) <- seu$CellType2
levels(seu) <- rev(c('ThymocyteDN_Thymus', 'ThymocyteDP_Thymus', 'ThymocyteDN_Cycling_Thymus', 'ThymocyteDP_Cycling_Thymus',
                     'T/ILC_Cycling_Spleen', 'TILC_Cycling_BoneMarrow', 'abT_Cycling_Thymus', 
                     'gdT_CD2pos_Cycling_Thymus', 'gdT_CD2neg_Cycling_Thymus', 
                     'gdT_CD2neg_BoneMarrow', 'gdT_CD2neg_Spleen', 'gdT_CD2neg_LymphNode', 'gdT_CD2neg_Thymus',
                     'gdT_CD2posCD8Apos_Thymus', 'gdT_CD2posCD8Aneg_Thymus', 'gdT_CD2pos_Spleen',  'gdT_CD2pos_LymphNode', 
                     'abT_Cytotoxic_Spleen', 'abT_CCL5pos_Spleen', 
                     'TILC_CCL5pos_BoneMarrow', 'TILC_CCL5pos_LymphNode', 'TILC_CCL5pos_Thymus',
                     'abT_CCL5neg_Spleen', 'T/ILC_CCL5neg_BoneMarrow', 
                     'CD8T_KLF2pos_LymphNode', 'CD8T_KLF2pos_Thymus', 
                     'CD4T_KLF2pos_LymphNode', 'CD4T_KLF2pos_Thymus',
                     'CD8T_KLF2neg_LymphNode', 'CD8T_KLF2neg_Thymus',
                     'CD4T_KLF2neg_LymphNode', 'CD4Tfh_LymphNode', 
                     'CD4T_KLF2neg_CCR9pos_Thymus', 'CD4T_KLF2neg_CCR9neg_Thymus', 
                     'abT_CD1DposCD1Epos_Thymus', 
                     'ILC_Cytotoxic_Spleen', 'ILC_NCR1posEOMESpos_Spleen', 'ILC_KITpos_LymphNode', 'ILC_KITpos_BoneMarrow')) # Reorder the clusters 

seu$CellType2 <- Idents(seu)
Idents(seu) <- seu$CellType2

### Plot canonical genes ----
DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = c('CD3E', 'CD3G', 'CD3D', 'CD247', 'ZAP70',
                     'RAG1', 'RAG2', 'DNTT', 
                     'CD4', 'CD8B', 'CD8A', 'TRDC', 'CD2', 
                     'PCLAF', 'BIRC5', 'TOP2A', 'STMN1',
                     'ENSSSCG00000034914', 'RHEX', # 'ENSSSCG00000034914' = CD163L
                     'CCL5', 'KLRK1', 
                     'GZMB', 'GNLY',
                     'S1PR1', 'CCR7', 
                     'KLF2', 
                     'IL7R', 
                     'LEF1', 'SELL',
                     'CXCR5', 'CD40LG', 
                     'ICOS', 'CTLA4', 'PDCD1', 'PRDM1', 'IL2RA', 
                     'CCR9',
                     'CD1D', 'CD1E', 'GATA3', 'CCR4',
                     
                     'KIT', 'CD69', 'EOMES', 'NCR1', 'KLRB1', 'RORC'),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()

### Pseudobulk visualization ----

av.exp <- AverageExpression(seu, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample
counts <- as.matrix(av.exp@assays$RNA@data)
cor.exp <- as.data.frame(cor(counts, method = 'spearman'))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, levels(Idents(seu)))
cor.df$x <- factor(cor.df$x,levels = levels(Idents(seu)))
cor.df$y <- factor(cor.df$y,levels = levels(Idents(seu)))
cor.df$correlation <- round(cor.df$correlation, digits = 2)
cor.df$x <- factor(cor.df$x, levels=rev(c('ThymocyteDN_Thymus', 'ThymocyteDN_Cycling_Thymus', 'ThymocyteDP_Thymus', 'ThymocyteDP_Cycling_Thymus',
                                      'T/ILC_Cycling_Spleen', 'TILC_Cycling_BoneMarrow', 'abT_Cycling_Thymus', 
                                      'gdT_CD2neg_Cycling_Thymus', 'gdT_CD2pos_Cycling_Thymus', 
                                      'gdT_CD2neg_BoneMarrow', 'gdT_CD2neg_Spleen', 'gdT_CD2neg_LymphNode', 'gdT_CD2neg_Thymus',
                                      'gdT_CD2pos_Spleen', 'gdT_CD2posCD8Apos_Thymus', 'gdT_CD2posCD8Aneg_Thymus', 'gdT_CD2pos_LymphNode',
                                      'abT_Cytotoxic_Spleen', 'abT_CCL5pos_Spleen', 'abT_CCL5neg_Spleen',
                                      'TILC_CCL5pos_BoneMarrow', 'TILC_CCL5pos_LymphNode', 'TILC_CCL5pos_Thymus',
                                      'CD8T_KLF2pos_Thymus', 'CD8T_KLF2pos_LymphNode', 'CD4T_KLF2pos_Thymus', 'CD4T_KLF2pos_LymphNode',
                                      'T/ILC_CCL5neg_BoneMarrow', 
                                      'CD8T_KLF2neg_Thymus', 'CD8T_KLF2neg_LymphNode', 'CD4T_KLF2neg_CCR9pos_Thymus', 'CD4T_KLF2neg_LymphNode',
                                      'CD4Tfh_LymphNode', 'CD4T_KLF2neg_CCR9neg_Thymus', 
                                      'abT_CD1DposCD1Epos_Thymus', 
                                      'ILC_Cytotoxic_Spleen', 'ILC_NCR1posEOMESpos_Spleen', 'ILC_KITpos_LymphNode', 'ILC_KITpos_BoneMarrow')))
cor.df$y <- factor(cor.df$y, levels=rev(c('ThymocyteDN_Thymus', 'ThymocyteDN_Cycling_Thymus', 'ThymocyteDP_Thymus', 'ThymocyteDP_Cycling_Thymus',
                                      'T/ILC_Cycling_Spleen', 'TILC_Cycling_BoneMarrow', 'abT_Cycling_Thymus', 
                                      'gdT_CD2neg_Cycling_Thymus', 'gdT_CD2pos_Cycling_Thymus', 
                                      'gdT_CD2neg_BoneMarrow', 'gdT_CD2neg_Spleen', 'gdT_CD2neg_LymphNode', 'gdT_CD2neg_Thymus',
                                      'gdT_CD2pos_Spleen', 'gdT_CD2posCD8Apos_Thymus', 'gdT_CD2posCD8Aneg_Thymus', 'gdT_CD2pos_LymphNode',
                                      'abT_Cytotoxic_Spleen', 'abT_CCL5pos_Spleen', 'abT_CCL5neg_Spleen',
                                      'TILC_CCL5pos_BoneMarrow', 'TILC_CCL5pos_LymphNode', 'TILC_CCL5pos_Thymus',
                                      'CD8T_KLF2pos_Thymus', 'CD8T_KLF2pos_LymphNode', 'CD4T_KLF2pos_Thymus', 'CD4T_KLF2pos_LymphNode',
                                      'T/ILC_CCL5neg_BoneMarrow', 
                                      'CD8T_KLF2neg_Thymus', 'CD8T_KLF2neg_LymphNode', 'CD4T_KLF2neg_CCR9pos_Thymus', 'CD4T_KLF2neg_LymphNode',
                                      'CD4Tfh_LymphNode', 'CD4T_KLF2neg_CCR9neg_Thymus', 
                                      'abT_CD1DposCD1Epos_Thymus', 
                                      'ILC_Cytotoxic_Spleen', 'ILC_NCR1posEOMESpos_Spleen', 'ILC_KITpos_LymphNode', 'ILC_KITpos_BoneMarrow')))

ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()+
  scale_fill_gradientn(colours = c('beige', 'yellow','orange','red', 'darkred'), oob = squish, limits = c(min(cor.df$correlation), 1))

DGE <- DGEList(counts = counts, genes = rownames(counts), group = colnames(counts)) # make into edgeR DGE object
plotMDS(DGE, 
        top = length(rownames(av.exp[["RNA"]])), #consider all genes 
        labels=levels(Idents(seu)))
plotMDS(DGE, 
        top = 500, 
        labels=levels(Idents(seu)))


### Save final dataset ----
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_AllTissues_CellTypeMerge_TILC.h5Seurat')

## B ----
unique(Idents(sp))
unique(Idents(ln))
unique(Idents(th))
unique(Idents(bm))

spB <- subset(sp, idents = c('Bcell', 'ASC'))
lnB <- subset(ln, idents = c('B_AICDAneg', 'B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'ASC'))
thB <- subset(th, idents = c('Bcell', 'ASC'))
bmB <- subset(bm, idents = c('B', 'ASC', 'B_DNTTpos'))

spB$CellType <- paste(spB$CellType, 'Spleen', sep = '_')
lnB$CellType <- paste(lnB$CellType, 'LymphNode', sep = '_')
thB$CellType <- paste(thB$CellType, 'Thymus', sep = '_')
bmB$CellType <- paste(bmB$CellType, 'BoneMarrow', sep = '_')

seu <- merge(x = spB, y = c(lnB, thB, bmB))

sub.list <- SplitObject(seu, split.by = "ID")
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
                                      k.score = min(table(seu$ID))-1, # fixes error: Error in idx[i, ] <- res[[i]][[1]] : number of items to replace is not a multiple of replacement length
                                      dims = 1:(min(table(seu$ID))-1)) #fixes error: Error in FindIntegrationAnchors(sub.list, normalization.method = "SCT",  : Max dimension too large: objects 4 contain fewer than 30 cells. Please specify a maximum dimensions that is less than the number of cells in any object (25).
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT",
                     k.weight = min(table(seu$ID)))

seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

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

seu$CellType <- factor(seu$CellType)
Idents(seu) <- seu$CellType
seu$CellType2 <- seu$CellType
Idents(seu) <- seu$CellType2
levels(seu) <- rev(c("ASC_BoneMarrow", "ASC_Thymus", "ASC_LymphNode", "ASC_Spleen", 
                     "B_BoneMarrow", "Bcell_Thymus", "Bcell_Spleen", 
                     "B_AICDAneg_LymphNode", "B_AICDApos_LymphNode", 
                     "B_AICDApos_Cycling_LymphNode", "B_AICDAneg_Cycling_LymphNode",
                     "B_DNTTpos_BoneMarrow")) # Reorder the clusters based on phylo order from hierarchical tree
seu$CellType2 <- Idents(seu)
Idents(seu) <- seu$CellType2

### Plot canonical genes ----
DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = c('JCHAIN', 'XBP1', 'PRDM1', 'IRF4',
                     'PAX5', 'CD79A', 'CD79B', 'CD19', 'MS4A1',
                     'GPR183', 'FCER2', 'KLF2', 
                     'CCR7', 'SELL', 'S1PR1', 
                     'AICDA', 'BCL7A', 'HMCES', 'NEIL1', 'RGS13', 'EAF2', 'DAPP1', 'S1PR2', 'CD86', # S1PR2 = GC confinement; BCL7A = chromatin remodeling
                     'BCL6', 'GCSAM', 
                     'PCLAF', 'BIRC5', 'TOP2A', 'STMN1',
                     'DNTT', 'LXN', 'ERG', 'EBF1', 'VPREB1', 'H2AFY', 'CD34', 'FLT3'),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) & RotatedAxis()

### Pseudobulk visualization ----

av.exp <- AverageExpression(seu, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample
counts <- as.matrix(av.exp@assays$RNA@data)
cor.exp <- as.data.frame(cor(counts, method = 'spearman'))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, levels(Idents(seu)))
cor.df$x <- factor(cor.df$x,levels = levels(Idents(seu)))
cor.df$y <- factor(cor.df$y,levels = levels(Idents(seu)))
cor.df$correlation <- round(cor.df$correlation, digits = 2)
cor.df$x <- factor(cor.df$x, levels=rev(c("ASC_BoneMarrow", "ASC_Thymus", "ASC_LymphNode", "ASC_Spleen", 
                                          "B_BoneMarrow", "Bcell_Thymus", "Bcell_Spleen", 
                                          "B_AICDAneg_LymphNode", "B_AICDApos_LymphNode", 
                                          "B_AICDApos_Cycling_LymphNode", "B_AICDAneg_Cycling_LymphNode",
                                          "B_DNTTpos_BoneMarrow")))
cor.df$y <- factor(cor.df$y, levels=rev(c("ASC_BoneMarrow", "ASC_Thymus", "ASC_LymphNode", "ASC_Spleen", 
                                          "B_BoneMarrow", "Bcell_Thymus", "Bcell_Spleen", 
                                          "B_AICDAneg_LymphNode", "B_AICDApos_LymphNode", 
                                          "B_AICDApos_Cycling_LymphNode", "B_AICDAneg_Cycling_LymphNode",
                                          "B_DNTTpos_BoneMarrow")))

ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()+
  scale_fill_gradientn(colours = c('beige', 'yellow','orange','red', 'darkred'), oob = squish, limits = c(min(cor.df$correlation), 1))

DGE <- DGEList(counts = counts, genes = rownames(counts), group = colnames(counts)) # make into edgeR DGE object
plotMDS(DGE, 
        top = length(rownames(av.exp[["RNA"]])), #consider all genes 
        labels=levels(Idents(seu)))
plotMDS(DGE, 
        top = 500, 
        labels=levels(Idents(seu)))

### Save final dataset ----
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_AllTissues_CellTypeMerge_B.h5Seurat')

## Myeloid ----
unique(Idents(sp))
unique(Idents(ln))
unique(Idents(th))
unique(Idents(bm))

spM <- subset(sp, idents = c('Mac/Mono', 'cDC'))
lnM <- subset(ln, idents = c('Mac/Mono/cDC'))
thM <- subset(th, idents = c('Mac/Mono/cDC', 'pDC'))
bmM <- subset(bm, idents = c('cDC', 'Mac/Mono_CD163pos_Cycling', 'Mac/Mono_CD163neg_CD14pos', 
                             'Mac/Mono_CD163neg_Cycling', 'Mac/Mono_CD163pos', 'Mac/Mono_CD163neg_CD14neg', 
                             'Early_myeloid_cycling', 'cDC_cycling', 'pDC', 'Osteoclast'))

spM$CellType <- paste(spM$CellType, 'Spleen', sep = '_')
lnM$CellType <- paste(lnM$CellType, 'LymphNode', sep = '_')
thM$CellType <- paste(thM$CellType, 'Thymus', sep = '_')
bmM$CellType <- paste(bmM$CellType, 'BoneMarrow', sep = '_')

seu <- merge(x = spM, y = c(lnM, thM, bmM))

sub.list <- SplitObject(seu, split.by = "ID")
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
                                      k.score = min(table(seu$ID))-1, # fixes error: Error in idx[i, ] <- res[[i]][[1]] : number of items to replace is not a multiple of replacement length
                                      dims = 1:(min(table(seu$ID))-1)) #fixes error: Error in FindIntegrationAnchors(sub.list, normalization.method = "SCT",  : Max dimension too large: objects 4 contain fewer than 30 cells. Please specify a maximum dimensions that is less than the number of cells in any object (25).
seu <- IntegrateData(seu.anchors, # integrate data
                     normalization.method = "SCT",
                     k.weight = min(table(seu$ID)))

seu <- RunPCA(seu, # run PCA analysis for 50 dimensions of the data
              npcs = 50, 
              verbose = TRUE) 

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

seu$CellType <- factor(seu$CellType)
Idents(seu) <- seu$CellType
seu$CellType2 <- seu$CellType
Idents(seu) <- seu$CellType2
levels(seu) <- rev(c('Osteoclast_BoneMarrow', 'Early_myeloid_cycling_BoneMarrow', 
                     'cDC_Spleen', 'cDC_BoneMarrow', 'cDC_cycling_BoneMarrow',
                     'Mac/Mono/cDC_Thymus', 'Mac/Mono/cDC_LymphNode', 'Mac/Mono_Spleen',
                     'Mac/Mono_CD163pos_BoneMarrow', 'Mac/Mono_CD163pos_Cycling_BoneMarrow',
                     'Mac/Mono_CD163neg_Cycling_BoneMarrow', 'Mac/Mono_CD163neg_CD14neg_BoneMarrow', 'Mac/Mono_CD163neg_CD14pos_BoneMarrow', 
                    "pDC_BoneMarrow", 'pDC_Thymus')) # Reorder the clusters based on phylo order from hierarchical tree
seu$CellType2 <- Idents(seu)
Idents(seu) <- seu$CellType2

### Plot canonical genes ----
DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = c('AIF1', 'CST3', 
                     'APOE', 'ATP6V0D2', 'NFATC1', 'MMP9', 'CALCR', 'SIGLEC15', 'ACP5', 'DCSTAMP', 'OCSTAMP', 'TNFRSF11A',
                     'KIT', 
                     'PCLAF', 'BIRC5', 'TOP2A', 'STMN1',
                     'SLA-DQB1', 'HLA-DRA',
                     'ENSSSCG00000036618', 'CD163', 'CD14', 'CD68', 'NLRP3', 'TLR4', 'CSF1R', 'CCR2',
                     'TCF4', 'IRF8', 'CD93', 'CLEC12A', 'XBP1', 'CD4', 'CD8B'),
        cols = c('gold', 'darkmagenta'), col.min = -1.5, col.max = 2) + RotatedAxis()

### Pseudobulk visualization ----

av.exp <- AverageExpression(seu, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample
counts <- as.matrix(av.exp@assays$RNA@data)
cor.exp <- as.data.frame(cor(counts, method = 'spearman'))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, levels(Idents(seu)))
cor.df$x <- factor(cor.df$x,levels = levels(Idents(seu)))
cor.df$y <- factor(cor.df$y,levels = levels(Idents(seu)))
cor.df$correlation <- round(cor.df$correlation, digits = 2)
cor.df$x <- factor(cor.df$x, levels=rev(c("pDC_BoneMarrow", 'pDC_Thymus', 'cDC_Spleen', 'cDC_BoneMarrow', 'cDC_cycling_BoneMarrow',
                                          'Mac/Mono/cDC_Thymus', 'Mac/Mono/cDC_LymphNode', 'Mac/Mono_Spleen',
                                          'Mac/Mono_CD163pos_BoneMarrow', 'Mac/Mono_CD163pos_Cycling_BoneMarrow',
                                          'Mac/Mono_CD163neg_Cycling_BoneMarrow', 'Mac/Mono_CD163neg_CD14neg_BoneMarrow', 'Mac/Mono_CD163neg_CD14pos_BoneMarrow', 
                                          'Early_myeloid_cycling_BoneMarrow', 'Osteoclast_BoneMarrow')))
cor.df$y <- factor(cor.df$y, levels=rev(c("pDC_BoneMarrow", 'pDC_Thymus', 'cDC_Spleen', 'cDC_BoneMarrow', 'cDC_cycling_BoneMarrow',
                                          'Mac/Mono/cDC_Thymus', 'Mac/Mono/cDC_LymphNode', 'Mac/Mono_Spleen',
                                          'Mac/Mono_CD163pos_BoneMarrow', 'Mac/Mono_CD163pos_Cycling_BoneMarrow',
                                          'Mac/Mono_CD163neg_Cycling_BoneMarrow', 'Mac/Mono_CD163neg_CD14neg_BoneMarrow', 'Mac/Mono_CD163neg_CD14pos_BoneMarrow', 
                                          'Early_myeloid_cycling_BoneMarrow', 'Osteoclast_BoneMarrow')))

ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()+
  scale_fill_gradientn(colours = c('beige', 'yellow','orange','red', 'darkred'), oob = squish, limits = c(min(cor.df$correlation), 1))

DGE <- DGEList(counts = counts, genes = rownames(counts), group = colnames(counts)) # make into edgeR DGE object
plotMDS(DGE, 
        top = length(rownames(av.exp[["RNA"]])), #consider all genes 
        labels=levels(Idents(seu)))
plotMDS(DGE, 
        top = 500, 
        labels=levels(Idents(seu)))

### Save final dataset ----
SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_AllTissues_CellTypeMerge_Myeloid.h5Seurat')

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
#[4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] edgeR_3.38.4          limma_3.52.4          cowplot_1.1.3         scales_1.3.0         
#[5] writexl_1.4.2         readxl_1.4.3          dplyr_1.1.4           SeuratDisk_0.0.0.9020
#[9] ggplot2_3.5.0         SeuratObject_4.1.3    Seurat_4.3.0.1       

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.16             colorspace_2.1-0       rjson_0.2.21           deldir_1.0-9          
#[5] ellipsis_0.3.2         ggridges_0.5.4         circlize_0.4.16        GlobalOptions_0.1.2   
#[9] clue_0.3-64            rstudioapi_0.15.0      spatstat.data_3.0-1    leiden_0.4.3          
#[13] listenv_0.9.1          ggrepel_0.9.5          bit64_4.0.5            fansi_1.0.6           
#[17] codetools_0.2-19       splines_4.2.2          doParallel_1.0.17      polyclip_1.10-4       
#[21] jsonlite_1.8.8         ica_1.0-3              cluster_2.1.4          png_0.1-8             
#[25] uwot_0.1.16            shiny_1.8.0            sctransform_0.3.5      spatstat.sparse_3.0-2 
#[29] compiler_4.2.2         httr_1.4.7             Matrix_1.6-1           fastmap_1.1.1         
#[33] lazyeval_0.2.2         cli_3.6.2              later_1.3.2            htmltools_0.5.7       
#[37] tools_4.2.2            igraph_1.5.1           gtable_0.3.4           glue_1.7.0            
#[41] RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.12            scattermore_1.2       
#[45] cellranger_1.1.0       vctrs_0.6.5            ape_5.7-1              spatstat.explore_3.2-3
#[49] nlme_3.1-163           progressr_0.14.0       iterators_1.0.14       lmtest_0.9-40         
#[53] spatstat.random_3.1-6  stringr_1.5.1          globals_0.16.2         mime_0.12             
#[57] miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3         
#[61] future_1.33.1          MASS_7.3-60            zoo_1.8-12             promises_1.2.1        
#[65] spatstat.utils_3.0-3   parallel_4.2.2         RColorBrewer_1.1-3     ComplexHeatmap_2.12.1 
#[69] reticulate_1.35.0      pbapply_1.7-2          gridExtra_2.3          stringi_1.8.3         
#[73] S4Vectors_0.34.0       foreach_1.5.2          BiocGenerics_0.42.0    shape_1.4.6.1         
#[77] rlang_1.1.3            pkgconfig_2.0.3        matrixStats_1.0.0      lattice_0.21-8        
#[81] ROCR_1.0-11            purrr_1.0.2            tensor_1.5             patchwork_1.2.0       
#[85] htmlwidgets_1.6.4      bit_4.0.5              tidyselect_1.2.0       parallelly_1.37.1     
#[89] RcppAnnoy_0.0.21       plyr_1.8.9             magrittr_2.0.3         R6_2.5.1              
#[93] IRanges_2.30.1         generics_0.1.3         pillar_1.9.0           withr_3.0.0           
#[97] fitdistrplus_1.1-11    survival_3.5-7         abind_1.4-5            sp_2.0-0              
#[101] tibble_3.2.1           future.apply_1.11.1    crayon_1.5.2           hdf5r_1.3.8           
#[105] KernSmooth_2.23-22     utf8_1.2.4             spatstat.geom_3.2-5    plotly_4.10.4         
#[109] GetoptLong_1.0.5       locfit_1.5-9.8         grid_4.2.2             data.table_1.15.2     
#[113] digest_0.6.34          xtable_1.8-4           tidyr_1.3.1            httpuv_1.6.14         
#[117] stats4_4.2.2           munsell_0.5.0          viridisLite_0.4.2
