library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(readxl)
library(future)
library(NMF) 
library(uwot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ComplexHeatmap)
library(viridis)
library(writexl)
library(CSOmapR)
library(Matrix)
library(plot3D)

#Color scheme to use in future plots:
cols = c('cyan4', 'gold3', 'chartreuse4', 'deeppink4', 'sandybrown',
         'cornflowerblue', 'navy', 'lightpink', 'salmon', 
         'deepskyblue2', 'tan4', 'mediumpurple1', 
         'darkgreen', 'gray50', 'darkmagenta', 'red', 'hotpink', 'khaki', 
         'orange3', 'limegreen', 'cadetblue3', 'firebrick', 'deepskyblue4', 'darkseagreen', 'burlywood3', 'black',
         'goldenrod3', 'blue', 'deeppink', 'lightskyblue3', 'mistyrose3', 'purple')

# Load Seurat object:
seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_LymphNode_CellTypeAnnotated.h5Seurat')
DefaultAssay(seu) <- 'RNA'

# Convert to human ortholog genes:
# Load in ortho genes
orthoGenes <- read.delim("/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/PigToHuman_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Human.homology.type == 'ortholog_one2one') # subset to only one to one orthologs

# Filter  data to include only one-to-one gene orthologs & convert to human gene symbols:
genes <- as.data.frame(rownames(seu[['RNA']]@data)) # extract pig gene names from dataset
colnames(genes) <- 'gene'
pigGenes <- read_excel('/home/Jayne.Wiarda/SI_PP_SC_ST/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalList <-gsub("_", "-", pigGenes$FinalList) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[pigGenes$FinalList %in% genes$gene, ] # slim down to only genes in our dataset
orthos <- intersect(pigGenes$ENSID, orthoGenes$Gene.stable.ID) # find which genes are one-to-one orthologs
#length(orthos) # how many genes are orthologs?
pigGenes <- pigGenes[pigGenes$ENSID %in% orthos, ]
#dim(pigGenes)
orthoGenes <- orthoGenes[orthoGenes$Gene.stable.ID %in% pigGenes$ENSID, ] # slim down to only ortho genes in our dataset
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Human.gene.stable.ID, .keep_all = TRUE)  # retain only unique combinations of pig & human Ensembl IDs, ignoring transcript IDs
#dim(orthoGenes) # should not have same number of rows as in pigGenes
counts <- seu[['RNA']]@data[rownames(seu[['RNA']]@data) %in% pigGenes$FinalList,]
pigGenes <- pigGenes %>% arrange(factor(FinalList, levels = rownames(counts))) # arrange pigGenes in same order as counts
orthoGenes <- orthoGenes %>% arrange(factor(Gene.stable.ID, levels = pigGenes$ENSID)) # arrange orthoGenes in same order as pigGenes (and consequently counts)
rownames(counts) <- orthoGenes$Human.gene.name # change pig genes to human gene names

# Do a custom bit here to incorporate more MHC-II SLA genes, since some are one-to-many orthologs but we are highly interested in MHC-II related pathways for this work...
add <- setdiff(rownames(seu[['RNA']]@data)[which(grepl("SLA-D", rownames(seu[['RNA']]@data)))],
               pigGenes$FinalList[which(grepl("SLA-D", pigGenes$FinalList))]) # ID MHC-II genes from Seurat object that weren't recognized as one-to-one orthologs
counts2 <- seu[['RNA']]@data[rownames(seu[['RNA']]@data) %in% add,]
rownames(counts2) <- gsub('SLA-','HLA-', rownames(counts2)) # change to human gene symbol IDs
counts <- rbind(counts, counts2)

meta <- seu@meta.data
meta$samples <- meta$SampleID

# Cell-cell signaling ----

# Create CellChat object:
cellchat <- createCellChat(object = counts, # use new Seurat object with human gene names
                           meta = meta,
                           group.by = "CellType") # set cell identities to cell type annotations

# Set cell interaction database:
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # look at only cell-cell contact
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocess expression data
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probabilities and network inferences:
#cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProb(cellchat, type =  "truncatedMean", trim = 0.1) # count as zero expression if in <10% of annotated cell type (default = 25%)

# Further process & visualize individual CellChat objects ----


# Infer cell signaling pathway communication & calculate aggregated communication network:
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow=c(1,1))
g1 <- netVisual_heatmap(cellchat, measure = 'count', color.heatmap = "YlOrRd", color.use = cols[1:length(unique(seu$CellType))])
g2 <- netVisual_heatmap(cellchat, measure = 'weight', color.heatmap = "YlOrRd", color.use = cols[1:length(unique(seu$CellType))])
g1 + g2
#levels(cellchat@idents)
#netVisual_bubble(cellchat, targets.use = c(1:5), remove.isolate = FALSE)
#netVisual_bubble(cellchat, sources.use = c(1:5), remove.isolate = FALSE)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
## B cells receiving only

ht1 <- netAnalysis_signalingRole_heatmap(cellchat,
                                         pattern = "incoming", color.heatmap = "YlOrRd", color.use = cols[1:length(unique(seu$CellType))])
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,
                                         pattern = "outgoing", color.heatmap = "YlOrRd", color.use = cols[1:length(unique(seu$CellType))])
ht1 + ht2
netAnalysis_signalingRole_heatmap(cellchat, pattern = 'all', color.heatmap = "YlOrRd", color.use = cols[1:length(unique(seu$CellType))])

# signaling contributions for specific pathways:
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10, color.heatmap = "YlOrRd", signaling = 'CD45', color.use = cols[1:length(unique(seu$CellType))])
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10, color.heatmap = "YlOrRd", signaling = 'MHC-II', color.use = cols[1:length(unique(seu$CellType))])
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10, color.heatmap = "YlOrRd", signaling = 'CD40', color.use = cols[1:length(unique(seu$CellType))])
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10, color.heatmap = "YlOrRd", signaling = 'CD86', color.use = cols[1:length(unique(seu$CellType))])
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10, color.heatmap = "YlOrRd", signaling = 'CD80', color.use = cols[1:length(unique(seu$CellType))])
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10, color.heatmap = "YlOrRd", signaling = 'PDL2', color.use = cols[1:length(unique(seu$CellType))])
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10, color.heatmap = "YlOrRd", signaling = 'PD-L1', color.use = cols[1:length(unique(seu$CellType))])

# L-R interactions for particular signaling pathways:
netVisual_chord_gene(cellchat, lab.cex = 0.5, slot.name = 'net', color.use = cols[1:length(unique(seu$CellType))], signaling = 'CD40')
netVisual_chord_gene(cellchat, lab.cex = 0.5, slot.name = 'net', color.use = cols[1:length(unique(seu$CellType))], signaling = 'MHC-II')
netVisual_chord_gene(cellchat, lab.cex = 0.5, slot.name = 'net', color.use = cols[1:length(unique(seu$CellType))], signaling = 'CD45')

netVisual_aggregate(cellchat, signaling = 'MHC-II', layout = "chord", color.use = cols[1:length(unique(seu$CellType))])
netVisual_aggregate(cellchat, signaling = 'MHC-II', layout = "circle", color.use = cols[1:length(unique(seu$CellType))])
netVisual_aggregate(cellchat, signaling = 'MHC-II', layout = "hierarchy", color.use = cols[1:length(unique(seu$CellType))], vertex.receiver = c(5:7))

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Save CellChat object:
saveRDS(cellchat, file = '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/CellChatLN.rds')

# 3D plot interaction space ----

# start by making a data frame of ligand/receptor pairs from CellChat database. Only use cell-cell contact interactions here.
#LR <- data.frame(CellChatDB.use$interaction$ligand, CellChatDB.use$interaction$receptor)
#colnames(LR) <- c('ligand', 'receptor')
#write_xlsx(LR, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/LigandReceptorPairs_CellChatCellCellContact.xlsx')
LR <- read_xlsx('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/41422_2020_353_MOESM14_ESM.xlsx') # read in CellChat human DB LR pairs from CSOmap manuscript
LR <- data.frame(LR[,1:2])
colnames(LR) <- c('ligand', 'receptor')
TPM <- counts

# Make affinity matrix of interactions
affinityMat = getAffinityMat(TPM, LR, verbose = T)
coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 500,
  verbose = T,
  df = 0.5)

coords = data.frame(coords_res$Y)
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
write_xlsx(coords, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/SpatialReconstructionCoordinates_LymphNode_df0.5_ManuscriptCellChatLRs.xlsx')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))
labelData <- data.frame(colnames(seu), seu$CellType)
colnames(labelData) <- c('Barcode', 'CellType')
identical(labelData$Barcode, rownames(coords)) # make sure cell barcode order matches!
join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)
write_xlsx(cellinfo_tbl, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/SpatialReconstructionCoordinatesAndAnnotation_LymphNode_df0.5_ManuscriptCellChatLRs.xlsx')
density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)
p_3Ddensity = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
p_3Ddensity
p_3Ddensity = plot3D(cellinfo_tbl, color_by = "CellType", title = "3D density",
                     colors = c(rep('red', 15)))
p_3Ddensity
density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

# replot in Seurat object
seu$CellType = factor(seu$CellType, 
                      levels=rev(c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                                   'Mac/Mono/cDC',
                                   'TILC_CCL5pos', 'ILC_KITpos',
                                   'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                                   'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg')))
rownames(cellinfo_tbl) <- cellinfo_tbl$cellName
identical(rownames(cellinfo_tbl), colnames(seu)) # TRUE
seu <- AddMetaData(seu, cellinfo_tbl) # add x,y,z coordinates as metadata
FeatureScatter(seu, 'x', 'y', group.by = 'CellType', cols = cols[1:length(unique(seu$CellType))], shuffle = TRUE, plot.cor = FALSE) + xlim(-0.02, 0.02) + ylim(-0.02, 0.02)
FeatureScatter(seu, 'x', 'z', group.by = 'CellType', cols = cols[1:length(unique(seu$CellType))], shuffle = TRUE, plot.cor = FALSE) + xlim(-0.02, 0.02) + ylim(-0.02, 0.02)
FeatureScatter(seu, 'z', 'y', group.by = 'CellType', cols = cols[1:length(unique(seu$CellType))], shuffle = TRUE, plot.cor = FALSE) + xlim(-0.02, 0.02) + ylim(-0.02, 0.02)

ggplot(cellinfo_tbl, aes(x=x, y=y)) + 
  geom_density_2d_filled(alpha = 0.5) + 
  xlim(-0.02, 0.02) +
  ylim(-0.02, 0.02) + 
  theme_bw() +
  geom_density_2d(linewidth = 0.25, colour = "black")
ggplot(cellinfo_tbl, aes(x=x, y=z)) + 
  geom_density_2d_filled(alpha = 0.5, contour_var	=) + 
  xlim(-0.02, 0.02) +
  ylim(-0.02, 0.02) + 
  theme_bw() +
  geom_density_2d(linewidth = 0.25, colour = "black")
ggplot(cellinfo_tbl, aes(x=z, y=y)) + 
  geom_density_2d_filled(alpha = 0.5) + 
  xlim(-0.02, 0.02) +
  ylim(-0.02, 0.02) + 
  theme_bw() +
  geom_density_2d(linewidth = 0.25, colour = "black")

# Get and plot colocalization significance ----
coords <- read_xlsx('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/SpatialReconstructionCoordinates_LymphNode_df0.5_ManuscriptCellChatLRs.xlsx')
cellinfo_tbl <- read_xlsx('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/SpatialReconstructionCoordinatesAndAnnotation_LymphNode_df0.5_ManuscriptCellChatLRs.xlsx')
coords <- as.matrix(coords)
rownames(coords) <- colnames(TPM)
identical(rownames(coords), cellinfo_tbl$cellName) # make sure this is TRUE!!
signif_results = getSignificance(coords, labels = cellinfo_tbl$CellType, verbose = T)
coloc <- signif_results$pvalue_tbl
coloc$CellType1 <- sub('---.*', '', coloc$cluster_pair)
coloc$CellType2 <- sub('.*---', '', coloc$cluster_pair)
write_xlsx(coloc, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/CellInteractions/SpatialReconstructionCoordinates_LymphNode_df0.5_ManuscriptCellChatColocalization.xlsx')

coloc$CellType1 = factor(coloc$CellType1, 
                         levels=rev(c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                                      'Mac/Mono/cDC',
                                      'TILC_CCL5pos', 'ILC_KITpos',
                                      'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                                      'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg')))
coloc$CellType2 = factor(coloc$CellType2, 
                         levels=c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                                  'Mac/Mono/cDC',
                                  'TILC_CCL5pos', 'ILC_KITpos',
                                  'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                                  'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg'))
coloc2 <- coloc[,c(1,2, 3, 5, 4)] # switch cell type columns
colnames(coloc2) <- colnames(coloc)
coloc <- rbind(coloc, coloc2)

ggplot(coloc, aes(CellType1, CellType2, fill= p.value)) + 
  geom_tile() +
  scale_fill_gradient2(low="cornflowerblue", mid = 'white', high="tomato", midpoint = 0.5)

coloc=within(coloc,{
  Proximity= 'NS'
  Proximity[ p.value < 0.05 | q.value < 0.05 ] = "High Proximity"
  Proximity[ p.value > 0.95 | q.value > 0.95 ] = "Low Proximity"
})

ggplot(coloc, aes(CellType1, CellType2, fill= Proximity)) + 
  geom_tile(color = 'black', size = 0.5) +
  scale_fill_manual(values = c('darkorange2', 'grey85')) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1))

# Calculate 3D distances from center ----
# center of 3D plot is 0,0,0. We can apply Pythagorean theorem with 3 dimensions to calculate 3D distance from each x,y,z coordinate as follows:
cellinfo_tbl$distance3D = sqrt(rowSums(
  (cellinfo_tbl[, c('x', 'y', 'z')] - 0)^2
))
# load library ggridges and ggplot2
library(ggridges)
cellinfo_tbl$CellType = factor(cellinfo_tbl$CellType, 
                               levels=c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                                        'Mac/Mono/cDC',
                                        'TILC_CCL5pos', 'ILC_KITpos',
                                        'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                                        'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg'))
ggplot(cellinfo_tbl, aes(x = distance3D, y = CellType, fill = CellType)) + 
  scale_fill_manual(values = alpha(c(rev(cols[1:15])), 0.6)) + 
  theme_bw() +
  geom_density_ridges(scale = 1.5, 
                      quantile_lines=TRUE, 
                      quantile_fun=function(distance3D,...)mean(distance3D))

# Calculate spherical compositions from center ----
#int <- max(cellinfo_tbl$distance3D)/5 # set an interval for radius to measure from central 0,0,0 coordinate so that we get 10 layers/intervals
int <- .01/5
cellinfo_tbl$layer3D <- with(
  cellinfo_tbl,
  as.factor(ifelse((distance3D <= (1*int)), 'Layer01', 
                   ifelse((distance3D <= 2*int) & (distance3D > (1*int)), 'Layer02',
                          ifelse((distance3D <= 3*int) & (distance3D > (2*int)), 'Layer03',
                                 ifelse((distance3D <= 4*int) & (distance3D > (3*int)), 'Layer04', 'Layer05'))))))

cellinfo_tbl$CellType = factor(cellinfo_tbl$CellType, 
                               levels=rev(c('B_AICDAneg_Cycling', 'B_AICDApos_Cycling', 'B_AICDApos', 'B_AICDAneg', 'ASC',
                                            'Mac/Mono/cDC',
                                            'TILC_CCL5pos', 'ILC_KITpos',
                                            'CD4Tfh', 'CD4T_KLF2neg', 'CD4T_KLF2pos',
                                            'gdT_CD2pos', 'CD8T_KLF2pos', 'CD8T_KLF2neg', 'gdT_CD2neg')))

ggplot(cellinfo_tbl, aes(x = layer3D, fill = CellType)) + 
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c(cols)) + 
  theme_classic()

newmetasub <- subset(cellinfo_tbl, CellType == 'CD4Tfh' | CellType == 'CD4T_KLF2neg' | CellType == 'CD4T_KLF2pos')
ggplot(newmetasub, aes(x = layer3D, fill = CellType)) + 
  geom_bar() +
  scale_fill_manual(values = c(cols[5:7])) + 
  theme_classic()

newmetasub <- subset(cellinfo_tbl, CellType == 'B_AICDAneg' | CellType == 'B_AICDAneg_Cycling' | CellType == 'B_AICDApos' | CellType == 'B_AICDApos_Cycling')
ggplot(newmetasub, aes(x = layer3D, fill = CellType)) + 
  geom_bar() +
  scale_fill_manual(values = c(cols[12:15])) + 
  theme_classic()

# Calculate spherical compositions from center using quantiles ----
cellinfo_tbl$quantile3D <- as.numeric(cut_number(cellinfo_tbl$distance3D, n = 5))
ggplot(cellinfo_tbl, aes(x = quantile3D, fill = CellType)) + 
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c(cols[1:15])) +
  theme_classic()

cellinfo_tbl$quantile3D <- as.character(cellinfo_tbl$quantile3D)
ggplot(cellinfo_tbl, aes(x = CellType, fill = quantile3D)) + 
  geom_bar(position = 'fill') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

newmetasub <- subset(cellinfo_tbl, CellType == 'CD4Tfh' | CellType == 'CD4T_KLF2neg' | CellType == 'CD4T_KLF2pos')
ggplot(newmetasub, aes(x = quantile3D, fill = CellType)) + 
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c(cols[5:7]))

newmetasub <- subset(cellinfo_tbl, CellType == 'B_AICDAneg' | CellType == 'B_AICDAneg_Cycling' | CellType == 'B_AICDApos' | CellType == 'B_AICDApos_Cycling')
ggplot(newmetasub, aes(x = quantile3D, fill = CellType)) + 
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c(cols[12:15]))

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.5 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggridges_0.5.4        plot3D_1.4            CSOmapR_1.0           writexl_1.4.2         viridis_0.6.4         viridisLite_0.4.2     ComplexHeatmap_2.12.1
#[8] scales_1.3.0          tidyr_1.3.1           uwot_0.1.16           Matrix_1.6-1          NMF_0.27              cluster_2.1.4         rngtools_1.5.2       
#[15] registry_0.5-1        future_1.33.1         readxl_1.4.3          patchwork_1.2.0       CellChat_2.1.2        ggplot2_3.5.0         igraph_1.5.1         
#[22] dplyr_1.1.4           SeuratDisk_0.0.0.9020 SeuratObject_4.1.3    Seurat_4.3.0.1       

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.4             spatstat.explore_3.2-3 reticulate_1.35.0      tidyselect_1.2.0       htmlwidgets_1.6.4      BiocParallel_1.30.4   
#[7] Rtsne_0.16             munsell_0.5.0          codetools_0.2-19       ica_1.0-3              miniUI_0.1.1.1         misc3d_0.9-1          
#[13] withr_3.0.0            spatstat.random_3.1-6  colorspace_2.1-0       progressr_0.14.0       Biobase_2.56.0         ggalluvial_0.12.5     
#[19] rstudioapi_0.15.0      stats4_4.2.2           ROCR_1.0-11            ggsignif_0.6.4         tensor_1.5             listenv_0.9.1         
#[25] labeling_0.4.3         polyclip_1.10-4        bit64_4.0.5            farver_2.1.1           coda_0.19-4.1          parallelly_1.37.1     
#[31] vctrs_0.6.5            generics_0.1.3         lambda.r_1.2.4         R6_2.5.1               doParallel_1.0.17      ggbeeswarm_0.7.2      
#[37] clue_0.3-64            RcppEigen_0.3.4.0.0    isoband_0.2.7          hdf5r_1.3.8            spatstat.utils_3.0-3   cachem_1.0.8          
#[43] promises_1.2.1         beeswarm_0.4.0         gtable_0.3.4           Cairo_1.6-1            globals_0.16.2         goftest_1.2-3         
#[49] rlang_1.1.3            systemfonts_1.0.5      GlobalOptions_0.1.2    splines_4.2.2          rstatix_0.7.2          lazyeval_0.2.2        
#[55] spatstat.geom_3.2-5    broom_1.0.5            BiocManager_1.30.22    yaml_2.3.8             reshape2_1.4.4         abind_1.4-5           
#[61] crosstalk_1.2.1        ggnetwork_0.5.13       backports_1.4.1        httpuv_1.6.14          tcltk_4.2.2            tools_4.2.2           
#[67] gridBase_0.4-7         statnet.common_4.9.0   ellipsis_0.3.2         jquerylib_0.1.4        RColorBrewer_1.1-3     BiocGenerics_0.42.0   
#[73] Rcpp_1.0.12            plyr_1.8.9             purrr_1.0.2            ggpubr_0.6.0           deldir_1.0-9           pbapply_1.7-2         
#[79] GetoptLong_1.0.5       cowplot_1.1.3          S4Vectors_0.34.0       zoo_1.8-12             ggrepel_0.9.5          magrittr_2.0.3        
#[85] data.table_1.15.2      RSpectra_0.16-1        futile.options_1.0.1   sna_2.7-2              scattermore_1.2        circlize_0.4.16       
#[91] lmtest_0.9-40          RANN_2.6.1             fitdistrplus_1.1-11    matrixStats_1.0.0      mime_0.12              xtable_1.8-4          
#[97] IRanges_2.30.1         gridExtra_2.3          shape_1.4.6.1          compiler_4.2.2         tibble_3.2.1           KernSmooth_2.23-22    
#[103] crayon_1.5.2           htmltools_0.5.7        later_1.3.2            formatR_1.14           MASS_7.3-60            car_3.1-2             
#[109] cli_3.6.2              parallel_4.2.2         pkgconfig_2.0.3        sp_2.0-0               plotly_4.10.4          spatstat.sparse_3.0-2 
#[115] foreach_1.5.2          svglite_2.1.3          vipor_0.4.5            bslib_0.6.1            stringr_1.5.1          digest_0.6.34         
#[121] sctransform_0.3.5      RcppAnnoy_0.0.21       spatstat.data_3.0-1    cellranger_1.1.0       leiden_0.4.3           shiny_1.8.0           
#[127] rjson_0.2.21           lifecycle_1.0.4        nlme_3.1-163           jsonlite_1.8.8         carData_3.0-5          network_1.18.2        
#[133] BiocNeighbors_1.14.0   fansi_1.0.6            pillar_1.9.0           lattice_0.21-8         ggrastr_1.0.2          fastmap_1.1.1         
#[139] httr_1.4.7             survival_3.5-7         glue_1.7.0             FNN_1.1.4              png_0.1-8              iterators_1.0.14      
#[145] bit_4.0.5              presto_1.0.0           stringi_1.8.3          sass_0.4.8             tidyverse_2.0.0        irlba_2.3.5.1         
#[151] future.apply_1.11.1    ape_5.7-1    
