library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(scales)

# Thymus-Thymus ----
## Load query data

#We will treat our scRNA-seq as the query dataset. 

#Load in h5 Seurat object of our data:

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/NormalizeIntegrateDimReducOutputs/PigImmuneMultiTissue_Thymus_CellTypeAnnotated.h5Seurat')

#Change default assay:

DefaultAssay(seu) <- 'SCT'

#How many cells, genes, gene assays, and dimensional reductions?

seu

#Read in reference dataset and update to Seurat v4:

ref <- readRDS('/home/Jayne.Wiarda/Thy_5999_rename.rds') # this file was graciously provided by authors of Gu et al.: doi.org/10.1016/j.celrep.2022.111050
ref # see we don't have SCT assay, so we create one now
ref <- SCTransform(ref, return.only.var.genes = FALSE, 
                   verbose = TRUE) 
DefaultAssay(ref) <- 'SCT'
ref$celltypes <- Idents(ref)
Idents(ref)
ref <- RenameIdents(object = ref, 
                    '0' = '0_DN_C',
                    '1' = '1_DN_Q',
                    '2' = '2_DP_C_1',
                    '3' = '3_DP_C_2',
                    '4' = '4_DP_Q',
                    '5' = '5_CD2pos_gd_T',
                    '6' = '6_CD2neg_gd_T',
                    '7' = '7_T_entry',
                    '8' = '8_Treg1',
                    '9' = '9_Treg2',
                    '10' = '10_CD8SP',
                    '11' = '11_CD4SP',
                    '12' = '12_ISG_CD8_T',
                    '13' = '13_cytotoxic_CD8_T',
                    '14' = '14_CD8aa',
                    '15' = '15_B_cells')
ref$celltypes <- Idents(ref)

cols <- c('#F8766D','#E68613','#CD9600','#ABA300','#7CAE00',
          '#0CB702','#D7EF16','#00C19A','#910C00','#ED68ED',
          '#00A9FF','#8494FF','#FFBDFF','#00B8E7', '#FF61CC','#FF68A1')
DimPlot(ref, cols = cols)

#Perform cell type predictions:
seu.list <- SplitObject(seu, split.by = "ID") # split by sample IDs
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

for(i in 1:length(seu.list)) {
  seu.list[[i]] <- RunPCA(seu.list[[i]], npcs = 50) # have to calculate and store PCA to run MappingScore() function later
}

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
                              refdata = list(cell_type = ref$celltypes), 
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
write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceMappingPrediction/ThymusToPigThymus_CellTypePredictions.xlsx')
write_xlsx(MappingScores, '/home/Jayne.Wiarda/scRNAseq_MultiTissueImmune/ReferenceMappingPrediction/ThymusToPigThymus_MappingScores.xlsx')

#Incorporate cell prediction & mapping scores into original Seurat object of query data:

DefaultAssay(seu) <- 'RNA'
seu <- AddMetaData(object = seu, 
                   metadata = c(CellTypePredictions, MappingScores))

#Visualize prediction & mapping results:

FeaturePlot(seu, reduction = 'tsne', features = 'MappingScores', min.cutoff = 0.5) & 
  scale_color_gradientn( colours = c('yellow', 'gold','orange', 'red')) # do we see high mapping scores, indicating good representation of query data by reference dataset? (Answer: yes)

VlnPlot(seu, split.by = 'CellType', features = 'MappingScores') # any major differences in mapping scores between ileum and jejunum data?
# Since jejunum and ileum samples have similar high mapping scores, suggests that the ileum atlas reference dataset annotations are cross-applicable to the jejunal samples used here
library(ggalluvial)


dat <- data.frame(seu$MappingScores)
colnames(dat) <- 'MappingScores'
dat <- dat %>% mutate(MappingScore = cut(MappingScores, 
                                            breaks=c(0,.2,.4,.6,.8,1)))
seu$MappingScoreBins <- dat$MappingScore
dat <- data.frame(table(seu$CellType, seu$predicted.id, seu$MappingScoreBins))
colnames(dat) <- c('CellType', 'PredictedID', 'MappingScore', 'Frequency')

ggplot(data = dat,
       aes(axis1 = CellType, axis2 = PredictedID, y = Frequency)) +
  geom_alluvium(aes(fill = MappingScore)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_brewer(type = "qual", palette = "YlOrRd")


DimPlot(seu, reduction = 'tsne', cols = cols, group.by = 'predicted.id')

# Show example of a few prediction scores:
FeaturePlot(seu, # plot all prediction scores
            features = c(colnames(CellTypePredictions[2:17])),
            reduction = 'tsne', # change to 'tsne' to overlay onto t-SNE plot
            ncol = 6) & 
  scale_color_gradientn( colours = c('grey90', 'blue3'),  limits = c(0, 1)) & 
  NoAxes() & NoLegend() 

VlnPlot(seu, split.by = 'CellType', features = 'prediction.score.max') # any major differences in prediction scores between ileum and jejunum data? 
VlnPlot(seu, group.by = 'predicted.id', features = 'prediction.score.max') # any major differences in prediction scores between ileum and jejunum data? 
# See mostly high prediction scores, giving confidence to cell type identities assigned via reference annotation-based prediction
DimPlot(seu, group.by = 'predicted.id', reduction = 'tsne', label = FALSE) # plot predicted cell IDs (based on highest prediction scores)

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

newdf$Annotation.x <- factor(newdf$Annotation.x, levels=c('Progenitor', 'ThymocyteDN_Cycling',
                                                          'ThymocyteDN', 'ThymocyteDP_Cycling',
                                                          'ThymocyteDP', 'gdT_CD2posCD8Apos',
                                                          'gdT_CD2posCD8Aneg', 'gdT_CD2pos_Cycling',
                                                          'gdT_CD2neg', 'gdT_CD2neg_Cycling',
                                                          'abT_Cycling', 
                                                          'abT_CD1DposCD1Epos', 'CD4T_KLF2neg_CCR9pos',
                                                          'CD4T_KLF2neg_CCR9neg', 'CD8T_KLF2neg',
                                                          'CD8T_KLF2pos', 'CD4T_KLF2pos',
                                                          'TILC_CCL5pos', 
                                                          'Bcell', 'ASC', 
                                                          'Mac/Mono/cDC', 'pDC'))
newdf$PredictedID.x <- factor(newdf$PredictedID.x, levels=c('0_DN_C',
                                                            '1_DN_Q',
                                                            '2_DP_C_1',
                                                            '3_DP_C_2',
                                                            '4_DP_Q',
                                                            '5_CD2pos_gd_T',
                                                            '6_CD2neg_gd_T',
                                                            '7_T_entry',
                                                            '8_Treg1',
                                                            '9_Treg2',
                                                            '10_CD8SP',
                                                            '11_CD4SP',
                                                            '12_ISG_CD8_T',
                                                            '13_cytotoxic_CD8_T',
                                                            '14_CD8aa',
                                                            '15_B_cells'))

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

seu$predicted.id <- factor(seu$predicted.id, levels=c('0_DN_C',
                                                                                '1_DN_Q',
                                                                                '2_DP_C_1',
                                                                                '3_DP_C_2',
                                                                                '4_DP_Q',
                                                                                '5_CD2pos_gd_T',
                                                                                '6_CD2neg_gd_T',
                                                                                '7_T_entry',
                                                                                '8_Treg1',
                                                                                '9_Treg2',
                                                                                '10_CD8SP',
                                                                                '11_CD4SP',
                                                                                '12_ISG_CD8_T',
                                                                                '13_cytotoxic_CD8_T',
                                                                                '14_CD8aa',
                                                                                '15_B_cells'))
cols <- c('#F8766D','#E68613','#CD9600','#ABA300','#7CAE00',
          '#0CB702','#D7EF16','#00C19A','#910C00','#ED68ED',
          '#00A9FF','#8494FF','#FFBDFF','#00B8E7', '#FF61CC','#FF68A1')
DimPlot(seu, group.by = 'predicted.id', cols = cols, reduction = 'tsne')
FeaturePlot(seu, reduction = 'tsne', 
            features = c('MappingScores'), 
            min.cutoff = 0.4) & 
  DarkTheme() &
  scale_color_gradientn( colours = c('white', 'yellow', 'gold','orange', 'red')) # do we see high mapping scores, indicating good representation of query data by reference dataset?

FeaturePlot(seu, features = c('prediction.score.0_DN_C', 'prediction.score.1_DN_Q',
                             'prediction.score.2_DP_C_1', 'prediction.score.3_DP_C_2',
                             'prediction.score.4_DP_Q', 'prediction.score.5_CD2pos_gd_T',
                             'prediction.score.6_CD2neg_gd_T', 'prediction.score.7_T_entry',
                             'prediction.score.8_Treg1', 'prediction.score.9_Treg2',
                             'prediction.score.10_CD8SP', 'prediction.score.11_CD4SP',
                             'prediction.score.12_ISG_CD8_T', 'prediction.score.13_cytotoxic_CD8_T', 
                             'prediction.score.14_CD8aa', 'prediction.score.15_B_cells'), 
            reduction = 'tsne', 
            ncol = 8, 
            cols = c('grey75', 'red3')) & NoAxes()

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