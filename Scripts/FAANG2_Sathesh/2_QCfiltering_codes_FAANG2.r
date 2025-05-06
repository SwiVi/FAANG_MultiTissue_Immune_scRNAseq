###########################################################
#### REMOVING NON_EXPRESSED GENES & POOR QUALITY CELLS ####
###########################################################
###### CONTINUATION OF 2_SoupX_codes.R. PRIOR OUTPUT FILES FROM SOUPX REQUIRED FOR EACH SAMPLE. Mito gene IDs also require.

pacman::p_load(dplyr, matrixStats, ggplot2, Rtsne, cowplot, Seurat, tidyr, DropletUtils)
setwd ("C:/Users/Sathesh.KumarSivasan/Documents/FAANG2/")

data_dir <- c(BM00 = "SoupX/BM00", BM98 = "SoupX/BM98", ICLN00 = "SoupX/ICLN00")
library_id <- c("BM00", "BM98", "ICLN00", "ICLN98")

lapply(data_dir, dir)
scRNA_data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = scRNA_data)

pDat <-data.frame(barcode = colnames(seurat_object))
pDat$SampleID <- gsub("([^_]+).+", "\\1", pDat$barcode, perl = TRUE)

pDat$BarBak <- pDat$barcode
pDat <- pDat %>% separate(BarBak, c("Sam","Loupe"))
pDat <- pDat[,-3] # remove Sam column
for (i in seq_along(library_id)){ 
  pDat$Loupe <- ifelse(pDat$SampleID == library_id[i], paste0(pDat$Loupe, paste0("-",i)), pDat$Loupe)
}
rownames(pDat) <- pDat$Loupe # make Loupe barcodes the rownames of pDat
head(pDat, n = 3); tail(pDat, n = 3); table(duplicated(pDat$Loupe))

table(pDat$SampleID)
#    BM00   BM98 ICLN00
#    2128   1448   9563

## Create count matrix (cDat)
cDat <-  as.matrix(GetAssayData(object = seurat_object, slot = 'counts'))
dim(cDat)	# 25880 Genes X 13139 Cells

## Filter out lowly expressed genes from the count matrix
keep <- rowSums(cDat) > 0
cDat <- cDat[keep,]
dim (cDat) # 17514 Genes X 13139 Cells

## Create a dataframe with feature data (fDat):
fDat <- data.frame(ID = rownames(cDat)) # create a dataframe of the filtered genes
rownames(fDat) <- fDat$ID # make genes the row names of dataframe
head(fDat, n = 3)
#                                                    ID
# TBP-ENSSSCG00000037372         TBP-ENSSSCG00000037372
# PSMB1-ENSSSCG00000027257     PSMB1-ENSSSCG00000027257
# FAM120B-ENSSSCG00000029697 FAM120B-ENSSSCG00000029697

# Add gene symbol information:
con <- gzfile(file.path(data_dir[1], "features.tsv.gz"))
ssc_genes <- read.delim(con, sep = "\t", header = FALSE, as.is = TRUE)[, 1:2]
colnames(ssc_genes) <- c("CellRanger", "Dashed")
ssc_genes$Dashed<- gsub("_", "-", ssc_genes$Dashed, perl = TRUE); head(ssc_genes, n = 3)

ssc_genes$Symbol <- sub("_.*", "", ssc_genes$CellRanger)
ssc_genes$EnsemblID <- sub(".*_", "", ssc_genes$CellRanger)
ssc_genes$Duplicated <- duplicated(ssc_genes$Symbol) | duplicated(ssc_genes$Symbol, fromLast = TRUE)
ssc_genes$Name <- ifelse(ssc_genes$Duplicated == "TRUE", ssc_genes$EnsemblID, ssc_genes$Symbol)
fDat <- merge(fDat, ssc_genes, by.x ="row.names", by.y = "Dashed", all.x =TRUE, sort = FALSE)
rownames(fDat) <- fDat[, 1]
fDat <- fDat[, -1]; head(fDat, n = 3)

mitoGenes <- read.table("MitoGenes.txt")
fDat$Mitochondrial <- fDat$EnsemblID %in% mitoGenes$V1
table(fDat$Mitochondrial); head(fDat, n = 3)
# FALSE  TRUE 
# 17479    35

## Modify count matrix (cDat) barcodes & gene names:
all(rownames(cDat) == rownames(fDat))# [1] TRUE
rownames(cDat) <- fDat$Name
head(rownames(cDat), n = 3) # [1] "TBP"     "PSMB1"   "FAM120B" 

all(colnames(cDat) == pDat$barcode) # [1] TRUE
colnames(cDat) <- pDat$Loupe 
tail(colnames(cDat), n = 3) # [1] "TTTGTTGGTATGGTAA-3" "TTTGTTGTCCACCCTA-3" "TTTGTTGTCTTGGCTC-3"

########################################
#### Filter out poor quality cells: ####
########################################
## Assess single-cell sequencing depths and number of genes detected:
pDat$UmiSums<- colSums(cDat)
pDat$GenesDetected <- colSums(cDat!=0)
ggplot(pDat, aes(x=SampleID, y=GenesDetected, fill= SampleID)) +
  geom_violin(draw_quantiles=0.5)+
  ylab("Total number of genes detected")

ggplot(pDat, aes(x=SampleID,y=UmiSums, fill=SampleID)) +
  geom_violin(draw_quantiles=0.5)+
  ylab("Total number of molecules(depth)")

## Plot 50 most highly expressed genes:
mRel <- t(t(cDat)/colSums(cDat)); table(colSums(mRel))

rownames(mRel)  <- fDat$ID
topExpressed <- rowMedians(mRel)
names(topExpressed) <- rownames(mRel)
topExpressed <- topExpressed %>% sort(.,decreasing=TRUE) %>% names
plotData <- t(mRel)[,topExpressed[1:50]] %>% reshape2::melt() %>%
  dplyr::rename(Cell=Var1, Gene=Var2, RelativeExpression=value)
plotData <- merge(plotData, fDat, by.x = "Gene", by.y = "ID", all.x = TRUE) # Add mitochondrial gene information for top 50 expressed genes
head(plotData, n = 3) # dataframe should look like this:

ggplot(plotData, aes(x=Gene, y=RelativeExpression, color= Mitochondrial)) +
  geom_boxplot() + coord_flip() + theme_bw()

## Plot 50 most frequently expressed genes:
freqOfExp <- cDat!=0
rownames(freqOfExp) <- fDat$ID
freqOfExp <- sort(rowSums(freqOfExp)/ncol(freqOfExp),decreasing=TRUE)
head(freqOfExp, n = 3)

plotData <- data.frame("Gene"=names(freqOfExp),"Frequency"=freqOfExp)
plotData <- merge(plotData, fDat, by.x = "Gene", by.y = "ID", all.x = TRUE, sorted =FALSE)
plotData <- plotData[order(plotData$Frequency, decreasing= TRUE), ]
head(plotData, n = 3)

## most frequently detected genes across all cells
ggplot(plotData[1:50,], aes(x=factor(Gene,levels=Gene), y=Frequency, color= Mitochondrial)) + geom_bar(stat="identity", fill ="white") +
  coord_flip() + xlab("Gene") + theme_bw()

## Assess cell quality:
theme_set(theme_grey())
mMito <- cDat[fDat$Mitochondrial,]
idtop <- fDat[fDat$Name %in% names(freqOfExp)[1:50],"ID"]
mTop <- cDat[idtop,]!=0
pDat$prcntTop <- colSums(mTop)/50
pDat$prcntMito <- colSums(mMito)/colSums(cDat)
ggplot(pDat, aes(x=prcntMito, y=GenesDetected, color = SampleID))+ geom_point() + 
  facet_wrap(~SampleID, nrow =1)+ theme_get() + ylab("#Genes detected per cell")

ggplot(pDat, aes(x=prcntTop, y=GenesDetected, color = SampleID))+ geom_point() + 
  facet_wrap(~SampleID, nrow =1) + theme_get()

## Check for barcode duplications:
pDat$DuplicatedBarcodes <- duplicated(rownames(pDat)) | duplicated(rownames(pDat), fromLast = TRUE)
table(pDat$DuplicatedBarcodes) # FALSE 13139

## Establish QC thresholds:
# Look at histograms to determine if there are obvious cutoff values:
ggplot(pDat, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) + ylim(0,30) + facet_wrap(~SampleID) +
  geom_vline(aes(xintercept=.10),color="red",lty="longdash") + # move this cutoff line where you see fit; 10.0% mitochondrial reads seems like a good cutoff
  RotatedAxis()

ggplot(pDat, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) + RotatedAxis() +
  geom_vline(aes(xintercept=300),color="red",lty="longdash") + # move this cutoff line where you see fit; 500 genes detected seems like a good cutoff
  facet_wrap(~SampleID) 

ggplot(pDat, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) + RotatedAxis() +
  geom_vline(aes(xintercept=500),color="red",lty="longdash") + # move this cutoff line where you see fit; 1000 UMIs seems like a good cutoff
  facet_wrap(~SampleID)

pDat <- mutate(pDat, 
    PassViability=prcntMito < 0.1, # only include cells with total mitochondrial reads under 10.0%
    PassGenesDet=GenesDetected > 300, # only consider cells with total genes detected more than 500
    PassLibSize=UmiSums > 500, # only consider cells with greater than 1000 total UMIs
    PassBarcodeFreq=DuplicatedBarcodes==FALSE, # only consider cells with non-repeated barcodes
    PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq) # list whether or not cells pass all filtering criteria
rownames(pDat) <- pDat$Loupe # make sure the pDat rownames correspond to Loupe barcodes!

## Overview over cells removed
table(pDat$SampleID,pDat$PassGenesDet)
         FALSE TRUE
  BM00     103 2025
  BM98     156 1292
  ICLN00    29 9534

table(pDat$SampleID,pDat$PassLibSize)
         FALSE TRUE
  BM00       0 2128
  BM98       0 1448
  ICLN00    13 9550

table(pDat$SampleID,pDat$PassViability)
         FALSE TRUE
  BM00     661 1467
  BM98     301 1147
  ICLN00   760 8803

table(pDat$SampleID,pDat$PassBarcodeFreq)
          TRUE
  BM00   2128
  BM98   1448
  ICLN00 9563

table(pDat$SampleID,pDat$PassAll)
          FALSE TRUE
  BM00     663 1465
  BM98     302 1146
  ICLN00   760 8803


## Save Data
stopifnot(identical(as.character(rownames(pDat)),colnames(cDat)))
stopifnot(identical(as.character(fDat$Name),rownames(cDat)))
out <- list()
out[["counts"]] <- cDat
out[["phenoData"]] <- pDat
out[["featureData"]] <- fDat
saveRDS(out,file=file.path("Filtered/", "FAANG2QC.rds")) # this saves all of our information before filtering out low quality cells

write.table(pDat, file = "PDAT_FAANG2_03312021.txt", sep = "\t", quote = FALSE)
write.table(fDat, file = "FDAT_FAANG2_03312021.txt", sep = "\t", quote = FALSE)

## Filter out poor quality cells:
dim(cDat) # 17514 Genes X 13139 Cells
cDat <- cDat[,pDat$PassAll]
dim(cDat) # 17514 Genes X 11414 Cells
pDat <- pDat[pDat$PassAll,] 
dim(pDat) # 11414    13

## subset data & save filtered data in CellRanger formats:
All <- CreateSeuratObject(counts = cDat, meta.data = pDat) # create Seurat object of counts & pheno data
Idents(All) <- "SampleID"

BM00 <- subset(All, ident = "BM00"); dim(BM00) # 17514  1465
write10xCounts(x = BM00@assays$RNA@counts, path = "Filtered/BM00", version = "3")
rm(BM00)

BM98 <- subset(All, ident = "BM98"); dim(BM98) # 17514  1146
write10xCounts(x = BM98@assays$RNA@counts, path = "Filtered/BM98", version = "3")
rm(BM98)

ICLN00 <- subset(All, ident = "ICLN00"); dim(ICLN00) # 17514  8803
write10xCounts(x = ICLN00@assays$RNA@counts, path = "Filtered/ICLN00", version = "3")
rm(ICLN00)

rm(All)
### END ###