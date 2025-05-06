###########################################################
#### REMOVING NON_EXPRESSED GENES & POOR QUALITY CELLS ####
###########################################################
###### CONTINUATION OF 2_SoupX_codes.R. PRIOR OUTPUT FILES FROM SOUPX REQUIRED FOR EACH SAMPLE. Mito gene IDs also require.

pacman::p_load(dplyr, scran, matrixStats, ggplot2, Rtsne, cowplot, Seurat, tidyr, knitr, DropletUtils)
setwd ("C:/Users/Sathesh.KumarSivasan/Documents/FAANG1/")

data_dir <- c(BM00 = "SoupX/BM00", BM98 = "SoupX/BM98", ICLN00 = "SoupX/ICLN00", ICLN98 = "SoupX/ICLN98", 
			SP00 = "SoupX/SP00", SP98 = "SoupX/SP98", THY00 = "SoupX/THY00", THY98 = "SoupX/THY98")
library_id <- c("BM00", "BM98", "ICLN00", "ICLN98", "SP00", "SP98", "THY00", "THY98")

lapply(data_dir, dir)
scRNA_data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = scRNA_data)

pDat <-data.frame(barcode = colnames(seurat_object))
pDat$SampleID <- gsub("([^_]+).+", "\\1", pDat$barcode, perl = TRUE)
dim(pDat) # 52643	2

pDat$BarBak <- pDat$barcode
pDat <- pDat %>% separate(BarBak, c("Sam","Loupe"))
pDat <- pDat[,-3] # remove Sam column
for (i in seq_along(library_id)){ 
  pDat$Loupe <- ifelse(pDat$SampleID == library_id[i], paste0(pDat$Loupe, paste0("-",i)), pDat$Loupe)
}
rownames(pDat) <- pDat$Loupe # make Loupe barcodes the rownames of pDat
head(pDat, n = 3); tail(pDat, n = 3); table(duplicated(pDat$Loupe))

table(pDat$SampleID)
#    BM00   BM98 ICLN00 ICLN98   SP00   SP98  THY00  THY98 
#    1444   2247   5069  11906   3737   3069  14961  10210

## Create count matrix (cDat)
cDat <-  as.matrix(GetAssayData(object = seurat_object, slot = 'counts'))
dim(cDat)	# 25880 Genes X 52643 Cells

## Filter out lowly expressed genes from the count matrix
keep <- rowSums(cDat) > 0
cDat <- cDat[keep,]
dim (cDat) #	18319 Genes X 52643 Cells

## Create a dataframe with feature data (fDat):
fDat <- data.frame(ID = rownames(cDat)) # create a dataframe of the filtered genes
rownames(fDat) <- fDat$ID # make genes the row names of dataframe
head(fDat, n = 3)

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
# 18284    35 

## Modify count matrix (cDat) barcodes & gene names:
all(rownames(cDat) == rownames(fDat))# [1] TRUE
rownames(cDat) <- fDat$Name
head(rownames(cDat), n = 3) # [1] "TBP"     "PSMB1"   "FAM120B" 

all(colnames(cDat) == pDat$barcode) # [1] TRUE
colnames(cDat) <- pDat$Loupe 
tail(colnames(cDat), n = 3) # [1] "TTTGTTGTCGCTAAAC-8" "TTTGTTGTCTAGTTCT-8" "TTTGTTGTCTGGTCAA-8"

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
table(pDat$DuplicatedBarcodes) # FALSE 52643


smryByGroup <- group_by(pDat, SampleID) %>%
    summarize(
        mPrcntMito = median(prcntMito),
        madPrcntMito = mad(prcntMito),
        threshold_PrcntMito = mPrcntMito + (2 * madPrcntMito),
        mGenesDetected = median(log10(GenesDetected)),
        madGenesDetected = mad(log10(GenesDetected)),
        threshold_GenesDetected = 10 ^ (mGenesDetected - (2 * madGenesDetected)),
        mUmiSums = median(log10(UmiSums)),
        madUmiSums = mad(log10(UmiSums)),
        threshold_UmiSums = 10 ^ (mUmiSums - (2 * madUmiSums))) %>%
    select(SampleID,starts_with("threshold")) 
kable(smryByGroup)
|SampleID | threshold_PrcntMito| threshold_GenesDetected| threshold_UmiSums|
|:--------|-------------------:|-----------------------:|-----------------:|
|BM00     |           0.0482118|                352.7713|          438.6381|
|BM98     |           0.0432851|                390.8827|          496.9991|
|ICLN00   |           0.0470202|                297.0174|          393.3472|
|ICLN98   |           0.1918148|                223.9471|          224.4235|
|SP00     |           0.0403860|                299.3544|          410.8096|
|SP98     |           0.0486567|                297.4626|          378.6679|
|THY00    |           0.1011791|                277.2525|          305.6773|
|THY98    |           0.1734174|                272.8829|          294.7839|

# 
|SampleID | threshold_PrcntMito| threshold_GenesDetected| threshold_UmiSums|
|:--------|-------------------:|-----------------------:|-----------------:|
|BM00     |           0.1	   |                400		|          500	   |
|BM98     |           0.1	   |                400		|          500	   |
|ICLN00   |           0.1	   |                300		|          400	   |
|ICLN98   |           0.15	   |                250		|          250	   |
|SP00     |           0.1	   |                300		|          500	   |
|SP98     |           0.1	   |                300		|          500	   |
|THY00    |           0.15	   |                300		|          500	   |
|THY98    |           0.15	   |                300		|          500	   |

### BM ### Mito - 0.1; Gene - 400; UMI - 500
BM00 <- pDat[which(pDat$SampleID == "BM00"),]
P1 <- ggplot(BM00, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) +
	ylim(0,30) + geom_vline(aes(xintercept=.1),color="red",lty="longdash") + RotatedAxis()
P2 <- ggplot(BM00, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) +
	geom_vline(aes(xintercept=400),color="red",lty="longdash") + RotatedAxis()
P3 <- ggplot(BM00, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) +
	geom_vline(aes(xintercept=500),color="red",lty="longdash") + RotatedAxis()

BM98 <- pDat[which(pDat$SampleID == "BM98"),]
P4 <- ggplot(BM98, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) +
	ylim(0,30) + geom_vline(aes(xintercept=.1),color="red",lty="longdash") + RotatedAxis()
P5 <- ggplot(BM98, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) +
	geom_vline(aes(xintercept=400),color="red",lty="longdash") + RotatedAxis()
P6 <- ggplot(BM98, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) +
	geom_vline(aes(xintercept=500),color="red",lty="longdash") + RotatedAxis()
plot_grid(P1, P2, P3, P4, P5, P6, nrow = 2)

BM <- pDat[which(pDat$SampleID == "BM00" | pDat$SampleID == "BM98"),]
BM <- mutate(BM, PassViability=prcntMito < 0.1, PassGenesDet=GenesDetected > 400, PassLibSize=UmiSums > 500, PassBarcodeFreq=DuplicatedBarcodes==FALSE,
    PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq)

P1 <- ggplot(BM, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) +
	ylim(0,30) + geom_vline(aes(xintercept=.1),color="red",lty="longdash") + facet_wrap(~SampleID) + RotatedAxis()
P2 <- ggplot(BM, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) +
	geom_vline(aes(xintercept=400),color="red",lty="longdash") + facet_wrap(~SampleID) + RotatedAxis()
P3 <- ggplot(BM, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) +
	geom_vline(aes(xintercept=500),color="red",lty="longdash") + facet_wrap(~SampleID) + RotatedAxis()

### ICLN ###
ICLN00 <- pDat[which(pDat$SampleID == "ICLN00"),] # Mito - 0.1; Gene - 300; UMI - 500
P1 <- ggplot(ICLN00, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) +
	ylim(0,30) + geom_vline(aes(xintercept=.1),color="red",lty="longdash") + RotatedAxis()
P2 <- ggplot(ICLN00, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) +
	geom_vline(aes(xintercept=300),color="red",lty="longdash") + RotatedAxis()
P3 <- ggplot(ICLN00, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) +
	geom_vline(aes(xintercept=500),color="red",lty="longdash") + RotatedAxis()

ICLN98 <- pDat[which(pDat$SampleID == "ICLN98"),] # Mito - 0.2; Gene - 300; UMI - 500
P4 <- ggplot(ICLN98, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) +
	ylim(0,30) + geom_vline(aes(xintercept=.15),color="red",lty="longdash") + RotatedAxis()
P5 <- ggplot(ICLN98, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) +
	geom_vline(aes(xintercept=300),color="red",lty="longdash") + RotatedAxis()
P6 <- ggplot(ICLN98, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) +
	geom_vline(aes(xintercept=500),color="red",lty="longdash") + RotatedAxis()
plot_grid(P1, P2, P3, P4, P5, P6, nrow = 2)

ICLN00 <- mutate(ICLN00, PassViability=prcntMito < 0.1, PassGenesDet=GenesDetected > 300, PassLibSize=UmiSums > 500, PassBarcodeFreq=DuplicatedBarcodes==FALSE,
    PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq)
ICLN98 <- mutate(ICLN98, PassViability=prcntMito < 0.15, PassGenesDet=GenesDetected > 300, PassLibSize=UmiSums > 500, PassBarcodeFreq=DuplicatedBarcodes==FALSE,
    PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq)

### SP ###
SP00 <- pDat[which(pDat$SampleID == "SP00"),] # Mito - 0.1; Gene - 300; UMI - 500
P1 <- ggplot(SP00, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) +
	ylim(0,30) + geom_vline(aes(xintercept=.1),color="red",lty="longdash") + RotatedAxis()
P2 <- ggplot(SP00, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) +
	geom_vline(aes(xintercept=300),color="red",lty="longdash") + RotatedAxis()
P3 <- ggplot(SP00, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) +
	geom_vline(aes(xintercept=500),color="red",lty="longdash") + RotatedAxis()

SP98 <- pDat[which(pDat$SampleID == "SP98"),] # Mito - 0.1; Gene - 300; UMI - 500
P4 <- ggplot(SP98, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) +
	ylim(0,30) + geom_vline(aes(xintercept=.1),color="red",lty="longdash") + RotatedAxis()
P5 <- ggplot(SP98, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) +
	geom_vline(aes(xintercept=300),color="red",lty="longdash") + RotatedAxis()
P6 <- ggplot(SP98, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) + scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) +
	geom_vline(aes(xintercept=500),color="red",lty="longdash") + RotatedAxis()
plot_grid(P1, P2, P3, P4, P5, P6, nrow = 2)

SP <- pDat[which(pDat$SampleID == "SP00" | pDat$SampleID == "SP98"),]
SP <- mutate(SP, PassViability=prcntMito < 0.1, PassGenesDet=GenesDetected > 300, PassLibSize=UmiSums > 500, PassBarcodeFreq=DuplicatedBarcodes==FALSE,
    PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq)


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
    PassGenesDet=GenesDetected > 300, # only consider cells with total genes detected more than 300
    PassLibSize=UmiSums > 500, # only consider cells with greater than 500 total UMIs
    PassBarcodeFreq=DuplicatedBarcodes==FALSE, # only consider cells with non-repeated barcodes
    PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq) # list whether or not cells pass all filtering criteria
rownames(pDat) <- pDat$Loupe # make sure the pDat rownames correspond to Loupe barcodes!

## Overview over cells removed
table(pDat$SampleID,pDat$PassGenesDet)
         FALSE  TRUE
  BM00       2  1442
  BM98       6  2241
  ICLN00     8  5061
  ICLN98    66 11840
  SP00      23  3714
  SP98      62  3007
  THY00    421 14540
  THY98    796  9414

table(pDat$SampleID,pDat$PassLibSize)
         FALSE  TRUE
  BM00       5  1439
  BM98      16  2231
  ICLN00   100  4969
  ICLN98   242 11664
  SP00      48  3689
  SP98      64  3005
  THY00    280 14681
  THY98    133 10077

table(pDat$SampleID,pDat$PassViability)
          FALSE  TRUE
  BM00      29  1415
  BM98      45  2202
  ICLN00    55  5014
  ICLN98  4655  7251
  SP00      18  3719
  SP98      75  2994
  THY00   2552 12409
  THY98   3882  6328

table(pDat$SampleID,pDat$PassBarcodeFreq)
          TRUE
  BM00    1444
  BM98    2247
  ICLN00  5069
  ICLN98 11906
  SP00    3737
  SP98    3069
  THY00  14961
  THY98  10210

table(pDat$SampleID,pDat$PassAll)
         FALSE  TRUE
  BM00      33  1411
  BM98      61  2186
  ICLN00   156  4913
  ICLN98  4684  7222
  SP00      86  3651
  SP98     182  2887
  THY00   2613 12348
  THY98   3883  6327


## Save Data
stopifnot(identical(as.character(rownames(pDat)),colnames(cDat)))
stopifnot(identical(as.character(fDat$Name),rownames(cDat)))
out <- list()
out[["counts"]] <- cDat
out[["phenoData"]] <- pDat
out[["featureData"]] <- fDat
saveRDS(out,file=file.path("Filtered/", "FAANG1QC.rds")) # this saves all of our information before filtering out low quality cells

## Filter out poor quality cells:
dim(cDat) # 18319 Genes X 52643 Cells
cDat <- cDat[,pDat$PassAll]
dim(cDat) # 18319 Genes X 40945 Cells
pDat <- pDat[pDat$PassAll,] 
dim(pDat) # 40945    13

## subset data & save filtered data in CellRanger formats:
All <- CreateSeuratObject(counts = cDat, meta.data = pDat) # create Seurat object of counts & pheno data
Idents(All) <- "SampleID"

BM00 <- subset(All, ident = "BM00"); dim(BM00) # 18319  1411
write10xCounts(x = BM00@assays$RNA@counts, path = "Filtered/BM00", version = "3")
rm(BM00)

BM98 <- subset(All, ident = "BM98"); dim(BM98) # 18319  2186
write10xCounts(x = BM98@assays$RNA@counts, path = "Filtered/BM98", version = "3")
rm(BM98)

ICLN00 <- subset(All, ident = "ICLN00"); dim(ICLN00) # 18319  4913
write10xCounts(x = ICLN00@assays$RNA@counts, path = "Filtered/ICLN00", version = "3")
rm(ICLN00)

ICLN98 <- subset(All, ident = "ICLN98"); dim(ICLN98) # 18319  7222
write10xCounts(x = ICLN98@assays$RNA@counts, path = "Filtered/ICLN98", version = "3")
rm(ICLN98)

SP00 <- subset(All, ident = "SP00"); dim(SP00) # 18319  3651
write10xCounts(x = SP00@assays$RNA@counts, path = "Filtered/SP00", version = "3")
rm(SP00)

SP98 <- subset(All, ident = "SP98"); dim(SP98) # 18319  2887
write10xCounts(x = SP98@assays$RNA@counts, path = "Filtered/SP98", version = "3")
rm(SP98)

THY00 <- subset(All, ident = "THY00"); dim(THY00) # 18319  12348
write10xCounts(x = THY00@assays$RNA@counts, path = "Filtered/THY00", version = "3")
rm(THY00)

THY98 <- subset(All, ident = "THY98"); dim(THY98) # 18319  6327
write10xCounts(x = THY98@assays$RNA@counts, path = "Filtered/THY98", version = "3")
rm(THY98)

rm(All)
### END ###