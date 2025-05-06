##################################################
#### REMOVING PROBABLE DOUBLETS FROM OUR DATA ####
##################################################

###### CONTINUATION OF doublet removal Python script
###### PRIOR OUTPUT .RDS FROM CELL/GENE FILTERING AND .CSV FROM SCRUBLET ANALYSIS REQUIRED AS INPUT FILES

pacman::p_load(dplyr, DropletUtils, Seurat)
setwd ("C:/Users/Sathesh.KumarSivasan/Documents/FAANG1/")

All <- readRDS("Filtered/FAANG1QC.rds")
dim(All$counts) # 18319 Genes X 52643 Cells

All$counts <- All$counts[,All$phenoData$PassAll] # keep only cells that passed all criteria in QC
dim(All$counts) # 18319 Genes X 40945 Cells

All$phenoData <- All$phenoData[All$phenoData$PassAll,]
dim(All$phenoData) # 40945    13

## Import our Scrublet scores:
BM00 <- read.csv("Scrublet/BM00_ScrubScore.csv")
BM98 <- read.csv("Scrublet/BM98_ScrubScore.csv")
ICLN00 <- read.csv("Scrublet/ICLN00_ScrubScore.csv")
ICLN98 <- read.csv("Scrublet/ICLN98_ScrubScore.csv")
SP00 <- read.csv("Scrublet/SP00_ScrubScore.csv")
SP98 <- read.csv("Scrublet/SP98_ScrubScore.csv")
THY00 <- read.csv("Scrublet/THY00_ScrubScore.csv")
THY98 <- read.csv("Scrublet/THY98_ScrubScore.csv")
scrub_combined <- rbind(BM00, BM98, ICLN00, ICLN98, SP00, SP98, THY00, THY98) 
dim(scrub_combined) # 40945     2
head(scrub_combined, n = 3)

All$phenoData$Scrublet <- scrub_combined$X0 # add column of scrublet data information for each cell to phenotype dataframe
All$phenoData <- mutate(All$phenoData, PassScrub = Scrublet < 0.25) # set the appropriate probability threshold... we decided on 0.25 for all samples here
rownames(All$phenoData) <- All$phenoData$Loupe # make sure the All$phenoData rownames correspond to Loupe barcodes!

table(All$phenoData$SampleID,All$phenoData$PassScrub) # how many cells passed scrubbing criteria?
         FALSE  TRUE
  BM00      20  1391
  BM98      29  2157
  ICLN00   165  4748
  ICLN98   239  6983
  SP00     177  3474
  SP98      95  2792
  THY00    543 11805
  THY98    192  6135

## Filter out probable doublets:
dim(All$counts) # 18319 Genes X 40945 Cells
All$counts <- All$counts[,All$phenoData$PassScrub] # keep only cells that passed all criteria in QC
dim(All$counts) # 18319 Genes X 39485 Cells

rownames(All$featureData) <- rownames(All$counts)
keep <- rowSums(All$counts) > 0
All$counts <- All$counts[keep,]
dim(All$counts) # 18062 Genes X 39485 Cells

All$phenoData <- All$phenoData[All$phenoData$PassScrub,]
dim(All$phenoData) # 39485    15

All$featureData <- All$featureData[rownames(All$counts),]
dim(All$featureData) # 18062     7

write.table(All$phenoData, file = "FAANG1_PDAT_FAANG1_03312021.txt", sep = "\t", quote = FALSE)

## Now save the data again in multiple formats:
stopifnot(identical(as.character(rownames(All$phenoData)),colnames(All$counts)))
stopifnot(identical(as.character(All$featureData$Name),rownames(All$counts)))
out <- list()
out[["counts"]] <- All$counts
out[["phenoData"]] <- All$phenoData
out[["featureData"]] <- All$featureData
saveRDS(out,file=file.path("QCed/", "FAANG1ScrubbedDFs.rds")) # this saves all of our information for counts, barcodes, and feature data as an .rds

seurat <- CreateSeuratObject(counts = All$counts, meta.data = All$phenoData) # create Seurat object of counts & pheno data
write10xCounts(x = seurat@assays$RNA@counts, path = "QCed/FAANG1ScrubbedCount", version = "3") # create CellRanger-like output files of our Seurat object
saveRDS(seurat,file=file.path("QCed/", "FAANG1ScrubbedSeurat.rds")) # this saves all of our information for counts, barcodes, and feature data as an .rds

## subset data & save filtered data in CellRanger formats:
Idents(seurat) <- "SampleID"

BM00 <- subset(seurat, ident = "BM00"); dim(BM00) # 18062  1391
write10xCounts(x = BM00@assays$RNA@counts, path = "QCed/BM00", version = "3")
rm(BM00)

BM98 <- subset(seurat, ident = "BM98"); dim(BM98) # 18062  2157
write10xCounts(x = BM98@assays$RNA@counts, path = "QCed/BM98", version = "3")
rm(BM98)

ICLN00 <- subset(seurat, ident = "ICLN00"); dim(ICLN00) # 18062  4748
write10xCounts(x = ICLN00@assays$RNA@counts, path = "QCed/ICLN00", version = "3")
rm(ICLN00)

ICLN98 <- subset(seurat, ident = "ICLN98"); dim(ICLN98) # 18062  6983
write10xCounts(x = ICLN98@assays$RNA@counts, path = "QCed/ICLN98", version = "3")
rm(ICLN98)

SP00 <- subset(seurat, ident = "SP00"); dim(SP00) # 18062  3474
write10xCounts(x = SP00@assays$RNA@counts, path = "QCed/SP00", version = "3")
rm(SP00)

SP98 <- subset(seurat, ident = "SP98"); dim(SP98) # 18062  2792
write10xCounts(x = SP98@assays$RNA@counts, path = "QCed/SP98", version = "3")
rm(SP98)

THY00 <- subset(seurat, ident = "THY00"); dim(THY00) # 18062  11805
write10xCounts(x = THY00@assays$RNA@counts, path = "QCed/THY00", version = "3")
rm(THY00)

THY98 <- subset(seurat, ident = "THY98"); dim(THY98) # 18062  6135
write10xCounts(x = THY98@assays$RNA@counts, path = "QCed/THY98", version = "3")
rm(THY98)

##### Woo-hoo! We've finished all of our QC analyses! Now we move on to the fun part!

