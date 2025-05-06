##################################################
#### REMOVING PROBABLE DOUBLETS FROM OUR DATA ####
##################################################

###### CONTINUATION OF doublet removal Python script
###### PRIOR OUTPUT .RDS FROM CELL/GENE FILTERING AND .CSV FROM SCRUBLET ANALYSIS REQUIRED AS INPUT FILES

pacman::p_load(dplyr, DropletUtils, Seurat)
setwd ("C:/Users/Sathesh.KumarSivasan/Documents/FAANG2/")

All <- readRDS("Filtered/FAANG2QC.rds")
dim(All$counts) # 17514 Genes X 13139 Cells

All$counts <- All$counts[,All$phenoData$PassAll] # keep only cells that passed all criteria in QC
dim(All$counts) # 17514 Genes X 11414 Cells

All$phenoData <- All$phenoData[All$phenoData$PassAll,]
dim(All$phenoData) # 11414    13

## Import our Scrublet scores:
BM00 <- read.csv("Scrublet/BM00_ScrubScore.csv")
BM98 <- read.csv("Scrublet/BM98_ScrubScore.csv")
ICLN00 <- read.csv("Scrublet/ICLN00_ScrubScore.csv")
scrub_combined <- rbind(BM00, BM98, ICLN00) 
dim(scrub_combined) # 11414     2
head(scrub_combined, n = 3)

All$phenoData$Scrublet <- scrub_combined$X0 # add column of scrublet data information for each cell to phenotype dataframe
All$phenoData <- mutate(All$phenoData, PassScrub = Scrublet < 0.25) # set the appropriate probability threshold... we decided on 0.25 for all samples here
rownames(All$phenoData) <- All$phenoData$Loupe # make sure the All$phenoData rownames correspond to Loupe barcodes!

table(All$phenoData$SampleID,All$phenoData$PassScrub) # how many cells passed scrubbing criteria?
         FALSE TRUE
  BM00      13 1452
  BM98       3 1143
  ICLN00   324 8479

## Filter out probable doublets:
dim(All$counts) # 17514 Genes X 11414 Cells
All$counts <- All$counts[,All$phenoData$PassScrub] # keep only cells that passed all criteria in QC
dim(All$counts) # 17514 Genes X 11074 Cells

rownames(All$featureData) <- rownames(All$counts)
keep <- rowSums(All$counts) > 0
All$counts <- All$counts[keep,]
dim(All$counts) # 17332 Genes X 11074 Cells

All$phenoData <- All$phenoData[All$phenoData$PassScrub,]
dim(All$phenoData) # 11074    15

All$featureData <- All$featureData[rownames(All$counts),]
dim(All$featureData) # 17332     7

write.table(All$phenoData, file = "FAANG2_PDAT_final_04012021.txt", sep = "\t", quote = FALSE)

## Now save the data again in multiple formats:
stopifnot(identical(as.character(rownames(All$phenoData)),colnames(All$counts)))
stopifnot(identical(as.character(All$featureData$Name),rownames(All$counts)))
out <- list()
out[["counts"]] <- All$counts
out[["phenoData"]] <- All$phenoData
out[["featureData"]] <- All$featureData
saveRDS(out,file=file.path("QCed/", "FAANG2ScrubbedDFs.rds")) # this saves all of our information for counts, barcodes, and feature data as an .rds

seurat <- CreateSeuratObject(counts = All$counts, meta.data = All$phenoData) # create Seurat object of counts & pheno data
write10xCounts(x = seurat@assays$RNA@counts, path = "QCed/FAANG2ScrubbedCount", version = "3") # create CellRanger-like output files of our Seurat object
saveRDS(seurat,file=file.path("QCed/", "FAANG2ScrubbedSeurat.rds")) # this saves all of our information for counts, barcodes, and feature data as an .rds

## subset data & save filtered data in CellRanger formats:
Idents(seurat) <- "SampleID"

BM00 <- subset(seurat, ident = "BM00"); dim(BM00) # 17332  1452
write10xCounts(x = BM00@assays$RNA@counts, path = "QCed/BM00", version = "3")
rm(BM00)

BM98 <- subset(seurat, ident = "BM98"); dim(BM98) # 17332  1143
write10xCounts(x = BM98@assays$RNA@counts, path = "QCed/BM98", version = "3")
rm(BM98)

ICLN00 <- subset(seurat, ident = "ICLN00"); dim(ICLN00) # 17332  8479
write10xCounts(x = ICLN00@assays$RNA@counts, path = "QCed/ICLN00", version = "3")
rm(ICLN00)


##### Woo-hoo! We've finished all of our QC analyses! Now we move on to the fun part!

