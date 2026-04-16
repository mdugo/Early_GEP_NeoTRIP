
library(GEOquery)
library(here)
library(singscore)
library(GSEABase)

data_dir <- here("data", "processed")

#*******************************
# DOWNLOAD PROCESSED DATA FROM GEO
#*******************************
message("Downloading supplementary files...")

getGEOSuppFiles(
  GEO = "GSE319641",
  baseDir = data_dir,
  makeDirectory = FALSE
)

message("Download completed!")

#*******************************
# FILE DECOMPRESSION
#*******************************

files <- list.files(data_dir, full.names = TRUE, pattern = "TPM")
message("Unzipping: ", f)
R.utils::gunzip(files, overwrite = TRUE)

#*******************************
# DATA IMPORT
#*******************************

expData <- read.table(file.path(data_dir, "NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt"),
                    header = T,
                    sep = "\t",
                    as.is = T)
rownames(expData) <- expData$HGNC_approved_symbol
fdata <- expData[, c(1:3)]
expData <- expData[, -c(1:3)]

# !!! additonal sample information to be required to the corresponding author
clinical <- read.table(file.path(here("data"), "NeoTRIP_baseline_D1C2_pheno_and_clinical.txt"),
                       header = T,
                       sep = "\t",
                       as.is = T)

pheno <- data.frame(Sample_ID = colnames(expData),
                    Patient_ID = gsub("_.*", "", colnames(expData)),
                    Timepoint = gsub(".*_", "", colnames(expData)),
                    stringsAsFactors = F)
pheno <- merge(pheno, clinical[,-c(2,3)], by = "Sample_ID")
rownames(pheno) <- pheno$Sample_ID
pheno <- pheno[order(pheno$Patient_ID, pheno$Timepoint), ]
expData <- expData[,rownames(pheno)]
identical(rownames(pheno), colnames(expData))
pheno$Early.pCR <- factor(pheno$Early.pCR, levels = c("Sustained e-pCR", "Transient e-pCR", "pCR", "RD"))

# create gene sets list

gs <- read.table(file.path(here("data"), "geneSets_Supplementary_Table1.txt"),
                 header = T,
                 sep = "\t",
                 as.is = T)

gsList <- strsplit(gs$Genes, ";")
names(gsList) <- gs$GeneSet
gsAnnot <- gs[,c(1,2)]
rownames(gsAnnot) <- gsAnnot$GeneSet
gsAnnot$Annotation<-factor(gsAnnot$Annotation, levels = sort(unique(gsAnnot$Annotation))[c(1:7, 9:10, 8)])

gsListGS <- lapply(names(gsList), function(n) {
  GeneSet(gsList[[n]], setName = n)
})
names(gsListGS) <- names(gsList)

#*******************************
# RUN SINGSCORE
#*******************************

expData <- expData[!is.na(fdata$entrez_gene_id),]
isExp <- rowSums(expData>=2) >= ncol(expData)/100*10
expData <- expData[isExp,]
rankMat <- rankGenes(expData)
sScore <- multiScore(
  rankData = rankMat,
  upSetColc = gsListGS,
  centerScore = TRUE,
  knownDirection = TRUE
)
sScore <- sScore$Scores * 10

#*******************************
# ASSEMBLING EXPRESSION SET
#*******************************

scoreSS<-ExpressionSet(assayData = sScore,
                       phenoData=new("AnnotatedDataFrame",pheno),
                       featureData = new("AnnotatedDataFrame",gsAnnot))
save(scoreSS, file = file.path(here("data", "processed"), "NeoTRIP_singscore_BL_D1C2.RData"))




