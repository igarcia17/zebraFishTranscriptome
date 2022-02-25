library(readr)
library(tibble)
library(SummarizedExperiment)
library(SEtools)
#library(org.Dr.eg.db)
#library(GOSemSim)
library(edgeR)
library(HDF5Array)

#Locate files
samplesfile <- './data/sra_samples.csv'
countsfile <- './data/final_counts_info.csv'

#Load data
colData <- read_csv(samplesfile, show_col_types = FALSE)
counts <- read_csv(countsfile, skip_empty_rows = TRUE, show_col_types = FALSE)
names(counts)[1] <- 'Ensembl_ID'

#Remove duplicates and missing values of gene_name in counts variable
counts <- counts[!duplicated(counts$gene_name), ]
counts <- counts[!is.na(counts$gene_name),]

#Add a column to colData with a meaningful tag, make Samples column a factor
colData$colData_tag <- as.factor(apply(colData[c('Sample','Time','Treatment','Replicate', 'Study')], 1, 
                                       paste, collapse = "_"))
colData$Sample <- as.factor(colData$Sample)

#Get df of counts, normalized counts and row information. Clean colData.
countsData <- subset(counts, select= -c(Ensembl_ID))

colData <- colData[match(as.vector(colnames(countsData)), colData$colData_tag),]
colData <- colData[!is.na(colData$colData_tag),]

normalizedCounts <- cpm(Filter(is.numeric, countsData), normalized.lib.sizes=TRUE, 
                        log=TRUE, prior.count=1)
normalizedCounts <- cbind(subset(counts, select = c(gene_name)), normalizedCounts)

rowsInfo <- subset(counts, select = c(gene_name, Ensembl_ID))

#Set names of rows
colData <- column_to_rownames(colData, "colData_tag")
countsData <- column_to_rownames(countsData, 'gene_name')
normalizedCounts <- column_to_rownames(normalizedCounts, 'gene_name')
rowsInfo <- column_to_rownames(rowsInfo, 'gene_name')

#Create summarized Experiment
summExp <- SummarizedExperiment(assays=list(raw = countsData, 
                                            normalized = normalizedCounts), 
                                colData = colData)
rowData(summExp) <- rowsInfo

#Remove intermediate variable 'counts'
rm(counts)

saveHDF5SummarizedExperiment(summExp, dir="./datasummExp", prefix="", replace=FALSE,
                             chunkdim=NULL, level=NULL)
#_______________________________________________________________________________



