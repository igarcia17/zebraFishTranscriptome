library(readr)
library(tibble)
library(SummarizedExperiment)
library(SEtools)
#library(org.Dr.eg.db)
#library(GOSemSim)
library(edgeR)
library(ggplot2)

#Load data
samplesfile <- './sra_samples.csv'
countsfile <- './final_counts.csv'

colData <- read_csv(samplesfile)
counts <- read_csv(countsfile, skip_empty_rows = TRUE)
names(counts)[1] <- 'Ensembl_ID'
#Remove duplicate genes and with missing value for gene_name
counts <- counts[!duplicated(counts$gene_name), ]
counts <- counts[!is.na(counts$gene_name),]
#Add a column to colData with a meaningful tag
colData$colData_tag <- as.factor(apply(colData[c('Sample','Time','Treatment')], 1, 
                             paste, collapse = "_"))

#uniqExp <- unique(colData$colData_tag)

#Get df of counts, normalized counts and row information
countsData <- subset(counts, select= -c(Ensembl_ID))
normalizedCounts <- cpm(Filter(is.numeric, countsData), normalized.lib.sizes=TRUE, 
                        log=TRUE, prior.count=1)
normalizedCounts <- cbind(subset(counts, select = c(gene_name)), normalizedCounts)

rowsData <- subset(counts, select = c(gene_name, Ensembl_ID))

#Set names of rows
colData <- column_to_rownames(colData, "SRR")
countsData <- column_to_rownames(countsData, 'gene_name')
normalizedCounts <- column_to_rownames(normalizedCounts, 'gene_name')
rowsInfo <- column_to_rownames(rowsData, 'gene_name')

#Create summarized Experiment
summExp <- SummarizedExperiment(assays=list(raw = countsData, 
                                            normalized = normalizedCounts), 
                                colData = colData)
rowData(summExp) <- rowsInfo

#Remove intermediate variable 'counts'
rm(counts)

trial <- aggSE(summExp, rowData(summExp)$gene_name)

#_______________________________________________________________________________
#To get only a subset of certain class (for bargraph)
se1 <- summExp[colData(summExp)$Sample == "Liver"]
se2 <- se1[colData(se1)$Time == '5dpf']
#For a particular gene
se3.2 <- se2[!is.na(rowData(se2)$gene_name == 'lrp6')]

mat <- as.matrix(assays(se2)$normalized)
pheatmap(mat[1:100,1:ncol(mat)/30])
sechm(se3.2)
#Access data frame of normalized data
temp <- assays(summExp)$normalized
se <- summExp[1:10, 1:10]
sehm(se, assayName="normalized", genes=rowData(summExp)$gene_name, do.scale=TRUE)
