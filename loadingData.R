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

#Get df of counts, normalized counts and row information
countsData <- subset(counts, select= -c(gene_name))
normalizedCounts <- cpm(Filter(is.numeric, countsData), normalized.lib.sizes=TRUE, 
                        log=TRUE, prior.count=1)
normalizedCounts <- cbind(subset(counts, select = c(Ensembl_ID)), normalizedCounts)

rowsData <- subset(counts, select = c(Ensembl_ID, gene_name))

#Set names of rows
colData <- column_to_rownames(colData, "SRR")
countsData <- column_to_rownames(countsData, "Ensembl_ID")
normalizedCounts <- column_to_rownames(normalizedCounts, 'Ensembl_ID')
rowsInfo <- column_to_rownames(rowsData, 'Ensembl_ID')

#Create summarized Experiment
summExp <- SummarizedExperiment(assays=list(raw = countsData, 
                                            normalized = normalizedCounts), 
                                colData = colData)
rowData(summExp) <- rowsInfo

assays(summExp)$normalized
se <- summExp[1:10, 1:10]
sehm(se, assayName="normalized", genes=rowData(summExp)$gene_name, do.scale=TRUE)
