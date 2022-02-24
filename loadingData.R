library(readr)
library(tibble)
library(SummarizedExperiment)
library(SEtools)
#library(org.Dr.eg.db)
#library(GOSemSim)
library(edgeR)
library(ggplot2)

#Load data
samplesfile <- './data/sra_samples.csv'
countsfile <- './data/final_counts.csv'

colData <- read_csv(samplesfile)
counts <- read_csv(countsfile, skip_empty_rows = TRUE)
names(counts)[1] <- 'Ensembl_ID'
#Remove duplicates and missing values of gene_name
counts <- counts[!duplicated(counts$gene_name), ]
counts <- counts[!is.na(counts$gene_name),]
#Add a column to colData with a meaningful tag
colData$colData_tag <- as.factor(apply(colData[c('Sample','Time','Treatment','Replicate', 'Study')], 1, 
                             paste, collapse = "_"))

#Get df of counts, normalized counts and row information
countsData <- subset(counts, select= -c(Ensembl_ID))
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

#_______________________________________________________________________________
#To get only a subset of certain class (for bargraph)
se1 <- summExp[,colData(summExp)$Sample == "Liver"]
se2 <- se1[colData(se1)$Time == '5dpf']
#For a particular gene
se3.2 <- se2[!is.na(rowData(se2)$gene_name == 'lrp6')]

assays(summExp)$normalized[rownames(summExp)[sample(1:10, 5, replace=FALSE)],tissueList()]

mat <- as.matrix(assays(se2)$normalized)
pheatmap(mat[1:10,1:20])
sechm(assay(se3.2)$normalized)
#Access data frame of normalized data
temp <- assays(summExp)$normalized
se <- summExp[1:10, 1:10]
sehm(se, assayName="normalized", genes=rowData(summExp)$gene_name, do.scale=TRUE)

tissuelist <- summExp$Sample[(sample(1:70, 2))]
genelist <- rownames(summExp)[sample(1:10, 5, replace=FALSE)]
trial <- assays(summExp)$normalized[genelist,
                           colData(summExp)$Sample == tissuelist]
       