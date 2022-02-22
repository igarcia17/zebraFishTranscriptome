
setwd("~/Desktop/inesGarcia/zebraFishTranscriptome")

#library(shiny)
library(readr)
library(tibble)
library(SummarizedExperiment)
#library(org.Dr.eg.db)
#library(GOSemSim)
#library(edgeR)
#library(mkAssay)

#Load data
#df <- read_csv('./final_counts_info.csv')
samplesfile <- './sra_samples.csv'
countsfile <- './final_counts.csv'

colData <- read_csv(samplesfile)
counts <- read_csv(countsfile, skip_empty_rows = TRUE)
names(counts)[1] <- 'Ensembl_ID'
#colnames(counts)

#Get df of counts and row information
countsData <- subset(counts, select= -c(gene_name))
rowsData <- subset(counts, select = c(Ensembl_ID, gene_name))

#Set names of rows
colData <- column_to_rownames(colData, "SRR")
countsData <- column_to_rownames(countsData, "Ensembl_ID")
rowsInfo <- column_to_rownames(rowsData, 'Ensembl_ID')

#Create summarized Experiment
summExp <- SummarizedExperiment(countsData, colData = colData)
rowData(summExp) <- rowsInfo

#Normalize matrix
librarySize <- colSums(Filter(is.numeric, summExp$rawcounts))


ui <- fluidPage(
  
)

server <- function(input, output, session) {
  
}

shinyApp(ui = ui, server = server)
'

