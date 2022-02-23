setwd("~/Desktop/inesGarcia/zebraFishTranscriptome")

library(shiny)
library(readr)
library(tibble)
library(SummarizedExperiment)
#library(org.Dr.eg.db)
#library(GOSemSim)
library(edgeR)
library(ggplot2)

#___________________________________DATA________________________________________
samplesfile <- './sra_samples.csv'
countsfile <- './final_counts.csv'

colData <- read_csv(samplesfile)
counts <- read_csv(countsfile, skip_empty_rows = TRUE)
names(counts)[1] <- 'Ensembl_ID'

#Add a column to colData with a meaningful tag
colData$colData_tag <- as.factor(apply(colData[c('Sample','Time','Treatment')], 1, 
                                       paste, collapse = "_"))

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

#Remove intermediate variable 'counts'
rm(counts)

#_________________________________________UI____________________________________
ui <- fluidPage(
  titlePanel(title = 
               'Gene expression in zebra fish across tissues and developmental stages'),
  sidebarLayout(
    sidebarPanel(
 ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Counts", tableOutput(outputId = "rawCounts"),
        tabPanel("Heatmap", plotOutput("heatmap")),
        tabPanel("Bar graph", plotOutput("bargraph"))
    )
  ))))

#______________________SERVER___________________________________________________
server <- function(input, output, session) {
    dataTable <- assays(summExp)$raw[]
    dataHM <- assays(summExp)$normalized[]
    output$rawCounts <- renderTable(datatohow, rownames = TRUE)
    
}

shinyApp(ui = ui, server = server)

