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

#_________________________________________UI____________________________________
ui <- fluidPage(
  titlePanel(title = 
               'Gene expression in zebra fish across tissues and developmental stages'),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(inputId = 'genes',
                     label = 'Choose genes:',
                     choices = rownames(summExp),
                     selected = rownames(summExp)[sample(1:31330, 
                                                                  10, replace=FALSE)],
                     options = list(maxItems = 100),
                     multiple = TRUE
                     )
 ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Counts", tableOutput(outputId = "rawCounts"),
        tabPanel("Heatmap", plotOutput(outputId = "heatmap")),
        tabPanel("Graph", plotOutput(outputId = "graph"))
    )
  ))))
#rowData(summExp)$gene_name %in% input$genes

#______________________SERVER___________________________________________________
server <- function(input, output, session) {
  genelist.selected <- reactive({input$genes})
  
  dataTable <- assays(summExp)$raw[genelist.selected,]
  dataHM <- assays(summExp)$normalized[]
    
  output$rawCounts <- renderTable(dataTable, rownames = TRUE)
    
}

shinyApp(ui = ui, server = server)

