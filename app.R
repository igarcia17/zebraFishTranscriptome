setwd("~/Desktop/inesGarcia/zebraFishTranscriptome")

library(shiny)
library(readr)
library(tibble)
library(SummarizedExperiment)
library(SEtools)
#library(org.Dr.eg.db)
#library(GOSemSim)
library(edgeR)
library(ggplot2)
library(pheatmap)
#Load data

##################is is really an advantage to have these two files at this point?
#####################isnt it better to load just one and make the colData from the headers?
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

#_________________________________________UI____________________________________
ui <- fluidPage(
  titlePanel(title = 
               'Gene expression in zebra fish across tissues and developmental stages'),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput(inputId = 'genes', choices = NULL,
                     label = 'Genes to work with:',
                     options = list(maxItems = 99),
                     multiple = TRUE
                     ),
      selectizeInput(inputId = 'tissue', choices = NULL,
                     label = 'Select tissue:',
                     multiple = FALSE
      )
 ),
    mainPanel(
      tabsetPanel(
        tabPanel("Counts", tableOutput(outputId = "rawCounts")),
        tabPanel("Heatmap", plotOutput(outputId = "heatmap")),
        tabPanel("Graph", plotOutput(outputId = "graph"))
    )
  )))
#rowData(summExp)$gene_name %in% input$genes

#______________________SERVER___________________________________________________
server <- function(input, output, session) {
  #summExp <- summExp[1:10, 1:10]
  updateSelectizeInput(session, 'genes', choices = rownames(summExp), 
                       selected = rownames(summExp)[sample(1:10, 5, replace=FALSE)],
                       server = TRUE)
  updateSelectizeInput(session, 'tissue', choices = summExp$Sample, 
                       selected = summExp$Sample[(sample(1:70, 2))],
                       server = TRUE)
  
  genelist <- reactive({input$genes})
  tissuelist <- reactive({input$tissue})
  
  dataTable <- reactive({assays(summExp)$raw[genelist(),]})
  dataHM <- reactive({assays(summExp)$normalized[genelist(),
                                                 colData(summExp)$Sample == tissuelist()]})
  
  output$rawCounts <- renderTable(dataTable(), rownames = TRUE, digits = 0)
  
  output$heatmap <- renderPlot(pheatmap(dataHM()))
  
  output$graph <- renderPlot(
    ggplot(data = dataHM(), aes(x = dataHM(), y = genelist()))
    )
    
}

shinyApp(ui = ui, server = server)

