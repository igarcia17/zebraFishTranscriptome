# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# zebraFishTranscriptome - A Shiny app for visualizing the gene expression in 
# zebra fish animal model
# Created by Ines Garcia (@igarcia17), first version February 2022
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# electronic mail address: ines #dot# garciaortiz99 #at# gmail #dot# com
# 
library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(bslib)
library(readr)
library(tibble)
library(SummarizedExperiment)
library(edgeR)
library(pheatmap)


#_________________________________DATA__________________________________________
#Locate files
.samplesfile <- './data/sra_samples.csv'
.countsfile <- './data/final_counts_info.csv'

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

#_________________________________________UI____________________________________

ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "united"),
  titlePanel( h1("Gene expression in zebrafish", align = "center"), windowTitle = 
               'zebrafishTranscriptomics'),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput(inputId = 'genes', choices = NULL,
                     label = 'Genes to work with:',
                     options = list(maxItems = 99),
                     multiple = TRUE
                     ),
      selectizeInput(inputId = 'tissue', choices = NULL,
                     label = 'Select tissue:',
                     multiple = TRUE
      )
 ),
    mainPanel(
      tabsetPanel(
        tabPanel("Counts", tableOutput(outputId = "rawCounts")),
        tabPanel("Heatmap", plotOutput(outputId = "heatmap")),
        tabPanel("Disclaimer", textOutput(outputId = 'disc'))
      )
  )))

#______________________SERVER___________________________________________________
server <- function(input, output, session) {
  
  #UPDATES SELECTIZE
  updateSelectizeInput(session, 'genes', choices = rownames(summExp), 
                       selected = rownames(summExp)[sample(1:10, 5, replace=FALSE)],
                       server = TRUE)
  updateSelectizeInput(session, 'tissue', choices = levels(summExp$Sample), 
                       selected = levels(summExp$Sample),
                       server = TRUE)
  
  #VARIABLES
  genelist <- reactive({input$genes})
  tissuelist <- reactive({input$tissue})
  
  dataTable <- reactive({assays(summExp)$raw[genelist(),
                                             colData(summExp)$Sample %in% tissuelist()]})
 
  #OUTPUTS
  output$rawCounts <- renderTable(dataTable(), rownames = TRUE, digits = 0)
  
  output$heatmap <- renderPlot(
    pheatmap(assays(summExp)$normalized[genelist(),colData(summExp)$Sample %in% tissuelist()],
                                        scale="row",cluster_rows = FALSE,
                                        cluster_cols = TRUE))
  output$disc <- renderText('Disclaimer: if the heatmap does not appear until 
                            the window is resized, try to turn the device off from the 
                            R terminal')
  
}

shinyApp(ui = ui, server = server)

