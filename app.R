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
library(pheatmap)
library(HDF5Array)
#_________________________________DATA__________________________________________
loadHDF5SummarizedExperiment(dir="./datasummExp", prefix="")
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

