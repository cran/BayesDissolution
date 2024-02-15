
loadDataUI <- function(id) {

  ns <- NS(id)

  tagList(
    wellPanel(
      helpText(HTML("Load file. First column contains Group (Ref/Test). Each subsequent column",
                    "contains dissolution data for increasing time points.")),
      fileInput(ns("input.file"),
                label="Choose *.csv or *.xlsx File",
                accept=c("csv", ".csv", "xlsx", ".xlsx")),
      downloadLink(ns("exampleData"), label="Example Data"),

    )


  )


}

############################################################################
############################################################################

summaryUI <- function(id) {

  ns <- NS(id)

  tagList(

    fluidRow(column(3, selectInput(ns("ciMethod"), label="CI Method",
                                   choices=c("boot", "boot BCA", "PBS", "GPQ", "Bayes"),
                                   selected="boot")),
             column(3, numericInput(inputId=ns("level"), label="Level/Prob", value=0.9, min=0.5, max=1, step=0.01)),
             column(3, numericInput(inputId=ns("B"), label="# Samples", value=10000, min=1000, max=100000, step=1000))),
    fluidRow(column(4, tableOutput(ns("results"))),
             column(8, plotOutput(ns("scatterplot"))))

  )


}

############################################################################
############################################################################

ui <- navbarPage( title="Dissolution: F2",

                  tabPanel(title="Main Page",
                           fluidRow(
                             column(6, loadDataUI("f2")),
                           ),
                           summaryUI("f2")

                  )

)
