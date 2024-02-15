loadDataServer <- function(id) {
  moduleServer(id, function(input, output, session) {


    ## Load the data
    filedata <- reactive({
      req( input$input.file )

      extension = tools::file_ext((input$input.file$datapath))

      if ( extension=="csv" ){
        d = read.csv(input$input.file$datapath)
      }else if ( extension=="xlsx"){
        d = as.data.frame(readxl::read_excel(input$input.file$datapath))
      }


      X <- cbind(2 - (d[, 1] == d[1, 1]), d[-1])
      names(X)[1] <- "Group"
      yRef <- X[X$Group == 1, -1]
      yTest <- X[X$Group == 2, -1]
      dis_data <- make_dis_data(yRef, yTest)


      return(dis_data)

    })



    ## Download an example data set
    output$exampleData <- downloadHandler(
      filename = function() {
        "dis_data_example.xlsx"
      },
      content = function(file) {

        myfile <- paste0(system.file("extdata", package = "BayesDissolution"), "/dis_data_example.xlsx")
        file.copy(myfile, file)
      }
    )




    return(filedata)
  })
}

############################################################################
############################################################################

summaryServer <- function(id, filedata) {
  moduleServer(id, function(input, output, session) {

    ## Create scatterplot of allocation
    output$scatterplot <- renderPlot({

      req(input$input.file)

      dis_data = filedata()

      p = ggdissplot(dis_data, show.mean=TRUE, show.SD=TRUE)


      print(p)
    })

    output$results <- renderTable({
      req(input$input.file)

      dis_data = filedata()

      f2 = f2calc(dis_data)

      f2.ci = switch(input$ciMethod, "boot"={f2boot(dis_data, level=input$level, B=input$B)},
                     "boot BCA"={f2bca(dis_data, level=input$level, B=input$B)},
                     "PBS"={f2pbs(dis_data, level=input$level, B=input$B)},
                     "GPQ"={f2gpq(dis_data, level=input$level, B=input$B)},
                     "Bayes"={f2bayes(dis_data, prob=input$level, B=input$B)})
      out = data.frame(y=c(round(f2, 1),
                           paste(round(f2.ci$ci.quantile[1], 1), "-", round(f2.ci$ci.quantile[2], 1))))
      rownames(out) = c("Est", "95% CI")
      return(out)
    }, bordered=TRUE, colnames=FALSE, rownames=TRUE)

  })
}

############################################################################
############################################################################

server <- function(input, output, session) {

  filedata = loadDataServer("f2")

  summaryServer("f2", filedata)

}
