#' Run BayesDissolution Shiny App
#'
#' Runs a shiny application for calculating the different confidence and credible intervals for the F2 parameter. The different intervals are
#' constructed using the \code{f2bayes()}, \code{f2bca()}, \code{f2boot()}, \code{f2gpq()}, and \code{f2pbs()} functions. The shiny application comes
#' preloaded with an example excel data set based on the \code{dis_data} data set.
#'
#' @examples
#' ### The function requires no input to run
#' if(FALSE){ ## Make me TRUE to run
#'   runExample()
#' }
#'
#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "BayesDissolution")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `BayesDissolution`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

