#' Class dis_data creation
#'
#' This function creates a data object of class dis_data.
#'
#' @param yRef A data frame or matrix containing the dissolution data for the reference group data. The rows of the data set correspond
#' to the individual dissolution runs. The columns of the data frame contains the individual run's dissolution
#' measurements sorted in time.
#' @param yTest A data frame or matrix containing the dissolution data for the test group data. The rows of the data set correspond
#' to the individual dissolution runs. The columns of the data frame contains the individual run's dissolution
#' measurements sorted in time.
#' @return The function returns a data object of class dis_data.
#' @examples
#' ### dis_data comes loaded with the package
#' ### but need to update dis_data to be an object of class dis_data
#' new_dis_data <- make_dis_data(yRef = dis_data[dis_data$group == "Reference", -1],
#'                               yTest = dis_data[dis_data$group == "Test", -1])
#'
#' @export
make_dis_data <- function(yRef, yTest){

  out <- list()

  if(is.data.frame(yRef) == TRUE){
    yRef <- data.matrix(yRef)
  }
  if(is.data.frame(yTest) == TRUE){
    yTest <- data.matrix(yTest)
  }
  out$yRef <- yRef
  out$yTest <- yTest


  out$nTab <- nrow(yRef)
  out$nTime <- ncol(yRef)

  out$ybarR <- colMeans(yRef)
  out$ybarT <- colMeans(yTest)

  ## Sum of squares due to Group:Tablet
  out$ssGroupTab <- out$nTime * (sum((rowMeans(yRef) - mean(yRef))^2) + sum((rowMeans(yTest) - mean(yTest))^2))
  out$degFreeGroupTab <- 2 * out$nTab - 2

  ## Sum of squares residuals
  out$ssRes <- sum(sweep(yRef, 2, out$ybarR, FUN = "-")^2) +
        sum(sweep(yTest, 2, out$ybarT, FUN = "-")^2) - out$ssGroupTab
  out$degFreeRes <- 2 * (out$nTab - 1) * (out$nTime - 1)

  class(out) = c("dis_data", "list")
  return(out)
}

