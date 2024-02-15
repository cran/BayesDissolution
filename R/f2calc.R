#' Calculation of the f2 Statistic
#'
#' This function calculates the f2 statistic as described in Moore and Flanner (1996).
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#' the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#' For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#' measurements sorted in time. Alternatively, the user may provide a data object of class dis_data containing the
#' dissolution data. See the \code{make_dis_data()} function for the particular structure of the data object.
#' @return The function returns the f2 statistic calculated from the observed dissolution data.
#' @note Use the \code{plotdiss()} or \code{ggplotdiss()} function to visually check if it's appropriate to calculate the f2 statistic.
#' @references Moore, J.W. and Flanner, H.H. (1996). Mathematical comparison of distribution profiles. Pharmaceutical Technology, 20(6):64-74.
#' @examples
#' ### dis_data comes loaded with the package
#' f2calc(dis_data)
#'
#' @export
f2calc <- function(dis_data){

  if(is.data.frame(dis_data)){
    X <- cbind(2-(dis_data[,1] == dis_data[1,1]), dis_data[-1])
    names(X)[1] <- "Group"
    yRef  <- X[X$Group == 1,-1]
    yTest <- X[X$Group == 2,-1]
    dis_data <- make_dis_data(yRef, yTest)
  }
  if(class(dis_data)[1] != "dis_data"){
     stop("Input must be of class 'dis_data'. See function make_dis_data()")
  }else{
    D <- mean((dis_data$ybarR - dis_data$ybarT)^2)
    f2 <- 100 - 25*log10(1 + D)
    return(f2)
  }
}
