#' Calculation of a parametric bootstrap 100*level% confidence interval for the F2 parameter
#'
#' This function calculates a 100*level% confidence interval for the F2 parameter using a parametric bootstrap based on a two variance component model with means for Time x Group, i.e.,  Dissolution ~ Time x Group + (1|Tablet:Group).
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#' the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#' For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#' measurements sorted in time. Alternatively, the user may provide a data object of class dis_data containing the
#' dissolution data. See the \code{make_dis_data()} function for the particular structure of the data object.
#' @param level The confidence level. A value between 0 and 1.
#' @param B A positive integer specifying the number of bootstrap samples.
#' @param ci.type The type of confidence interval to report. Specifying \code{quantile} returns a bootstrap confidence interval based on the sample quantiles. Specifying \code{HPD} returns a highest density region interval.
#' @param get.dist logical; if \code{TRUE}, returns the posterior samples of the F2 distribution.
#' @return The function returns a 100*level% confidence interval for the F2 parameter calculated from the observed dissolution data.
#' @note Use the \code{plotdiss()} or \code{ggplotdiss()} function to visually check if it's appropriate to calculate the f2 statistic.
#' @examples
#' ### dis_data comes loaded with the package
#' f2pbs(dis_data, level = 0.9, B = 1000)
#'
#' @export

f2pbs <- function(dis_data, level = 0.9, B = 1000, ci.type = c("quantile", "HPD"), get.dist = FALSE){

  if(is.data.frame(dis_data)){
    X <- cbind(2-(dis_data[,1] == dis_data[1,1]), dis_data[-1])
    names(X)[1] <- "Group"
    yRef  <- X[X$Group == 1,-1]
    yTest <- X[X$Group == 2,-1]
    dis_data <- make_dis_data(yRef, yTest)
  }
  if(class(dis_data)[1] != "dis_data"){
     stop("Input must be of class 'dis_data'. See function make_dis_data()")
  }else if(B <= 0 | !is.numeric(B)){
    stop("B must be a positive integer.")
  }else if(level <= 0 | level >= 1 | !is.numeric(level)){
    stop("level must be a value between 0 and 1.")
  }else{

    K <- B
    sigmaSqE <- dis_data$ssRes / dis_data$degFreeRes
    sigmaSqTab <- (dis_data$ssGroupTab / dis_data$degFreeGroupTab - sigmaSqE) / dis_data$nTime
    sigmaSqTab[sigmaSqTab < 0] <- 0
    sigmaSqTot <- sigmaSqTab + sigmaSqE

    Sigma <- matrix(sigmaSqTab, dis_data$nTime, dis_data$nTime)
    diag(Sigma) <- sigmaSqTot
    Sigma.diff <- 2 * Sigma / dis_data$nTab

    ybarDiff.pbs <- mnormt::rmnorm(K, mean = dis_data$ybarR - dis_data$ybarT, varcov = Sigma.diff)
    D.pbs <- rowMeans(ybarDiff.pbs^2)

    f2.pbs <- 100 - 25 * log10(1 + D.pbs)

    ## Get PBS 95% CI for f2 (there are other CI methods)
    out <- process_results("pbs", f2.pbs, ci.type, level, get.dist)
  }
  return(out)
}
