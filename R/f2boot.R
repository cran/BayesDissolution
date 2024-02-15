#' Calculation of a wild bootstrap 100*level% confidence interval for the F2 parameter
#'
#' This function calculates a 100*level% confidence interval for the F2 parameter using wild bootstrap
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
#' @references Liu, S. and Cai, X. and Shen, M. and Tsong, Y. (2023). In vitro dissolution profile comparison using bootstrap
#'  bias corrected similarity factor, f2. Journal of Biopharmaceutical Statistics, 34(1):78-89.
#' @examples
#' ### dis_data comes loaded with the package
#' f2boot(dis_data, level = 0.9, B = 1000)
#'
#' @export
f2boot <- function(dis_data, level = 0.9, B = 1000, ci.type = c("quantile", "HPD"), get.dist = FALSE){

  ## Wild bootstrap of f2 and bias-corrected f2
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

    ## Step 1.  Resample from the population by sampling with replacement from the nTab tablets of each group
    f2.bs <- rep(NA, K)
    for (k in 1:K){
      y1 <- dis_data$yRef[sample(1:dis_data$nTab, replace = TRUE),]
      y2 <- dis_data$yTest[sample(1:dis_data$nTab, replace = TRUE),]
      ybar1 <- colMeans(y1)
      ybar2 <- colMeans(y2)

      D <- mean((ybar1 - ybar2)^2)

      f2.bs[k] <- 100 - 25 * log10(1 + D)
    }

   ## Step 2.  Get 100*level% CI for f2
   out <- process_results(name = "bs", f2.dist = f2.bs, ci.type = ci.type, level = level, get.dist = get.dist)
  }

  return(out)
}
