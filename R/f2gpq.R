#' Calculation of a generalized pivotal quantity 100*level% confidence interval for the F2 parameter
#'
#' This function calculates a 100*level% confidence interval for the F2 parameter using generalized pivotal quantity methods based on a two variance component model with means for Time x Group, i.e.,  Dissolution ~ Time x Group + (1|Tablet:Group).
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#' the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#' For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#' measurements sorted in time. Alternatively, the user may provide a data object of class dis_data containing the
#' dissolution data. See the \code{make_dis_data()} function for the particular structure of the data object.
#' @param level The confidence level. A value between 0 and 1.
#' @param B The number of generalized pivotal quantity samples.
#' @param ci.type The type of confidence interval to report. Specifying \code{quantile} returns a bootstrap confidence interval based on the sample quantiles. Specifying \code{HPD} returns a highest density region interval.
#' @param get.dist logical; if \code{TRUE}, returns the posterior samples of the F2 distribution.
#' @return The function returns a 100*level% confidence interval for the F2 parameter calculated from the observed dissolution data.
#' @note Use the \code{plotdiss()} or \code{ggplotdiss()} function to visually check if it's appropriate to calculate the f2 statistic.
#' @examples
#' ### dis_data comes loaded with the package
#' f2gpq(dis_data, level = 0.9, B = 10000)
#'
#' @export
f2gpq <- function(dis_data, level = 0.9, B = 10000, ci.type = c("quantile", "HPD"), get.dist = FALSE){
 ## Generalized pivotal quantity analysis for 2 variance components problem
 ## Y[ijk] = theta[ij] + T[k(i)] + eps[ijk], T[k(i)] ~ N( 0, sigmaSqT ), eps[ijk] ~ N( 0, sigmaSqE )
 ##   i = 1,2 (groups); j=1, ..., P (time points), k = 1, 2, .., n (tablets).  Data must be balanced for this function!

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

    ybar_Diff <- dis_data$ybarR - dis_data$ybarT

    Z <- matrix(rnorm(B * dis_data$nTime), B, dis_data$nTime)
    U1 <- stats::rchisq(B, dis_data$degFreeGroupTab)
    U2 <- stats::rchisq(B, dis_data$degFreeRes)
    sigmaSqE <- dis_data$ssRes / U2
    sigmaSqT <- (dis_data$ssGroupTab / U1 - sigmaSqE) / dis_data$nTime
    sigmaSqT[sigmaSqT < 0] <- 0
    sigmaSqTot <- sigmaSqT + sigmaSqE

    muDiff <- matrix(NA, B, dis_data$nTime)
    for (j in 1:B){
      V <- matrix(sigmaSqT[j], dis_data$nTime, dis_data$nTime)
      diag(V) <- sigmaSqTot[j]
      Vchol <- chol(V)
      muDiff[j,] <- ybar_Diff - sqrt(2 / dis_data$nTab) * t(Vchol) %*% Z[j,]
    }

    D_gpq <- rowMeans(muDiff^2)
    f2.gpq <- 100 - 25 * log10(1 + D_gpq)
    out <- process_results("gpq", f2.gpq, ci.type, level, get.dist)
   }
   return(out)
}
