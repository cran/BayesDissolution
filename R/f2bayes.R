#' Calculation of a Bayesian 100*prob% credible interval for the F2 parameter
#'
#' This function calculates a 100*prob% credible interval for the F2 parameter using Bayesian methods. The model assumes a
#' version of the Jerffreys' prior with a pooled variance-covariance matrix from based on the reference and test data sets.
#' See Novick (2015) for more details of the model.
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#' the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#' For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#' measurements sorted in time. Alternatively, the user may provide a data object of class dis_data containing the
#' dissolution data. See the \code{make_dis_data()} function for the particular structure of the data object.
#' @param prob The probability associated with the credible interval. A value between 0 and 1.
#' @param B A positive integer specifying the number of Monte Carlo samples.
#' @param ci.type The type of credible interval to report. Specifying \code{quantile} returns a credible interval based on the posterior sample quantiles of the F2 distribution. Specifying \code{HPD} returns a highest posterior density interval.
#' @param get.dist logical; if \code{TRUE}, returns the posterior samples of the F2 distribution.
#' @return The function returns a 100*prob% credible interval for the F2 parameter calculated from the observed dissolution data.
#' @note Use the \code{plotdiss()} or \code{ggplotdiss()} function to visually check if it's appropriate to calculate the f2 statistic.
#' @references Novick, S., Shen, Y., Yang, H., Peterson, J., LeBlond, D., and Altan, S. (2015). Dissolution Curve Comparisons Through the F2 Parameter, a Bayesian Extension of the f2 Statistic. Journal of Biopharmaceutical Statistics, 25(2):351-371.
#' @references Pourmohamad, T., Oliva Aviles, C.M., and Richardson, R. Gaussian Process Modeling for Dissolution Curve Comparisons. Journal of the Royal Statistical Society, Series C, 71(2):331-351.
#' @examples
#' ### dis_data comes loaded with the package
#' f2bayes(dis_data, prob = 0.9, B = 1000)
#'
#' @export


f2bayes <- function(dis_data, prob = 0.9, B = 1000, ci.type = c("quantile", "HPD"), get.dist = FALSE){

  ## Bayesian version with Jeffreys' prior (pooled variance-covariance)
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
  }else if(prob <= 0 | prob >= 1 | !is.numeric(prob)){
    stop("prob must be a value between 0 and 1.")
  }else{

    K <- B
    V1 <- var(dis_data$yRef)
    V2 <- var(dis_data$yTest)
    Vhat <- 0.5 * (V1 + V2)         ## Pooled variance estimator

    degFree <- 2 * (dis_data$nTab - 1)

    Scale.hat <- solve(degFree * Vhat)
    Omega.post <- stats::rWishart(K, df = degFree, Sigma = Scale.hat )
    V.post <- lapply(1:K, function(k){ solve(Omega.post[,,k]) })
    Vhalf.post <- lapply(1:K, function(k){ solve(chol(Omega.post[,,k])) })

    muT.post <- matrix(rnorm(K * dis_data$nTime, sd = sqrt(1 / dis_data$nTab)), K, dis_data$nTime)
    muR.post <- matrix(rnorm(K * dis_data$nTime, sd=sqrt(1/dis_data$nTab)), K, dis_data$nTime)

    for (k in 1:K){
      muR.post[k,] <- dis_data$ybarR + Vhalf.post[[k]] %*% muR.post[k,]
      muT.post[k,] <- dis_data$ybarT + Vhalf.post[[k]] %*% muT.post[k,]
    }
    D.post <- rowMeans(sapply(1:dis_data$nTime, function(j){  (muT.post[,j] - muR.post[,j])^2 }))
    f2.post <- 100 - 25*log10(1 + D.post)
    out <- process_results("bayes", f2.post, ci.type, prob, get.dist)
  }
  return(out)
}
