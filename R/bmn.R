#' Bayesian Multivariate Normal Model for Dissolution Profile Modeling
#'
#' This function implements the Bayesian multivariate normal model described in Pourmohamad et al (2022).
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#'  the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#'  For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#'  measurements sorted in time.
#' @param B A positive integer specifying the number of posterior samples to draw. By default \code{B} is set to 10000.
#' @return The function returns a list of B posterior samples for the following parameters:
#' \itemize{
#'   \item delta: A vector of posterior samples of delta as defined in Novick et. al 2015
#'   \item f2: A vector of posterior values of f2
#'   \item muR: A matrix of posterior samples for the reference group mean. Each row of the matrix corresponds to an observed time point, and each column of the matrix corresponds to a posterior sample.
#'   \item muT: A matrix of posterior samples for the test group mean. Each row of the matrix corresponds to an observed time point, and each column of the matrix corresponds to a posterior sample.
#' }
#' @note You should always check MCMC diagnostics on the posterior samples before drawing conclusions.
#' @references Novick, S., Shen, Y., Yang, H., Peterson, J., LeBlond, D., and Altan, S. (2015). Dissolution Curve Comparisons Through the F2 Parameter, a Bayesian Extension of the f2 Statistic. Journal of Biopharmaceutical Statistics, 25(2):351-371.
#' @references Pourmohamad, T., Oliva Aviles, C.M., and Richardson, R. Gaussian Process Modeling for Dissolution Curve Comparisons. Journal of the Royal Statistical Society, Series C, 71(2):331-351.
#' @examples
#' ### dis_data comes loaded with the package
#' ### We set B = 1000 to obtain 1000 posterior samples
#' B <- 1000
#' post <- bmn(dis_data, B = B)
#'
#' ### We can check how well the posterior samples of the means are mixing by
#' ### plotting the individual chains by time point
#' burnin <- B * 0.1     # A 10% burn-in
#' post$mu_R <- post$muR[,-(1:burnin)]
#' post$mu_T <- post$muT[,-(1:burnin)]
#'
#' N <- B - burnin      # Number of posterior samples after burn-in
#' chains <- data.frame(samples = rep(c(1:N, 1:N), each = ncol(dis_data) - 1),
#'                      group = rep(c("muR", "muT"), each = (ncol(dis_data) - 1) * N),
#'                      timepoint = paste("Time Point", rep(1:(ncol(dis_data) - 1), 2 * N)),
#'                      values = c(c(post$mu_R), c(post$mu_T)))
#'
#' g <- ggplot2::ggplot(chains, ggplot2::aes(samples, values)) +
#'                      ggplot2::geom_line() +
#'                      ggplot2::labs(x = "Iterations", y = "Posterior Sample Values") +
#'                      ggplot2::facet_wrap(group ~ timepoint) +
#'                      ggplot2::theme(text = ggplot2::element_text(size = 16))
#'
#' ### If we want to calculate the Pr(f2 > 50)
#' post$f2<- post$f2[-(1:burnin)]
#' prob <- sum(post$f2 > 50) / (B - burnin)
#'
#' ### Or if we want calculate a 95% credible interval for f2
#' alpha <- 0.05
#' f2_cred <- c(quantile(post$f2, alpha / 2),quantile(post$f2, 1 - alpha / 2))
#'
#' @export
bmn <- function(dis_data, B = 10000){

  if(B <= 0 | !is.numeric(B)){
    stop("B must be a positive integer.")
  }else if(!is.data.frame(dis_data)){
    stop("The dissolution data must be stored in a data frame.")
  }else if(length(unique(dis_data[,1])) != 2){
    stop("The dissolution data must contain 2 groups.")
  }else if(ncol(dis_data) < 3){
    stop("The dissolution data must contain at least 2 time points (but you probably intended for more than that).")
  }else{

    X <- cbind(2-(dis_data[,1] == dis_data[1,1]), dis_data[-1])
    names(X)[1] <- "Group"
    dat_R <- X[X$Group == 1,-1]
    dat_T <- X[X$Group == 2,-1]
    nlocs <- ncol(dat_R)
    nreps <- nrow(dat_R)

    Ybar_R <- apply(dat_R, 2, mean)
    Ybar_T <- apply(dat_T, 2, mean)

    anoms <- as.matrix(t(t(dat_R) - Ybar_R))
    S2_R <- matrix(0, nlocs, nlocs)
    for(i in 1:nreps){
      S2_R <- S2_R + (anoms[i,]) %*% t(anoms[i,])
    }

    anoms <- as.matrix(t(t(dat_T) - Ybar_T))
    S2_T <- matrix(0, nlocs, nlocs)
    for(i in 1:nreps){
      S2_T <- S2_T + (anoms[i,]) %*% t(anoms[i,])
    }

    nrun <- B
    f2 <- rep(NA, nrun)
    delta <- rep(NA, nrun)

    muR <- matrix(NA, nrow = nlocs, ncol = B)
    muT <- matrix(NA, nrow = nlocs, ncol = B)

    for(i in 1:nrun){

      V_R <- MCMCpack::riwish(nreps - 1, S2_R)
      mu_R <- mnormt::rmnorm(1, Ybar_R, V_R)

      V_T <- MCMCpack::riwish(nreps - 1, S2_T)
      mu_T <- mnormt::rmnorm(1, Ybar_T, V_T)

      delta[i] <- max(abs(mu_R - mu_T))
      f2[i] <- 50 * log10(1 / sqrt(1 + (1 / length(mu_R)) * sum((mu_R - mu_T)^2)) * 100)

      muR[,i] <- mu_R
      muT[,i] <- mu_T
    }

    return(list(delta = delta, f2 = f2, muR = muR, muT = muT))
  }
}
