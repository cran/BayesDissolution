#' Hierarchical Gaussian Process Model for Dissolution Profile Modeling
#'
#' This function implements the Bayesian hierarchical Gaussian process model described in Pourmohamad et al (2022).
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#'  the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#'  For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#'  measurements sorted in time. Alternatively, the user may provide a data object of class dis_data containing the
#'  dissolution data. See the \code{make_dis_data()} function for the particular structure of the data object.
#' @param B A positive integer specifying the number of posterior samples to draw. By default \code{B} is set to 10000.
#' @param locs A vector in ascending order that corresponds to each time point the dissolution data was measured at.
#' @param B A positive integer specifying the number of posterior samples to draw. By default \code{B} is set to 10000.
#' @param n_interp An integer value specifying the number of time points to interpolate at. This sets the interploated points to be to \code{seq(1st time point, last time point, length = n_interp)}.
#' @param control An optional list of priors and initial values, otherwise the default values/strategies found in Pourmohamad et al (2022) will be used. More specifically, \code{control} can be used to define the following settings:
#' \itemize{
#'    \item \code{sigma2_starting}: starting value for sigma^2
#'    \item \code{tau2_starting}: starting value for tau^2
#'    \item \code{phi_starting}: starting value for phi
#'    \item \code{psi_starting}: starting value for psi
#'    \item \code{sigma2_alpha} and \code{sigma2_beta}: parameters for the inverse gamma prior for sigma^2
#'    \item \code{tau2_alpha} and \code{tau2_beta}: parameters for the inverse gamma prior for tau^2
#'    \item \code{phi_alpha} and \code{phi_beta}: parameters for the gamma prior for phi
#'    \item \code{psi_alpha} and \code{psi_beta}: parameters for the gamma prior for psi
#'    \item \code{prop_phi}: proposal variance for the parameter phi
#'    \item \code{prop_psi}: proposal variance for the parameter psi
#' }
#' @param adaptive logical; an option for using adaptive MCMC. If \code{adaptive = TRUE}, this will replace both \code{prop_phi} and \code{prop_psi} by using past MCMC draws to inform the proposal variance.
#' @return The function returns a list of summary statistics and B posterior samples for parameters of the model. More specifically it returns:
#' \itemize{
#'   \item delta: The average delta value over the posterior samples of delta. The definition of delta is given in Novick et. al 2015.
#'   \item f2: The average f2 value over the posterior samples of f2.
#'   \item mcmc_chains: A list of posterior samples for delta, f2, the mean parameters (\code{mu_pars}), and the covariance parameters (\code{cov_pars}).
#' }
#' @note You should always check MCMC diagnostics on the posterior samples before drawing conclusions. Likewise, plots of the predicted dissolution curves should also be checked to evaluate if the model fit to the observed data seems reasonable.
#' @references Novick, S., Shen, Y., Yang, H., Peterson, J., LeBlond, D., and Altan, S. (2015). Dissolution Curve Comparisons Through the F2 Parameter, a Bayesian Extension of the f2 Statistic. Journal of Biopharmaceutical Statistics, 25(2):351-371.
#' @references Pourmohamad, T., Oliva Aviles, C.M., and Richardson, R. Gaussian Process Modeling for Dissolution Curve Comparisons. Journal of the Royal Statistical Society, Series C, 71(2):331-351.
#' @examples
#' ### dis_data comes loaded with the package
#' ### We set B = 100 to obtain 100 posterior samples, you probably want to run it
#' ### longer for say, B = 100000, but B = 100 runs fast for illustrative purposes
#' ### and passing CRAN checks
#' B <- 100
#'
#' tp <- seq(10, 80, 10) # Time points
#' burnin <- B * 0.1     # A 10% burn-in
#' thin <- 10            # Keep every 10th sample, i.e., thinning
#' post <- hgp(dis_data, locs = tp, B = B, n_interp = 100)
#'
#' ### Example: Removing burn-in and then thinning the posterior samples for the covariance parameters
#' ###          and then plotting the chains
#' phi <- post$mcmc_chains$cov_pars$phi[-c(1:burnin)]
#' phi <- phi[seq(1, (B-burnin), thin)]
#' psi <- post$mcmc_chains$cov_pars$psi[-c(1:burnin)]
#' psi <- psi[seq(1, (B-burnin), thin)]
#' sigma_R <- post$mcmc_chains$cov_pars$sigma_R[-c(1:burnin)]
#' sigma_R <- sigma_R[seq(1, (B-burnin), thin)]
#' sigma_T <- post$mcmc_chains$cov_pars$sigma_T[-c(1:burnin)]
#' sigma_T <- sigma_T[seq(1, (B-burnin), thin)]
#' tau <- post$mcmc_chains$cov_pars$tau[-c(1:burnin)]
#' tau <- tau[seq(1, (B-burnin), thin)]
#'
#' chains <- data.frame( # Data frame holding posterior samples
#' samples = rep(1:((B-burnin)/thin), times = 5),
#' parameter = rep(c("phi", "psi", "sigma_R", "sigma_T", "tau"),
#'                 each = (B-burnin)/thin),
#' values = c(phi, psi, sigma_R, sigma_T, tau))
#' chains$parameter <- factor(chains$parameter,
#'                            labels = c(expression(phi),
#'                                       expression(psi),
#'                                       expression(sigma[R]),
#'                                       expression(sigma[T]),
#'                                       expression(tau)))
#' ggplot2::ggplot(chains, ggplot2::aes(samples, values)) +
#'  ggplot2::geom_line() +
#'  ggplot2::labs(x = "Iterations", y = "Posterior Sample Values") +
#'  ggplot2::facet_wrap(~parameter, scales = "free",
#'              labeller = ggplot2::label_parsed) +
#'  ggplot2::theme(text = ggplot2::element_text(size = 22))
#'
#' ggplot2::ggplot(chains, ggplot2::aes(values)) +
#'  ggplot2::geom_density() +
#'  ggplot2::labs(x = "Values", y = "Posterior Density") +
#'  ggplot2::facet_wrap(~parameter, scales = "free",
#'             labeller = ggplot2::label_parsed) +
#'  ggplot2::theme(text = ggplot2::element_text(size = 22))
#'
#' ### Plotting the predicted dissolution profiles
#' dissplot(dis_data, tp)
#' grid <- sort(c(tp, seq(min(tp), max(tp), length = 100)))
#' grid1 <- (1:B)[-(1:burnin)][seq(1, (B-burnin), thin)]
#' grid2 <- ((B+1):(2*B))[-(1:burnin)][seq(1, (B-burnin), thin)]
#' lines(grid, apply(post$mcmc_chains$mu_pars[,grid1], 1, mean),
#'       col = "gray65", lwd = 2)
#' lines(grid, apply(post$mcmc_chains$mu_pars[,grid2], 1, mean),
#'       col = "black", lwd = 2)
#' lower <- apply(post$mcmc_chains$mu_pars[,grid1], 1,
#'                quantile, prob = 0.025)
#' upper <- apply(post$mcmc_chains$mu_pars[,grid1], 1,
#'                quantile, prob = 0.975)
#' polygon(c(grid, rev(grid)), c(lower, rev(upper)),
#'         col = scales::alpha("gray65",.2), border = NA)
#' lower <- apply(post$mcmc_chains$mu_pars[,grid2], 1,
#'               quantile, prob = 0.025)
#' upper <- apply(post$mcmc_chains$mu_pars[,grid2], 1,
#'                quantile, prob = 0.975)
#' polygon(c(grid, rev(grid)), c(lower, rev(upper)),
#'         col = scales::alpha("black",.2), border = NA)
#'
#' ### If we want to calculate the Pr(f2 > 50 & delta < 15)
#' prob <- sum(post$mcmc_chains$f2[grid1] > 50 &
#'             post$mcmc_chains$delta[grid1] < 15) / ((B - burnin)/thin)
#'
#' @importFrom stats dgamma rnorm runif sd var
#' @export
hgp <- function(dis_data, locs, B = 1000, n_interp = 30, control = list(), adaptive = FALSE){

  if(class(dis_data)[1] == "dis_data"){
    dis_data <- data.frame(rbind(data.frame(group = "Reference", dis_data$yRef), data.frame(group = "Test", dis_data$yTest)))
  }
  if(B <= 0 | !is.numeric(B)){
    stop("B must be a positive integer.")
  }else if(!is.data.frame(dis_data)){
    stop("The dissolution data must be stored in a data frame.")
  }else if(length(unique(dis_data[,1])) != 2){
    stop("The dissolution data must contain 2 groups.")
  }else if(ncol(dis_data) < 3){
    stop("The dissolution data must contain at least 2 time points (but you probably intended for more than that).")
  }else if(length(locs) != (ncol(dis_data) - 1)){
    stop("The number of time points does not match the dissolution data set.")
  }else if(!is.numeric(locs)){
    stop("The time point vector needs to be a numeric vector of length ncol(dis_data) - 1.")
  }else if(n_interp <= 0 | !is.numeric(n_interp)){
    stop("n_interp must be a positive integer.")
  }else if(!is.list(control)){
    stop("control needs to be a list of priors and/or initial values")
  }else{

    draw_mu = function(dat,n,Sigma){
      Xbar <- apply(dat,2,mean)
      Sstar <- solve(solve(S)+n*solve(Sigma))
      mstar <- Sstar%*%(solve(S)%*%m+n*solve(Sigma)%*%(Xbar))
      mnormt::rmnorm(1,mstar,Sstar)
    }

    draw_phiS <- function(){
      alpha <- ifelse(is.null(control$phiS_alpha),(max(locs)-min(locs))/(nloc/2),control$phiS_alpha)
      beta <- ifelse(is.null(control$phiS_beta),1,control$phiS_beta)
      phiSnew <- exp(rnorm(1,log(phiS),prop_phiS))
      Snew <- tau2S*covf2(D,phiSnew)+diag(nloc)*.0001
      alph <- mnormt::dmnorm(mu_A,m,Snew,log=TRUE)+mnormt::dmnorm(mu_B,m,Snew,log=TRUE)+dgamma(phiSnew,alpha,beta,log=TRUE)-
        mnormt::dmnorm(mu_A,m,S,log=TRUE)-mnormt::dmnorm(mu_B,m,S,log=TRUE)-dgamma(phiS,alpha,beta,log=TRUE)
      if(log(runif(1)) < alph){
        phiS <- phiSnew
        S <- Snew
        acceptS <- acceptS + 1
      }
      phiS
    }

    draw_phi <- function(){
      phinew <- exp(rnorm(1,log(phi),prop_phi))
      Sigma_Anew <- tau2_A*covdelta(D,phinew)
      Sigma_Bnew <- tau2_B*covdelta(D,phinew)

      alpha <- ifelse(is.null(control$phi_alpha),(max(locs)-min(locs))/(nloc/2),control$phi_alpha)
      beta <- ifelse(is.null(control$phi_beta),1,control$phi_beta)

      alph <- sum(mnormt::dmnorm(dat_A,mu_A,Sigma_Anew,log=TRUE))+sum(mnormt::dmnorm(dat_B,mu_B,Sigma_Bnew,log=TRUE))+dgamma(phinew,alpha,beta,log=TRUE)-
        sum(mnormt::dmnorm(dat_A,mu_A,Sigma_A,log=TRUE)+mnormt::dmnorm(dat_B,mu_B,Sigma_B,log=TRUE))-dgamma(phi,alpha,beta,log=TRUE)
      if(log(runif(1)) < alph){
        phi <- phinew
        Sigma_A <- Sigma_Anew
        Sigma_B <- Sigma_Bnew
        accept <- accept + 1
      }
      phi
    }

    draw_tau2S <- function(){
      alpha <- ifelse(is.null(control$tau2S_alpha),mean(diag(var(X[,-1]))),control$tau2S_alpha)
      beta <- ifelse(is.null(control$tau2S_beta),mean(diag(var(X[,-1]))),control$tau2S_beta)
      pscl::rigamma(1,alpha+nloc*2/2,beta+1/2*(mu_A-m)%*%solve(covf2(D,phiS),mu_A-m)+
                1/2*(mu_B-m)%*%solve(covf2(D,phiS),mu_B-m))
    }

    draw_tau2<- function(dat,mu,n){
      alpha <- ifelse(is.null(control$tau2_alpha),mean(diag(var(X[,-1]))),control$tau2_alpha)
      beta <- ifelse(is.null(control$tau2_beta),mean(diag(var(X[,-1]))),control$tau2_beta)
      wsse <- 0
      anom <- as.matrix(t(t(dat)-mu))
      wsse <- sum(diag(anom%*%solve(covdelta(D,phi))%*%t(anom)))
      pscl::rigamma(1,alpha+(nloc*n)/2,beta+1/2*wsse)
    }

    mu_post <- function(mu,tau2){
      S22 <- tau2S*covf2(D22,phiS)
      S12 <- tau2S*covf2(D12,phiS)
      mnew <- 50+t(S12)%*%solve(S,mu-50)
      signew <- S22 - t(S12)%*%solve(S,(S12))
      munew <- mnormt::rmnorm(1,mnew,signew)
      c(mu,munew)
    }

    # Initialize Data
    X <- cbind(2-(dis_data[,1]==dis_data[1,1]),dis_data[-1])
    names(X)[1] <- "Group"
    nloc <- length(locs)
    otherlocs <- seq(min(locs),max(locs),length=n_interp)
    alllocs <- c(locs,otherlocs)
    nreps_A <- sum(X[,1]==1)
    nreps_B <- sum(X[,1]==2)
    D <- abs(outer(locs,locs,"-"))
    D22 <- abs(outer(otherlocs,otherlocs,"-"))
    D12 <- abs(outer(locs,otherlocs,"-"))
    dat_A <- X[X$Group==1,-1]
    dat_B <- X[X$Group==2,-1]
    dat <- X

    # Initialize mu, sigma2,tau2,phi,tau2S
    # 1st group is denoted with _A and second group with _B
    mu_A <- rep(0,nloc)
    mu_B <- rep(0,nloc)
    tau2 <- ifelse(is.null(control$tau2_starting),mean(diag(var(X[,-1]))),control$tau2_starting)
    tau2_A <- tau2 # Covariance of 1st group mean
    tau2_B <- tau2 #Covariance of 2nd group mean
    phi <- ifelse(is.null(control$phi_starting),(max(locs)-min(locs))/(nloc/2),control$phi_starting)
    phiS <- phi
    tau2S <- ifelse(is.null(control$tau2S_starting),tau2*5,control$tau2S_starting) # Covariance of mean function
    covdelta <- function(D,phi)geoR::matern(D,phi,3/2)
    Sigma_A <- tau2_A*covdelta(D,phi) # 1st group
    Sigma_B <- tau2_B*covdelta(D,phi) # second group
    m <- rep(50,nloc)
    covf2 <- function(D,phiS)exp(-abs(D/phiS))
    S <- tau2S*covf2(D,phiS)
    # Metropolis Hastings proposal parameters
    prop_phiS <- ifelse(is.null(control$phiS_prop),1,control$phiS_prop)
    prop_phi <- ifelse(is.null(control$phi_prop),1,control$phi_prop)

    # initialize collection vectors
    accept <- 0
    acceptS <- 0
    nrun <- B
    f2 <- rep(NA,nrun)
    delta <- rep(NA,nrun)
    mucollect_A <- matrix(0,length(alllocs),nrun)
    mucollect_B <- matrix(0,length(alllocs),nrun)
    phivec <- matrix(0,nrun,5)

    # Run chain
    for(i in 1:nrun){
      # Draw mu_A|.
      mu_A <- draw_mu(dat_A,nreps_A,Sigma_A)
      # Draw mu_B|.
      mu_B <- draw_mu(dat_B,nreps_B,Sigma_B)
      # Draw tau2S
      tau2S <- draw_tau2S()
      S <- tau2S*covf2(D,phiS)

      # Draw phiS
      phiS <- draw_phiS()
      S <- tau2S*covf2(D,phiS)+diag(nloc)*.0001

      # Draw tau2_A
      tau2_A <- draw_tau2(dat_A,mu_A,nreps_A)
      Sigma_A <- tau2_A*covdelta(D,phi)
      # Draw tau2_B
      tau2_B <- draw_tau2(dat_B,mu_B,nreps_B)
      Sigma_B <- tau2_B*covdelta(D,phi)

      # Draw phi
      phi <- draw_phi()
      phivec[i,] <- c(phi,phiS,sqrt(tau2_A),sqrt(tau2_B),sqrt(tau2S))

      # posterior predictive for f2
      mustar_A <- mu_post(mu_A,tau2_A)
      mustar_B <- mu_post(mu_B,tau2_B)
      mucollect_A[,i] <- mustar_A[order(alllocs)]
      mucollect_B[,i] <- mustar_B[order(alllocs)]
      f2[i] <- 50*log10(1/sqrt(1+(1/length(mustar_A))*sum((mustar_A-mustar_B)^2))*100)
      delta[i] <- max(abs(mustar_A-mustar_B))

      if(adaptive & i>25){
        prop_phi = sd(phivec[1:i,1])/3
        prop_phiS = sd(phivec[1:i,2])/3
      }
    }
    cov_pars <- data.frame(phivec)
    names(cov_pars) <- c("phi","psi","sigma_R","sigma_T","tau")
    mucollect_A <- data.frame(mucollect_A)
    names(mucollect_A) <- rep("muR", ncol(mucollect_A))
    mucollect_B <- data.frame(mucollect_B)
    names(mucollect_B) <- rep("muT", ncol(mucollect_B))
    mu_pars = cbind(mucollect_A,mucollect_B)


    return(list(delta = mean(delta), f2 = mean(f2), mcmc_chains = list(delta=delta,f2=f2,mu_pars = mu_pars,cov_pars = cov_pars)))
  }
}
