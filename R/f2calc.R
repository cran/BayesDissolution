#' Calculation of the f2 Statistic
#'
#' This function calculates the f2 statistic as described in Moore and Flanner (1996).
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#'  the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#'  For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#'  measurements sorted in time.
#' @return The function returns the f2 statistic calculated from the observed dissolution data.
#' @note Use the plotdiss() function to visually check if it's appropriate to calculate the f2 statistic.
#' @references Moore, J.W. and Flanner, H.H. (1996). Mathematical comparison of distribution profiles. Pharmaceutical Technology, 20(6):64-74.
#' @examples
#' ### dis_data comes loaded with the package
#' f2calc(dis_data)
#'
#' @export
f2calc <- function(dis_data){

  if(!is.data.frame(dis_data)){
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

    Ybar_R <- apply(dat_R, 2, mean)
    Ybar_T <- apply(dat_T, 2, mean)

    f2 <- 50 * log10(1 / sqrt(1 + (1 / nlocs) * sum((Ybar_R - Ybar_T)^2)) * 100)

    return(f2)
  }
}
