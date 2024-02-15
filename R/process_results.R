#' Helper function for processing results
#'
#' This function helps process the final results for the different f2 functions (e.g., f2bayes).
#'
#' @param name A character string denoting the type of method used to calculate the interval.
#' @param f2.dist A vector of samples for the F2 parameter or f2 statistic.
#' @param ci.type The type of confidence, or credible, interval to return. The option \code{quantile} returns sample quantile based intervals, while the option \code{HPD} returns intervals based on highest density regions.
#' @param level The confidence level or probability associated with the confidence or credible interval, respectively. Must be a value between 0 and 1.
#' @param get.dist logical; if \code{TRUE}, returns the samples of the distribution.
#' @return The function returns a data object of class dis_data.
#' @examples
#' ### dis_data comes loaded with the package
#' out1 <- f2bayes(dis_data, prob = 0.9, B = 1000, get.dist = TRUE)
#'
#' out2 <- process_results("bayes", out1$f2.dist, level = 0.9)
#'
#' ### out1 and out2 should contain the results for the info and intervals
#'
#' @export
process_results = function(name, f2.dist, ci.type = c("quantile", "HPD"), level, get.dist = FALSE){

  out <- list()
  if (length(level) == 1){
    level <- c(0.5 * (1 - level), 1 - 0.5 * (1 - level))
  }

  out$info <- data.frame(type = name, K = length(f2.dist), level = level)

  ci.type <- match.arg(ci.type, several.ok = TRUE)

  if(any(ci.type == "quantile")){
    out$ci.quantile <- stats::quantile(f2.dist, p = level)
  }
  if(any(ci.type == "HPD")){
    out$ci.HPD <- coda::HPDinterval(coda::as.mcmc(f2.dist), prob = diff(level))
  }
  if(get.dist){
    out$f2.dist <- f2.dist
  }
  return(out)

}
