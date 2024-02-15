#' Dissolution Data Plot
#'
#' Minimalist ggplot function for plotting dissolution data sets.
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#' the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#' For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#' measurements sorted in time. Alternatively, the user may provide a data object of class dis_data containing the
#' dissolution data. See the \code{make_dis_data()} function for the particular structure of the data object.
#' @param show.mean logical; if \code{TRUE}, plot the connected mean dissolution values for each group.
#' @param show.SD logical; if \code{TRUE}, calculate the variance of the dissolution data at each time point for each group. The values are placed at the top of the plot over the corresponding time point.
#' @return The function returns a plot of the dissolution data.
#' @examples
#' ### dis_data comes loaded with the package
#' ggdissplot(dis_data)
#'
#' @importFrom graphics legend lines matplot mtext par
#' @export
ggdissplot <- function(dis_data, show.mean = FALSE, show.SD = FALSE){

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

    Group <- Time <- txt.ref <- txt.test <- y <- y.ref <- y.test <- NULL ## Hack to silence notes
    d <- data.frame(y = c(as.vector(unlist(dis_data$yRef)), as.vector(unlist(dis_data$yTest))),
                  Group = factor(rep(c("Reference", "Test"), each = dis_data$nTab * dis_data$nTime), levels = c("Reference", "Test")),
                  Time = factor(rep(1:dis_data$nTime, each = dis_data$nTab)),
                  Tablet = factor(rep(1:dis_data$nTab, dis_data$nTime))
    )

    p <- ggplot2::ggplot(d, ggplot2::aes(Time, y, col = Group, group = Group)) +
        ggplot2::geom_jitter(height = 0, width = 0.05) +
        ggplot2::xlab("Time point") +
        ggplot2::ylab("Dissolution (%)") +
        ggplot2::theme(legend.position = "top")

    if ( isTRUE(show.mean) ){
        d.mean <- data.frame(Time = factor(rep(1:dis_data$nTime, 2)),
                            Group = factor(rep(c("Reference", "Test"), each = dis_data$nTime)),
                            y = c(dis_data$ybarR, dis_data$ybarT)
                            )
        p <- p + ggplot2::geom_line(data = d.mean, linetype="dashed")
    }
    if ( isTRUE(show.SD) ){

        SD <- tapply(d$y, list(d$Time, d$Group), sd)
        d.sd <- data.frame(Time = c(0.7, 1:dis_data$nTime),
                            txt.ref =c("Ref", round(SD[,1], 1)),
                            txt.test = c("Test", round(SD[,2], 1)),
                            y.ref = max(d$y) + 0.1 * diff(range(d$y)),
                            y.test = max(d$y) + 0.05 * diff(range(d$y))
                            )


        p <- p + ggplot2::geom_text(data = d.sd, ggplot2::aes(Time, y.ref, label = txt.ref), inherit.aes=FALSE) +
          ggplot2::geom_text(data = d.sd, ggplot2::aes(Time, y.test, label = txt.test), inherit.aes = FALSE)
    }


  }
  return(p)

}
