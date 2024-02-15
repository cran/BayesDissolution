#' Dissolution Data Plot
#'
#' This function plots dissolution data sets.
#'
#' @param dis_data A data frame containing the dissolution data. The first column of the data frame should denote
#'  the group labels identifying whether a given dissolution belongs to the "reference" or "test" formulation group.
#'  For a given dissolution run, the remaining columns of the data frame contains the individual run's dissolution
#'  measurements sorted in time. Alternatively, the user may provide a data object of class dis_data containing the
#'  dissolution data. See the \code{make_dis_data()} function for the particular structure of the data object.
#' @param tp An optional vector of time points at which the dissolution data is measured at.
#' @param pch A vector of two elements specifying the plotting character to use for each group. If only one value is passed then the plotting character is the same for both groups.
#' @param color  A vector of two elements specifying the color in the plot to associate with each group. If only one value is passed then the color choice is the same for both groups.
#' @param groups A vector of two elements specifying the name to use for each group in the plot.
#' @param legend_location A string that denotes the location of where the legend should appear. Possible options are "left", "top", "bottom", "right", and any logical combination of the four, e.g., "bottomright" or "topleft".
#' @param xlab A string specifying the x-axis label.
#' @param ylab A string specifying the y-axis label.
#' @param mean logical; if \code{TRUE}, plot the connected mean dissolution values for each group
#' @param var logical; if \code{TRUE}, calculate the variance of the dissolution data at each time point for each group. The values are placed at the top of the plot over the corresponding time point.
#' @param var_label logical; if \code{TRUE}, use the group labels when printing out the variances.
#' @param ... other graphical parameters commonly found in \link[graphics]{plot.default}
#' @return The function returns a plot of the dissolution data.
#' @examples
#' ### dis_data comes loaded with the package
#' dissplot(dis_data)
#'
#' @importFrom graphics legend lines matplot mtext par
#' @export
dissplot <- function(dis_data, tp = NULL, pch = c(19, 17), color = c("gray65", "black"), groups = c("Reference", "Test"), legend_location = "bottomright",
                     xlab = "Time Points", ylab = "Percentage Dissolved", mean = FALSE, var = FALSE, var_label = TRUE, ...){

  if(class(dis_data)[1] == "dis_data"){
    dis_data <- data.frame(rbind(data.frame(group = "Reference", dis_data$yRef), data.frame(group = "Test", dis_data$yTest)))
  }
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
    nreps <- nrow(dat_R)

    Ybar_R <- apply(dat_R, 2, mean)
    Ybar_T <- apply(dat_T, 2, mean)

    s2_R <- round(apply(dat_R, 2, var), 2)
    s2_T <- round(apply(dat_T, 2, var), 2)

    if(is.null(tp)){
      tp <- 1:nlocs
    }

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(ps=15)
    matplot(tp, t(X[,-1]), pch = rep(pch, each = nreps),
            col = rep(color, each = nreps),
            xlab = xlab, ylab = ylab, ...)
    if(isTRUE(mean)){
      lines(tp, Ybar_R, col = color[1], lty = 1)
      lines(tp, Ybar_T, col = color[2], lty = 2)
    }
    if(isTRUE(var)){
      if(isTRUE(var_label)){
        mtext(c(groups[1], as.character(s2_R)), at = c(0.5, tp), side = 3, line = 1, col = color[1], cex = .75)
        mtext(c(groups[2], as.character(s2_T)), at = c(0.5, tp), side = 3, line = 2, col = color[2], cex = .75)
      }else{
        mtext(as.character(s2_R), at = tp, side = 3, line = 1, col = color[1], cex = .75)
        mtext(as.character(s2_T), at = tp, side = 3, line = 2, col = color[2], cex = .75)
      }
    }
    legend(legend_location, groups, pch = pch,col = color, bty = "n")
  }
}
