#' Contour plot of EDR
#'
#' @param x a MethylSeqDesign object
#' @return a contour plot of EDR, x axls is the sample size and y axls is the propotion of subsampling.
#' @export plotContour
#' 
#' @examples
plotContour <- function(x) {
    if (!"plotly" %in% rownames(installed.packages())) {
        install.packages("plotly")
    }
    library(plotly)
    N <- as.numeric(sapply(strsplit(colnames(x$EDR), " vs "), function(x) x[1]))
    prop <- as.numeric(rownames(x$EDR))
    plot2d.data <- data.frame(sample.size = as.numeric(sapply(N, function(x) rep(x, length(prop)))), depth.per.sample = rep(prop, length(N)), 
        estimated.EDR = as.numeric(x$EDR))
    p <- plot_ly(plot2d.data, x = ~sample.size, y = ~depth.per.sample, z = ~estimated.EDR, type = "contour")
    return(p)
}


