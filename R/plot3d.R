#' 3D plot of EDR
#'
#' @param x a MethylSeqDesign object
#' @return a 3D plot of EDR, x axls is the sample size, y axls is the propotion of subsampling, z axls is the EDR.
#' @export plot3d
#'
#' @examples
plot3d <- function(x) {
    if (!"plotly" %in% rownames(installed.packages())) {
        install.packages("plotly")
    }
    library(plotly)
    estimated.EDR <- as.matrix(x$EDR)
    p <- plot_ly(z = ~estimated.EDR) %>% add_surface() %>% layout(scene = list(xaxis = list(title = "Sample size"), yaxis = list(title = "Seqencing depth per sample"), 
        zaxis = list(title = "estimated EDR")))
    return(p)
}
