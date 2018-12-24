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
    N <- as.numeric(sapply(strsplit(colnames(x$EDR), " vs "), function(x) x[1]))
    prop <- as.numeric(rownames(x$EDR))
    plot3d.data <- list(sample.size = N, depth.per.sample = prop, estimated.EDR = x$EDR)
    
    p <- plot_ly(x = plot3d.data$sample.size, y = plot3d.data$depth.per.sample, z = plot3d.data$estimated.EDR) %>% add_surface() %>% 
        layout(scene = list(xaxis = list(title = "N"), yaxis = list(title = "R"), zaxis = list(title = "EDR")))
    return(p)
}
