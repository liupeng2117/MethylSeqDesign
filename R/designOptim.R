#' Optimize N and R
#'
#' @param EDR.result a MethylSeqDesign object
#' @param pilot.R sequencing depth of pilot data in millions
#' @param prop a numeric vector of propotion of subsampling
#' @param target.N a numeric vector or a matrix of target sample size. If it is a matrix, the number of columns is 2, and the number of rows
#' equals the number of target sample size.
#' @param targetEDR the target EDR want to achieve
#' @param budget the budget limit
#' @param price_per_lane the price per lane for Methyl-Seq
#' @param fix_price_per_sample the fixed part of the price per sample 
#' @depth_per_lane total number of reads per lane in millions
#' @return a list of the optimal EDR/budget, and the corresponding sample size in two groups and sequencing depth per sample in million reads.
#' @export designOptim
#'
#' @examples
designOptim <- function(EDR.result, pilot.R, prop, target.N, targetEDR = NULL, budget = NULL, price_per_lane = 2000, fix_price_per_sample = 100, 
    depth_per_lane = 150) {
    if (is.null(targetEDR) & is.null(budget)) {
        stop("Either targetEDR or buget should be given")
    }
    if (!class(target.N) %in% c("matrix", "numeric")) {
        stop("Argument target.N is not correctly specified")
    } else if (class(target.N) == "numeric") {
        target.N = matrix(rep(target.N, 2), ncol = 2)
    }
    
    sample.size <- t(matrix(rep(apply(target.N, 1, sum), length(prop)), ncol = length(prop)))
    cost.matrix <- ceiling(sweep(sample.size, 1, pilot.R * prop, FUN = "*")/depth_per_lane) * price_per_lane + sample.size * fix_price_per_sample
    # Given Power, optimize budget
    if (is.null(budget) & !is.null(targetEDR)) {
        cost.matrix[EDR.result$EDR <= targetEDR] <- Inf
        min.budget <- min(cost.matrix)
        optim.R <- prop[which(cost.matrix == min.budget, arr.ind = TRUE)[1]] * pilot.R
        optim.N <- target.N[which(cost.matrix == min.budget, arr.ind = TRUE)[2], ]
        optim.res <- list(min.budget = min.budget, optim.N = optim.N, optim.R = optim.R)
    }
    # Given budget, optimize power
    if (!is.null(budget) & is.null(targetEDR)) {
        EDR.result$EDR[cost.matrix >= budget] <- 0
        max.EDR <- max(EDR.result$EDR)
        optim.R <- prop[which(EDR.result$EDR == max.EDR, arr.ind = TRUE)[1]] * pilot.R
        optim.N <- N[which(EDR.result$EDR == max.EDR, arr.ind = TRUE)[2]]
        optim.res <- list(max.EDR = max.EDR, optim.N = optim.N, optim.R = optim.R)
    }
    return(optim.res)
}
