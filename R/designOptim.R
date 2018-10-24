#' Optimize N and R
#'
#' @param EDR.result a MethylSeqDesign object
#' @param pilot.R sequencing depth of pilot data in number of lanes
#' @param R The sequencing depth of interest in number of lanes, input is a vector.
#' @param N The sample size of interest per group, input is a vector.
#' equals the number of target sample size.
#' @param targetEDR the target EDR to achieve
#' @param budget the budget limit
#' @param price_per_lane the price per lane for Methyl-Seq
#' @param fix_price_per_sample the fixed part of the price per sample 
#' @depth_per_lane total number of reads per lane in millions
#' @return a list of three elements. The first element gives the best N and R combination, and optimal EDR/cost given budget/target EDR. 
#' The second element is a plot of admissible/inadmissible N and R combinations. The third element is a plot of cost and EDR.
#' @export designOptim
#'
#' @examples
designOptim <- function(EDR.result, pilot.R, R = c(1/10, 1/8, 1/6, 1/4, 1/2, 1, 1.5, 2), N = seq(4, 50, 2), targetEDR = NULL, budget = NULL, price_per_lane = 1500, 
    fix_price_per_sample = 300, depth_per_lane = 250) {
    if (!"ggplot2" %in% rownames(installed.packages())) {
        install.packages("ggplot2")
    }
    if (is.null(targetEDR) & is.null(budget)) {
        stop("Either targetEDR or buget should be given")
    }
    R = R[pilot.R >= R]
    sample.size <- t(matrix(rep(2 * N, length(R)), ncol = length(R)))
    cost.matrix <- ceiling(sweep(sample.size, 1, depth_per_lane * R, FUN = "*")/depth_per_lane) * price_per_lane + sample.size * fix_price_per_sample
    # Given Power, minimize budget
    if (is.null(budget) & !is.null(targetEDR)) {
        # Find the minimun budget
        min.budget <- min(cost.matrix[EDR.result$EDR >= targetEDR])
        EDR.achieve <- max(EDR.result$EDR[cost.matrix <= min.budget])
        optim.R <- R[which(cost.matrix == min.budget & EDR.result$EDR == EDR.achieve, arr.ind = TRUE)[1]]
        optim.N <- N[which(cost.matrix == min.budget & EDR.result$EDR == EDR.achieve, arr.ind = TRUE)[2]]
        optim.res = list(optim.N = optim.N, optim.R = optim.R, EDR.target = targetEDR, EDR.achive = EDR.achieve, min.budget = min.budget)
        # Find all admissible N and R
        admis.matrix <- cost.matrix
        for (i in 1:nrow(admis.matrix)) {
            for (j in 1:ncol(admis.matrix)) {
                admis.matrix[i, j] <- ifelse(sum(cost.matrix < cost.matrix[i, j] & EDR.result$EDR > EDR.result$EDR[i, j]) > 0, 0, 1)
            }
        }
        # plot1
        cols <- ifelse(as.numeric(admis.matrix) == 0, "grey", "black")
        if (0) {
            plot(x = rep(N, each = length(R)), y = rep(R, length(N)), pch = 4, col = cols, xlab = "N", ylab = "R", main = "case2")
            points(x = optim.N, y = optim.R, col = 2, pch = 4)
            
        }
        dat <- data.frame(N = rep(N, each = length(R)), R = rep(R, length(N)))
        library(ggplot2)
        p1 <- ggplot(dat, aes(x = N, y = R)) + geom_point(col = cols, size = 3, shape = 4) + geom_point(aes(x = optim.N, y = optim.R), colour = "red", 
            size = 3, shape = 4) + theme_bw()
        
        # plot2
        cols <- ifelse(as.numeric(admis.matrix) == 0, "grey", "black")
        if (0) {
            plot(x = as.numeric(cost.matrix), y = as.numeric(EDR.result$EDR), pch = 4, col = cols, xlab = "budget", ylab = "EDR", xlim = c(min.budget - 
                10000, min.budget + 10000), ylim = c(targetEDR - 0.15, min(targetEDR + 0.15, 1)), main = "case2")
            points(x = min.budget, y = EDR.achieve, pch = 4, col = 2)
            abline(h = targetEDR, lty = 2)
        }
        dat <- data.frame(cost = as.numeric(cost.matrix), EDR = as.numeric(EDR.result$EDR))
        library(ggplot2)
        p2 <- ggplot(dat, aes(x = cost, y = EDR)) + geom_point(col = cols, size = 3, shape = 4) + geom_point(aes(x = min.budget, y = EDR.achieve), 
            colour = "red", size = 3, shape = 4) + geom_hline(yintercept = targetEDR, colour = "blue", size = 1, linetype = "dashed") + annotate("text", 
            median(dat$cost), targetEDR, vjust = -1, label = paste0("target EDR = ", targetEDR)) + theme_bw()
    }
    # Given budget, maximize power
    if (!is.null(budget) & is.null(targetEDR)) {
        max.EDR <- max(EDR.result$EDR[cost.matrix <= budget])
        budget.achieve <- min(cost.matrix[EDR.result$EDR >= max.EDR])
        optim.R <- R[which(EDR.result$EDR == max.EDR & cost.matrix == budget.achieve, arr.ind = TRUE)[1]]
        optim.N <- N[which(EDR.result$EDR == max.EDR & cost.matrix == budget.achieve, arr.ind = TRUE)[2]]
        optim.res <- list(max.EDR = max.EDR, optim.N = optim.N, optim.R = optim.R)
        optim.res = list(optim.N = optim.N, optim.R = optim.R, budget = budget, cost = budget.achieve, max.EDR = max.EDR)
        # Find all admissible N and R
        admis.matrix <- cost.matrix
        for (i in 1:nrow(admis.matrix)) {
            for (j in 1:ncol(admis.matrix)) {
                admis.matrix[i, j] <- ifelse(sum(cost.matrix < cost.matrix[i, j] & EDR.result$EDR > EDR.result$EDR[i, j]) > 0, 0, 1)
            }
        }
        # plot1
        cols <- ifelse(as.numeric(admis.matrix) == 0, "grey", "black")
        if (0) {
            plot(x = rep(N, each = length(R)), y = rep(R, length(N)), pch = 4, col = cols, xlab = "N", ylab = "R", main = "case1")
            points(x = optim.N, y = optim.R, col = 2, pch = 4)
        }
        dat <- data.frame(N = rep(N, each = length(R)), R = rep(R, length(N)))
        library(ggplot2)
        p1 <- ggplot(dat, aes(x = N, y = R)) + geom_point(col = cols, size = 3, shape = 4) + geom_point(aes(x = optim.N, y = optim.R), colour = "red", 
            size = 3, shape = 4) + theme_bw()
        
        
        # plot2
        cols <- ifelse(as.numeric(admis.matrix) == 0, "grey", "black")
        if (0) {
            plot(x = as.numeric(cost.matrix), y = as.numeric(EDR.result$EDR), pch = 4, col = cols, ylim = c(max.EDR - 0.15, min((max.EDR + 0.15), 
                1)), xlim = c(budget - 10000, budget + 10000), xlab = "budget", ylab = "EDR", main = "case1")
            points(x = budget.achieve, y = max.EDR, pch = 4, col = 2)
            abline(v = budget, lty = 2)
        }
        dat <- data.frame(cost = as.numeric(cost.matrix), EDR = as.numeric(EDR.result$EDR))
        library(ggplot2)
        p2 <- ggplot(dat, aes(x = cost, y = EDR)) + geom_point(col = cols, size = 3, shape = 4) + geom_point(aes(x = budget.achieve, y = max.EDR), 
            colour = "red", size = 3, shape = 4) + geom_vline(xintercept = budget, colour = "blue", size = 1, linetype = "dashed") + annotate("text", 
            budget, median(dat$EDR), vjust = -1, label = paste0("budget = ", budget)) + theme_bw()
    }
    return(list(res = optim.res, plot1 = p1, plot2 = p2))
}
