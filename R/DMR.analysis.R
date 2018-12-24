#' Differential methylation analysis on pilot data
#'
#' @param N0  the sample size of the pilot data per group. If the two groups of the pilot data have different sample size, input can be a vector of length 2 
#' specifying the sample size for each group. 
#' @param cov.matrix a numeric matrix or data frame of coverage, the rows are genes and columns are samples.
#' @param methyl.matrix a numeric matrix or data frame  of methylated read count, the rows are genes and columns are samples. A one to one correspondence between the matrix of methylated read count and the matrix of coverage is required.
#' @param R the targeted sequencing depth in number of lanes. Input can be a vector. R is no larger than pilot.R.
#' @param piot.R Sequencing depth of pilot data in number of lanes.
#' @return a list of three element, the first element is the p values from DM analysis, the second element is a matrix of 
#' GLM model parameters and estimated tagwise dispersion for each region. The third element is a matrix of ratios for different sequencing depth. 
#' @export DMR.analysis
#'
#' @examples
DMR.analysis <- function(N0, cov.matrix, methyl.matrix, R, pilot.R) {
    if (class(N0) != "numeric" | length(N0) > 2) {
        stop("Argument N0 is not correctly specified")
    } else if (length(N0) == 1) {
        N0 = rep(N0, 2)
    }
    if (class(cov.matrix) == "data.frame") {
        cov.matrix = as.matrix(cov.matrix)
    }
    if (class(methyl.matrix) == "data.frame") {
        methyl.matrix = as.matrix(methyl.matrix)
    }
    if (class(cov.matrix) != "matrix") {
        stop("input coverage matrix should be either matrix or data frame")
    }
    if (class(methyl.matrix) != c("matrix")) {
        stop("input coverage matrix should be either matrix or data frame")
    }
    if (!"DropletUtils" %in% rownames(installed.packages())) {
        install.packages("DropletUtils")
    }
    library(DropletUtils)
    
    print(paste(N0, sep = "", collapse = " vs "))
    fai.est <- estimate.fai(x = N0, cov.matrix, methyl.matrix)
    test.result <- Z.wald(x = N0, cov.matrix, methyl.matrix, fai = fai.est, a = 1)
    p.values <- test.result[[1]]
    # calculate ratio of psai
    factorR <- matrix(, nrow = length(p.values), ncol = length(R))
    colnames(factorR) <- R
    for (i in 1:length(R)) {
        cov.matrix1 <- downsampleMatrix(cov.matrix, prop = R[i]/pilot.R)
        data0 = data1 = cov.matrix
        for (j in 1:sum(N0)) {
            m0 <- cov.matrix[, j]
            m1 <- cov.matrix1[, j]
            data0[, j] <- m0/(1 + (m0 - 1) * fai.est)
            data1[, j] <- m1/(1 + (m1 - 1) * fai.est)
        }
        
        mean.group1 <- apply(data0[, 1:N0[1], drop = F], 1, mean)
        mean.group2 <- apply(data0[, (N0[1] + 1):sum(N0), drop = F], 1, mean)
        psai0 <- (mean.group1 + mean.group2)/(mean.group1 * mean.group2)
        
        mean.group1 <- apply(data1[, 1:N0[1], drop = F], 1, mean)
        mean.group2 <- apply(data1[, (N0[1] + 1):sum(N0), drop = F], 1, mean)
        psai1 <- (mean.group1 + mean.group2)/(mean.group1 * mean.group2)
        
        factorR[, i] = psai0/psai1
    }
    
    model <- cbind(beta0 = test.result[[2]], beta1 = test.result[[3]], fai = fai.est)
    return(list(p.values = p.values, model = model, factorR = factorR))
}



# Z-value + wald-test
Z.wald <- function(x, cov.matrix, sim.meth, fai, a = 1) {
    D <- sum(x)
    G <- dim(cov.matrix)[1]
    if (dim(cov.matrix)[2] != D) {
        stop("sample size and coverage data matrix columns differ")
    }
    if (dim(sim.meth)[2] != D) {
        stop("sample size and methylation data matrix columns differ")
    }
    n.groups <- length(x)
    if (n.groups > 2) {
        stop("Currently more than 2 groups is not allowed")
    }
    
    # Assign each group size
    for (i in 1:n.groups) {
        assign(paste0("D", i), x[i])
    }
    
    ### calculate beta1-hat
    chi2 <- rep(0, G)
    beta0.hat <- rep(0, G)
    beta1.hat <- rep(0, G)
    
    for (i in 1:G) {
        tryCatch({
            X <- matrix(c(rep(1, D), rep(1, D1), rep(0, D2)), ncol = 2)
            
            if (n.groups > 2) {
                c <- D1
                for (j in 2:(n.groups - 1)) {
                  col <- c(rep(0, c), rep(1, get(paste0("D", j))), rep(0, D - c - get(paste0("D", j))))
                  X <- cbind(X, col)
                  c = c + get(paste0("D", j))
                }
            }
            
            V.inverse <- diag(cov.matrix[i, ]/(1 + (cov.matrix[i, ] - 1) * fai))
            Z <- asin(2 * (sim.meth[i, ] + a)/(cov.matrix[i, ] + 2 * a) - 1)
            sigma <- solve(t(X) %*% V.inverse %*% X)
            beta.hat <- sigma %*% t(X) %*% V.inverse %*% Z
            beta0.hat[i] <- beta.hat[1]
            beta1.hat[i] <- beta.hat[2]
            
            C <- diag(rep(1, n.groups))[, -1]
            var <- t(C) %*% sigma %*% C
            chi2[i] <- t(t(C) %*% beta.hat) %*% solve(var) %*% (t(C) %*% beta.hat)
        }, error = function(err) {
            # error handler picks up where error was generated
            print(paste("ERROR:  ", err))
        })
    }
    # Calculate test statistic, p and q values
    p.Zw <- pchisq(chi2, df = n.groups - 1, lower.tail = F)
    return(list(p.Zw = p.Zw, beta0 = beta0.hat, beta1 = beta1.hat, z.statistic = chi2))
}



## Estimate fai using park's method
estimate.fai <- function(cov.matrix, methyl.matrix, a = 1, x) {
    G <- dim(cov.matrix)[1]
    D <- dim(cov.matrix)[2]
    n.groups <- length(x)
    # Assign each group size
    for (i in 1:n.groups) {
        assign(paste0("D", i), x[i])
    }
    
    if (D == 2) 
        fai.est = 0.001 else {
        fai <- rep(0, G)
        for (i in 1:G) {
            tryCatch({
                V0.inverse <- diag(cov.matrix[i, ])
                # X matrix
                X <- matrix(c(rep(1, D), rep(1, D1), rep(0, D - D1)), ncol = 2)
                if (n.groups > 2) {
                  c <- D1
                  for (j in 2:(n.groups - 1)) {
                    col <- c(rep(0, c), rep(1, get(paste0("D", j))), rep(0, D - c - get(paste0("D", j))))
                    X <- cbind(X, col)
                    c = c + get(paste0("D", j))
                  }
                }
                Z <- asin(2 * (methyl.matrix[i, ] + a)/(cov.matrix[i, ] + 2 * a) - 1)
                sigma <- solve(t(X) %*% V0.inverse %*% X)
                beta.hat0 <- sigma %*% t(X) %*% V0.inverse %*% as.matrix(Z)
                chi2 <- sum(cov.matrix[i, ] * (Z - X %*% beta.hat0)^2)
                sigma2.hat <- chi2/(D - n.groups)
                fai[i] <- (D * (sigma2.hat - 1))/sum(cov.matrix[i, ] - 1)
            }, error = function(err) {
                # error handler picks up where error was generated
                print(paste("MY_ERROR:  ", err))
            })
        }
        fai <- fai[complete.cases(fai)]
        fai <- ifelse(fai < 0, 10^(-6), fai)
        fai <- ifelse(fai > 1, 1 - 10^(-6), fai)
        fai.est <- mean(fai)
    }
    return(fai.est)
}

