#' Estimate EDR based on p values
#'
#' @param res the result from DMR analysis
#' @param N0 the sample size of the pilot data per group. If the two groups of the pilot data have different sample size, input can be a vector of length 2
#' specifying the sample size for each group.
#' @param target.N the target sample size per group. Input can be a vector. If the targeted sample size is unbalanced, the input is a matrix of two columns
#' @param thresh.p the threshold for CBUM mixture model fitting, default is set to 0.005
#' @param FDR the level of False discovery rate to control, by default it is set to 0.05
#' @param M the time of parametric bootstrap, default is 10
#' @param restrict restrict lambda value, default is FALSE
#'
#' @return A list of 3 element. The first element is a matrix of EDR, the second element is a matrix of declared
#' Number of DMR, the third element is a matrix of empirical FDR.
#' @export Estimate.EDR.from.pilot
#'
#' @examples
Estimate.EDR.from.pilot <- function(res, N0, target.N, thresh.p = 0.005, FDR = 0.05, M = 10, restrict = FALSE) {
    p.values <- res$p.values
    # model<-res$model
    factorR <- res$factorR

    if (all(is.na(p.values))) {
        stop("all p-values were NA; nothing to compute")
    }
    if (any(is.na(p.values))) {
        p.values <- p.values[!is.na(p.values)]
    }
    if (min(p.values) == 0) {
        min.nonzero <- min(p.values[p.values > 0])
        p.values[p.values == 0] <- min.nonzero/2
    }
    if (max(p.values) == 1) {
        max.nonone <- max(p.values[p.values < 1])
        p.values[p.values == 1] <- 1 - (1 - max.nonone)/2
    }
    if (class(N0) != "numeric" | length(N0) > 2) {
        stop("Argument N0 is not correctly specified")
    } else if (length(N0) == 1) {
        N0 = rep(N0, 2)
    }
    if (!class(target.N) %in% c("matrix", "numeric")) {
        stop("Argument target.N is not correctly specified")
    } else if (class(target.N) == "numeric") {
        target.N = matrix(rep(target.N, 2), ncol = 2)
    }

    ngenes = length(p.values)
    # step 2 fit mixture model using CBUM+CDD method
    Fitted.Model = MixtureModel.Fittting.pi0(p.values, thresh.p = thresh.p, restrict = restrict)

    # step3 caluclate posterior probability of being DMR
    lambda = as.numeric(Fitted.Model[1])
    r = as.numeric(Fitted.Model[2])
    s = as.numeric(Fitted.Model[3])
    posterior = lambda/(dbeta(p.values, r, s) * (1 - lambda) + lambda)  #prob of being non-DE

    if (lambda > 0.99 | lambda < 0.5) {
        stop(paste("lambda >0.99 or <0.5,lambda=", Fitted.Model[1]))
    }

    parameter = list(n.old = N0, n.new = target.N)

    # step 4: calculate estimated EDR using parametric bootstrap method
    EDR = DeclareDMR = FDR.matrix = matrix(, nrow = ncol(factorR), ncol = nrow(target.N))
    rownames(EDR) = rownames(DeclareDMR) = rownames(FDR.matrix) = round(as.numeric(colnames(factorR)),
        3)
    colnames(EDR) = colnames(DeclareDMR) = colnames(FDR.matrix) = apply(target.N, 1, function(x) paste(x,
        sep = "", collapse = " vs "))
    for (i in 1:ncol(factorR)) {
        Result = lapply(1:M, function(x) Resampling(target.N, posterior, p.values, parameter, ngenes,
            FDR, factorR[, i]))
        ave.result <- Reduce("+", Result)/length(Result)
        EDR[i, ] <- ave.result["EDR", ]
        DeclareDMR[i, ] <- round(ave.result["DeclareDMR", ])
        FDR.matrix[i, ] <- ave.result["FDR", ]
    }
    pred <- list(EDR = EDR, DeclareDMR = DeclareDMR, FDR.matrix = FDR.matrix)
    class(pred) <- append(class(pred), "MethylSeqDesign")
    return(pred)
}



# Use censored beta-uniform mixture model to estimate pi0(lambda)
MixtureModel.Fittting.pi0 <- function(p.values, restrict = FALSE, l.upper = 0.99, thresh.p = 0.005) {
    if (!"limma" %in% rownames(installed.packages())) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("limma")
    }
    if (!"pi0" %in% rownames(installed.packages())) {
        devtools::install_github("gitlongor/pi0")
    }

    library(limma)
    library(pi0)

    fBUM <- function(z, d, lambda) {
        p.values = d
        r = z[1]
        s = z[2]
        Target.Func = -sum(log(lambda + (1 - lambda) * p.values^(r - 1) * (1 - p.values)^(s - 1)/beta(r,
            s)))
        return(Target.Func)
    }

    ## inital value(MME and non-parametric)
    mean.p = mean(p.values, na.rm = T)
    var.p = var(p.values, na.rm = T)
    init.r = ((1 - mean.p) * mean.p^2 - mean.p * var.p)/var.p
    init.s = ((1 - mean.p)^2 * mean.p - (1 - mean.p) * var.p)/var.p
    init = c(max(0, min(init.r, 0.9)), max(1, init.s))

    lambda0 = CBUM(p.values, thresh.censor = thresh.p, niter = 1000)

    if (attributes(lambda0)$converged == FALSE) {
        print("Fails to estimate lambda, CBUM model does not converge. Use CDD instead")
        lambda0 = convest(p.values)[[1]]
    }

    if (restrict) {
        lambda = max(min(lambda0, l.upper), 0.8)
        r.upper = 0.9
        s.upper = Inf
        r.lower = 0
        s.lower = 1
    } else {
        lambda = min(lambda0, 0.99)
        r.upper = 1
        s.upper = Inf
        r.lower = 0
        s.lower = 1
    }
    tryCatch({
        pars = optim(init, fn = fBUM, d = p.values, lambda = lambda, method = "L-BFGS-B", upper = c(r.upper,
            s.upper), lower = c(r.lower, s.lower))$par
        LL = sum(log(lambda + (1 - lambda) * p.values^(pars[1] - 1) * (1 - p.values)^(pars[2] - 1)/beta(pars[1],
            pars[2])))
        out.par <- c(lambda, pars, LL)
        names(out.par) <- c("lambda", "r", "s", "LL")
        return(out.par)
    }, error = function(e) {
        print("Fail to optimize r and s, return initial values.")
        LL = sum(log(lambda + (1 - lambda) * p.values^(init[1] - 1) * (1 - p.values)^(init[2] - 1)/beta(init[1],
            init[2])))
        out.par <- c(lambda, init, LL)
        names(out.par) <- c("lambda", "r", "s", "LL")
        return(out.par)
    })
}



# Get DMR when controling true FDR at 0.05 level
get.dm.regions <- function(p, level = 0.05, dmr) {
    p.de <- p[dmr]
    p.ordered <- p[order(p)]
    tdr <- rep(0, length(p))
    for (i in 1:length(p)) {
        tdr[i] = sum(p.de <= p.ordered[i])/sum(p.ordered <= p.ordered[i])
    }
    id <- which(tdr >= 1 - level)
    if (length(id) == 0) {
        return(print("0 DMR found"))
    } else num.de <- max(id)
    claim.dmr <- order(p)[1:num.de]
    return(claim.dmr)
}

# calculate estimated EDR using parametric bootstrap method
Resampling <- function(target.N, posterior, p.values, parameter, ngenes, FDR, factorR) {

    DE_status_posterior = sapply(1:length(posterior), function(x) sample(c(FALSE, TRUE), 1, prob = c(posterior[x],
        1 - posterior[x]), replace = TRUE))

    transform.p.values.norm <- function(p.values.each, DE_status_each, parameter, factorR.each) {
        n.old <- parameter[[1]]
        n.new <- parameter[[2]]

        if (DE_status_each) {
            statistic.old = qnorm(p.values.each/2, lower.tail = FALSE)
            statistic.new = statistic.old * sqrt((apply(n.new, 1, prod)/prod(n.old)) * (sum(n.old)/apply(n.new,
                1, sum))) * sqrt(factorR.each)
            p.values.star.star = (1 - pnorm(abs(statistic.new))) * 2
            return(p.values.star.star)
        } else return(rep(p.values.each, nrow(n.new)))
    }

    p.values.star.posterior = sapply(1:ngenes, function(x) transform.p.values.norm(p.values[x], DE_status_posterior[x],
        parameter = parameter, factorR.each = factorR[x]))

    p.values.star.posterior = matrix(p.values.star.posterior, ncol = nrow(target.N), byrow = T)

    Estimate_Posterior <- function(p.values.star.star) {
        if (min(p.values.star.star, na.rm = T) == 0) {
            min.nonzero <- min(p.values.star.star[p.values.star.star > 0])
            p.values.star.star[p.values.star.star == 0] <- min.nonzero/2
        }
        if (max(p.values.star.star, na.rm = T) == 1) {
            max.non1 <- max(p.values.star.star[p.values.star.star < 1])
            p.values.star.star[p.values.star.star == 1] <- (max.non1 + 1)/2
        }
        ## empirical FDR control
        p.DE = p.values.star.star[DE_status_posterior]
        p.nonDE = p.values.star.star[!DE_status_posterior]
        p.unique = sort(unique(p.values.star.star))

        FDR.unique <- vector(length = length(p.unique))
        for (i in 1:length(p.unique)) {
            FDR.unique[i] = sum(p.nonDE <= p.unique[i])/sum(p.values.star.star <= p.unique[i])
        }
        index = which(FDR.unique <= FDR)
        p.values.cut = if (length(index)) {
            p.unique[max(index)]
        } else 0  #all p value cutoffs exceed the specified FDR level, fails to control FDR
        FDR.cut <- if (length(index)) {
            FDR.unique[max(index)]
        } else 0

        if (min(p.values.star.star) > p.values.cut | is.na(p.values.cut)) {
            Declare_status = rep("nonDMR", ngenes)
        } else {
            Declare_status = DE_status_posterior
            Declare_status[which(p.values.star.star <= p.values.cut)] = "DMR"
            Declare_status[-which(p.values.star.star <= p.values.cut)] = "nonDMR"
        }
        A = sum((Declare_status == "nonDMR") * (!DE_status_posterior))
        B = sum((Declare_status == "nonDMR") * (DE_status_posterior))
        C = sum((Declare_status == "DMR") * (!DE_status_posterior))
        D = sum((Declare_status == "DMR") * (DE_status_posterior))
        if ((C + D) == 0) {
            TP_hat_post = 0
            ## no declared DMR
        } else {
            TP_hat_post = D/(C + D)
        }
        if ((A + B) == 0) {
            TN_hat_post = 0
            ## no declared DMR
        } else {
            TN_hat_post = A/(A + B)
        }
        EDR_post = D/(B + D)
        Declare_post = sum(Declare_status == "DMR")
        return(c(TP_hat_post, TN_hat_post, EDR_post, Declare_post, A = A, B = B, C = C, D = D, p.values.cut = p.values.cut,
            FDR = FDR.cut))
    }


    Estimate.Posterior.Result = matrix(apply(p.values.star.posterior, 2, Estimate_Posterior), ncol = nrow(target.N))

    Estimate.Posterior.Result <- round(as.matrix(Estimate.Posterior.Result), 3)

    row.names(Estimate.Posterior.Result) = c("TP", "TN", "EDR", "DeclareDMR", "A", "B", "C", "D", "FDRcut",
        "FDR")
    colnames(Estimate.Posterior.Result) = apply(target.N, 1, function(x) paste(x, sep = "", collapse = " vs "))

    return(Estimate.Posterior.Result)
}
