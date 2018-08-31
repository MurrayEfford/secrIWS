## Inputs 

## Number of coefficients nterms
## Number of occasions    t
## Number of individuals  Mt1
## Number of covariates   ncov
## Covariate values       zi[i,k] 
## Model types            Mt, Mb
## Time covariate         tj[j]
## Previously marked?     zij, Mij
## Capture data           cij


huggins.fit <- function (CH, 
                         Mb = c('None', 'LearnedResponse', 'MarkovResponse'),
                         Mt = c('None', 'TimeSpecific', 'TimeCovariate'),
                         zi = NULL,
                         tj = NULL,
                         ...) {
    
    negLLHuggins <- function (beta) {
        
        ## beta[1] is always beta0 (constant)
        ## beta[nterms] is always alpha (coefficient of zij OR Mij)
        
        result <- 0
        alpha <- beta[nterms]
        pij <- matrix(beta[1], nrow = Mt1, ncol = t)
        if (ncov>0) {
            zcov <- zi %*% beta[2:(ncov+1)]
            pij <- sweep (pij, FUN = '+', MARGIN = 1, STATS = zcov)
        }
        if (Mt == "TimeSpecific") {
            tspecific <- c(beta[(ncov+2):(ncov+t)], 0)
            pij <- sweep (pij, FUN = '+', MARGIN = 2, STATS = tspecific)        
        }
        if (Mt == "TimeCovariate") {
            tcov <- beta[ncov+2] * tj
            pij <- sweep (pij, FUN = '+', MARGIN = 2, STATS = tcov)        
        }
        
        ## naive p
        pijstar <- pij   
        ## learned response p
        pij <- switch (Mb,
                       LearnedResponse = pij + alpha * zij,
                       MarkovResponse  = pij + alpha * Mij,
                       pij)
        ## untransform
        pij <- exp(pij) / (1 + exp(pij))
        pijstar <- exp(pijstar) / (1 + exp(pijstar))
        
        ## should this have been Mij? bug in Density?
        if (Mb=="Removal") pij[zij>0.5] <- 0
        
        for (i in 1:Mt1) {
            for (j in 1:t) {
                prodx <- prod(1-pijstar[i,j:t])            
                gamij <- pij[i,j] / (1 - (1-zij[i,j]) * prodx)
                if (cij[i,j]> 0.5)
                    result <- result - log ( gamij)         ## caught
                else 
                    result <- result - log ( 1 - gamij)   ## not caught       
            }
        }
        result
    }
    Mb <- match.arg(Mb)
    Mt <- match.arg(Mt)
    
    Mt1 <- nrow(CH)
    t <- ncol(CH)
    ncov <- if (is.null(zi)) 0 else ncol(zi)
    nterms <- 1 + ncov + (Mb != 'None')
    if (!is.null(tj)) nterms <- nterms+1
    if (Mt == 'TimeCovariate') nterms <- nterms+1
    if (Mt == 'TimeSpecific') nterms <- nterms+(t-1)
    cij <- (apply(abs(CH), 1:2, sum) > 0) * 1
    zij <- t(apply(cij, 1, function(x) as.numeric(cumsum(x)>0)))
    zij <- cbind(rep(0,Mt1), zij[,-t])
    Mij <- cbind(rep(0,Mt1), cij[,1:(t-1)])
    start <- c(logit(mean(cij)*0.9), rep(0,nterms-1))
    fit <- suppressWarnings(optim (par = start, fn = negLLHuggins, hessian = TRUE, ...))
    covar <- try(MASS::ginv(fit$hessian))
    if (inherits(covar, 'try-error'))
        covar <- matrix(NA, nterms, nterms)
    list(capthist = CH, fit = fit, covar = covar, ncov = ncov, Mb = Mb, Mt = Mt, zi = zi, tj = tj)         
}

pibeta <- function (i, fit) {
    ## animal-specific probability of detection (animal i)
    beta <- fit$fit$par
    nocc <- ncol(fit$capthist)
    ncov <- fit$ncov
    p <- rep(beta[1], nocc)                  ## constant beta0
    if (ncov>0)
        p <- p + fit$zi[i,] %*% beta[2:(ncov+1)]  ## individual covariates (const)    
    if (fit$Mt == 'TimeSpecific') {
        p[1:(t-1)] <- p[1:(nocc-1)] + beta[ncov+(2:nocc)]
    }
    if (fit$Mt == 'TimeCovariate') {
        p <- p + beta[ncov+2] * fit$tj
    }
    p <- exp(p) / (1 + exp(p))
    1 - prod(1-p)
}

hugginsN <- function (CH) {
    zi <- NULL
    if (!is.null(covariates(CH))) {
        if ('distancetoedge' %in% names (covariates(CH)))
            zi <- matrix(covariates(CH)$distancetoedge, ncol=1)
    }
    fit <- huggins.fit(CH, zi = zi)
    
    phat <- sapply(1:nrow(CH), pibeta, fit=fit)
    Nhat <- sum (1 / phat)
    ## no variance yet s2 <- sum((1-phat) / phat^2)
    Nhat
}

hugginsp <- function (CH, dte = FALSE) {
    zi <- NULL
    if (!is.null(covariates(CH)) & dte) {
        if ('distancetoedge' %in% names (covariates(CH)))
            zi <- matrix(covariates(CH)$distancetoedge, ncol=1)
    }
    fit <- huggins.fit(CH, zi = zi)
    sapply(1:nrow(CH), pibeta, fit=fit)
}
