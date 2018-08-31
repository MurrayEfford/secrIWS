Px <- function(trapsites, template) {
    BVN <- !is.null(attr(template, 'sigmaX'))
    if (BVN) {
        sigmaX <- attr(template, 'sigmaX')
        sigmaY <- attr(template, 'sigmaY')
        spc <- spacing(trapsites)
        ## for each cell, compute BVN probability
        p <- numeric(nrow(trapsites))
        low <- apply(trapsites, 1, '-', spc/2)
        upp <- apply(trapsites, 1, '+', spc/2)
        cov <- diag(c(sigmaX, sigmaY)^2)
        centre <- unlist(template)  ## just x,y for one point
        for (j in 1:nrow(trapsites)) {
            p[j] <- mvtnorm::pmvnorm (lower = low[,j], upper = upp[,j],
                                      mean = centre, sigma = cov)
        }
        tol <- 1e-8  ## to catch integration errors?
        p[p<tol] <- 0
    }
    else{
        x <- suppressWarnings(addCovariates(trapsites, template, 'p', strict = T))  ##
        p <- covariates(x)$p
        p[is.na(p)] <- 0
    }
    p
}

## build one capture history
rCH <- function (px, p1, noccasions, hazard = FALSE) {
    ntrap <- length(px)
    pcapt <- sum(px, na.rm = TRUE)
    if (hazard) pcapt <- 1-exp(-pcapt)
    pcapt <- pcapt * p1
    if (pcapt<1e-8)
        rep(0,noccasions)
    else {
        prob <- c(1-pcapt, pcapt * px/sum(px))
        sample.int(ntrap+1, noccasions, prob = prob, replace = TRUE) - 1
    }
}

HRcentres <- function(pop) {
    getone <- function (mask) {
        wts <- covariates(mask)$p
        c(sum(mask$x*wts), sum(mask$y*wts))
    }
    HRlist <- attr(pop, 'HRlist')
    BVN <- !is.null(attr(HRlist[[1]], 'sigmaX'))
    if (BVN) ## symmetrical, no change
        as.matrix(pop)
    else
        t(sapply(HRlist, getone))
}

ellipse <- function (center, shape, radius, segments = 51) {
    ## from ellipse function in 'car' package
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    Q <- chol(shape, pivot = TRUE)
    order <- order(attr(Q, "pivot"))
    ellipse <- t(center + radius * t(unit.circle %*% Q[, order]))
    colnames(ellipse) <- c("x", "y")
    ellipse
}

sim.HRCH <- function (popn, trapsites, p1, noccasions, plt = TRUE, pltHR = TRUE, plttr = TRUE,
                      pltfirst = FALSE) {
    ## for each animal determine its location, home range and CH
    ## assume independence between animals
    np <- nrow(popn)
    HRlist <- attr(popn, 'HRlist')
    BVN <- !is.null(attr(HRlist[[1]], 'sigmaX'))
    spc2 <- spacing(trapsites) / 2
    trappoly <- expand.grid(x = range(trapsites$x)+c(-1,+1)*spc2,
                            y = range(trapsites$x)+c(-1,+1)*spc2)[c(1,3,4,2,1),]
    if (plt) plot(popn)
    if (pltHR) {

        if (BVN) {
            for (i in 1:np) {
                sigmaX <- attr(HRlist[[i]], 'sigmaX')
                sigmaY <- attr(HRlist[[i]], 'sigmaY')
                if (is.null(sigmaY)) sigmaY <- sigmaX
                cov <- diag(c(sigmaX, sigmaY)^2)
                centre <- unlist(popn[i,])
                ellipsi <- ellipse(centre, shape = cov, radius = 2.45)
                polygon (ellipsi, lwd=1.5)
            }
        }
        else {
            lapply (HRlist, plot, add=T, dots=F, cov='p',
                col=rev(grey(seq(0,0.99,0.06))), legend = F)
            if (pltfirst)
                lapply(HRlist, function(x) points(x[1,], pch=4))
            lapply (HRlist, plotmaskedge, col='black')
        }
        centresxy <- HRcentres(popn)
        points(centresxy, col='green',pch=16)
        polygon(trappoly, border='red')
    }
    out <- matrix(0, np, noccasions)
    px <- numeric(np)
    ## loop over animals
    for (i in 1:np) {
        pxi <- Px(trapsites, HRlist[[i]])
        out[i,] <- rCH(pxi, noccasions, p1 = p1)
        px[i] <- sum(pxi) ## sum over 'trap' cells
    }
    if (plttr) {
        plot(trapsites, add = TRUE)
        polygon(trappoly, border='red')
    }
    nonzeroCH <- apply(out,1,sum)>0
    out <- out[nonzeroCH,,drop=FALSE]  ## drop animals not caught
    class(out) <- 'capthist'
    traps(out) <- trapsites
    covariates(out) <- data.frame(truep = px[nonzeroCH])
    
    ## 2018-07-11
    ## convert old 2-D format to 3-D format
    rownames(out) <- 1:nrow(out)
    out <- updateCH(out)
    
    attr(out, 'trueN') <- sum(px)  ## remember sum over whole population
    out
}
