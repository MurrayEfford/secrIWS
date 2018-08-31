## 2016-01-24
## assume single session for now
## only allow square region aligned with grid
## use mask in object
## For each mask point, what proportion of a circular bivariate normal density
## lies within the buffered square bounding the traps?

predictPG <- function (object, buffer = 0, plt = FALSE, ...) {
    if (!(object$detectfn %in% c(0,14)))
        stop ("predictPG at present requires halfnormal detectfn ('HN','HHN')")
    X <- object$mask
    x <- range(traps(object$capthist)[,1]) + c(-1,1) * buffer
    y <- range(traps(object$capthist)[,2]) + c(-1,1) * buffer
    dp <- detectpar(object)
    cov <- diag(rep(dp$sigma^2,2))
    out <- apply(X,1, mvtnorm::pmvnorm, lower = c(x[1],y[1]), upper = c(x[2],y[2]), sigma = cov)
    if (plt) {
        covariates(X) <- data.frame(predictPG = out)
        plot (X, covariate = 'predictPG', ...)
        polygon(cbind(x[c(1,1,2,2)], y[c(1,2,2,1)]))
    }
    out
}


## following code is exploratory only -- not used in simulations

## expected value of first term in Efford MS eqn 2
term1 <- function (object, buffer = 0) {
    D <- predict(object)['D','estimate']
    w <- predictPG(object, plt = FALSE, buffer = buffer)
    pd <- pdot(object$mask, traps(object$capthist), detectfn = object$detectfn, 
               detectpar=detectpar(object), noccasions = ncol(object$capthist))
    sum(w*pd) * D * attr(object$mask, 'area')
}

## second term in Efford MS eqn 2
term2 <- function (object, buffer = 0) {
    D <- predict(object)['D','estimate']
    w <- predictPG(object, plt = FALSE, buffer = buffer)
    pd <- pdot(object$mask, traps(object$capthist), detectfn = object$detectfn, 
               detectpar=detectpar(object), noccasions = ncol(object$capthist))
    sum(w*(1-pd)) * D * attr(object$mask, 'area')
}

## using secrIWS 1.0.2 2016-01-24
## requires full likelihood model fit
innerN <- function(iwsobject, scennum = 1, buffer = 0) {
    if (class(iwsobject$fit[[1]]) != 'secr')
        stop ("require output from secrIWS 'runsim' with extractfn = identity")
    trueN <-iwsobject[[scennum]]$Telem.N
    t1 <- iwsobject[[scennum]]$Telem.n
    t2 <- sapply(iwsobject[[scennum]]$fit, term2, buffer)
    list(innerN = t1+t2, trueN = trueN)
}

#############################################################################

## fit0 <- secr.fit(captdata, detectfn = 14)
## pPG <- predictPG(fit0, plt=T, dots=F, buffer = 30)
## plot(traps(captdata), add=T)

# (term1(fit0, 0) + term2(fit0, 0)) / 7.29
# [1] 5.484453
# predict(fit0)['D',]
# link estimate SE.estimate      lcl     ucl
# D  log 5.484609   0.6470283 4.355869 6.90584