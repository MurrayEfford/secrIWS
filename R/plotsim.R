plotsim <- function (simout, add = FALSE, estimator = c('SECR','TELEM','MMDM2','MMDM', 'SECRT'),
                     type = c('inner', 'total', 'telem'), p = seq(0,1,0.01), 
                     ylim = c(-100,100), xlab = 'Percentile', ylab = 'quantile(PE)', 
                     label = 'xy', ...) {
    ## assumes extractfn has returned redicted or derived values
    estimator <- match.arg(estimator)
    type <- match.arg(type)
    D <- switch(estimator, 
                SECR = sapply(simout$fit,'[[','D','estimate'),
                TELEM = simout$Telem.Dhat, 
                MMDM2 = simout$Nhat/as.numeric(simout$MMDM2area),
                MMDM  = simout$Nhat/as.numeric(simout$MMDMarea),
                SECRT = simout$inner.Dhat
    )
    nNA <- sum(is.na(D))
    if (nNA>0) warning(nNA, " estimate(s) could not be calculated")
    arenaD <- simout$Nr / simout$mask.area
    innerD <- simout$N.inner / simout$inner.area
    telemD <- simout$Telem.N / simout$inner.area
    tD <- switch (type, total = arenaD, inner = innerD, telem = telemD)
    PE <- (D - tD)/ tD * 100
    if (!add) {
        if (!(label %in% c('x','xy'))) xlab <- ''
        if (!(label %in% c('y','xy'))) ylab <- ''
        plot(p, quantile(PE, p, na.rm = TRUE), type='l', ylim = ylim, 
             xaxs = 'i', yaxs = 'i', las = 1, axes = FALSE, xlab='', ylab='',...)
        if (label %in% c('x','xy')) {
            axis(1, at = seq(0,1,0.2))
            mtext (side = 1, line = 2.2, xlab)
        }
        if (label %in% c('y','xy')) {
            axis(2, at=seq(-100,100,20), las = 1)
            mtext (side = 2, line = 2.5, ylab)
        }
        box()
        abline(h=0, col = 'grey')
        abline(v=0.5, col = 'grey')
    }
    else
        lines(p, quantile(PE, p, na.rm = TRUE), ...)
    invisible(data.frame(Dhat = D, arenaD = arenaD, innerD = innerD))
}
#########################################

