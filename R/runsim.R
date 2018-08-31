## potential development 2015-05-20

## Major
## hugginsN formula specification of model?

## Minor

## tweak secr::PG to allow rectangular (half trap spacing) buffer
## plotsim draw header info from plotted object
## runsim return object of class c('secrIWS', 'list')?
## summary.secrIWS method
## plot.secrIWS method
## cite IWS in secr PG.Rd

runsim <- function (nrepl, trapsites, template, N = 20, popnbuffer = 35, 
                    allinside = FALSE, gridaligned = FALSE, IWSlowerleft = FALSE, 
                    recentre = FALSE, rotateHR = FALSE, p1 = 0.4, noccasions = 7,
                    pTelem = 0.75, NTelem = 10, maskbuffer, maskspacing, 
                    extractfn = derived, plt = FALSE, SECR = TRUE, seed = 123,
                    ncores = 1, ... ) {

    ## preliminaries
    desc <- packageDescription("secrIWS")  ## for version number
    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    fit <- vector('list')
    spacing <- spacing(trapsites)
    if (missing(maskbuffer)) maskbuffer <- popnbuffer
    if (missing(maskspacing)) maskspacing <- spacing/3
    if (ncores>1) plt <- FALSE
    BVN <- !is.null(attr(template, 'sigmaX'))
    TELEM <- (NTelem>0) & (pTelem>0)
    if (TELEM & !(gridaligned | BVN)) {
        warning ("TELEM estimator only implemented for grid aligned or BVN HR")
        TELEM <- FALSE
    }
        
    ## mask is used only by secr::secr.fit()
    mask <- make.mask(trapsites, buffer = maskbuffer, spacing = maskspacing, type = 'traprect')

    ## outline of 'study site' (trapsites plus half-spacing strip)
    dxy <- c(-1,1) * spacing/2
    trappoly <- expand.grid(x = range(trapsites$x) + dxy,
                            y = range(trapsites$x) + dxy)[c(1,3,4,2,1),]
    innerarea <- polyarea(trappoly)

    onereplicate <- function(r) {
        ## simulate placement and shape of home ranges; assign probabilities
        pop <- placeanimals (N, trapsites, popnbuffer, template, Ndist = 'fixed',
                             allinside = allinside, gridaligned = gridaligned,
                             IWSlowerleft = IWSlowerleft, recentre = recentre, rotateHR = rotateHR)
        ## count centroids on 'study site'
        centresxy <- HRcentres(pop)
        N.inner <- sum(pointsInPolygon(centresxy, trappoly))
        ## plot and generate capture histories
        CH <- sim.HRCH(pop, trapsites, p1 = p1, noccasions = noccasions,
                       plt = plt, pltHR = plt, plttr = plt, pltfirst = TRUE)
        Telem.N <- attr(CH, 'trueN')        

        ## statistics computed directly on CH
        n <- nrow(CH)
        MMDM <- MMDM(CH)
        ## try to dodge 'finite coordinates' problem in chull within buffer.contour...
        MMDMarea <- try(sapply(buffer.contour(trapsites, buffer = MMDM, plt = FALSE, convex = TRUE), polyarea))
        if (inherits(MMDMarea, 'try-error')) MMDMarea <- NA
        MMDM2area <- try(sapply(buffer.contour(trapsites, buffer = MMDM/2, plt = FALSE, convex = TRUE), polyarea))
        if (inherits(MMDMarea, 'try-error')) MMDM2area <- NA
        
        ## Huggins closed population estimate for MMDM and MMDM2
        p <- hugginsp(CH, dte = FALSE)
        Nhat <- sum(1/p)
        Telem.Dhat <- Telem.n <- Dfit <- inner.Dhat <- NA
        
        ## estimate proportion of overlap (animal equivalents)
        PG <- estimatePG(CH, NTelem, pTelem)
        Telem.n <- sum(PG)

        ## IWS estimator
        if (TELEM & (n>0)) {
            Telem.Dhat <- TelemDhat (CH, PG, innerarea)
        }

        ## optionally fit SECR model and extract summary
        if (SECR & (n>0)) {
            fit <- try(secr.fit(CH, mask = mask, biasLimit = NA, verify = FALSE, ...))
            if (inherits(fit, 'try-error')) {
                warning ("secr model failed")
                fit <- NA; 
            }
            else {
                ## summarize model; must produce dataframe like predict.secr
                Dfit <- extractfn(fit)
                D <- Dfit['D','estimate']

                ## SECR realized Nhat on inner zone
                ## (pass full model fit to predictPG)
                w <- try(predictPG(fit, plt = FALSE, buffer = spacing/2))
                if (inherits(w, 'try-error')) w <- NA
                pd <- try(pdot(mask, traps(CH), detectfn = fit$detectfn, 
                           detectpar = detectpar(fit), noccasions = ncol(CH)))
                if (inherits(pd, 'try-error')) pd <- NA
                inner.Dhat <- (Telem.n + sum(w*(1-pd)) * D * attr(mask, 'area')) / innerarea
            }
        }
        if (nrepl>1) cat ('Completed replicate ', r, '\n')

        list(Nhat = Nhat, MMDM = MMDM, MMDMarea=MMDMarea, MMDM2area=MMDM2area, fit = Dfit, Telem.N = Telem.N, 
             Telem.n = Telem.n, Telem.Dhat = Telem.Dhat, N.inner = N.inner, Nr = nrow(pop), n = nrow(CH), 
             inner.Dhat = inner.Dhat)
    }
    ## loop over replicates
    if (ncores == 1) {
        set.seed(seed)        
        output <- lapply(1:nrepl, onereplicate)
    }
    else {
        ## use 'parallel'
        list(...)    ## ensures promises evaluated see parallel vignette 2015-02-02        
        clust <- makeCluster(ncores, outfile = 'logfile.txt')
        clusterSetRNGStream(clust, seed)                
        output <- parLapply (clust, 1:nrepl, onereplicate)
        stopCluster(clust)

    }
    ## output
    list(
         fit = lapply(output, '[[', 'fit'),
         Nr = sapply(output, '[[', 'Nr'),
         n = sapply(output, '[[', 'n'),
         mask.area = maskarea(mask),   ## aka 'arena'
         N.inner = sapply(output, '[[', 'N.inner'),
         Nhat = sapply(output, '[[', 'Nhat'),
         MMDM = sapply(output, '[[', 'MMDM'),
         MMDMarea = sapply(output, '[[', 'MMDMarea'),
         MMDM2area = sapply(output, '[[', 'MMDM2area'),
         Telem.N = sapply(output, '[[', 'Telem.N'),
         Telem.n = sapply(output, '[[', 'Telem.n'),
         Telem.Dhat = sapply(output, '[[', 'Telem.Dhat'),
         inner.Dhat = sapply(output, '[[', 'inner.Dhat'),
         N = N,
         inner.area = innerarea,
         pTelem = pTelem,
         NTelem = NTelem,
         popnbuffer = popnbuffer,
         allinside = allinside,
         gridaligned = gridaligned,
         IWSlowerleft = IWSlowerleft,
         recentre = recentre,
         p1 = p1,
         noccasions = noccasions,
         maskbuffer = maskbuffer,
         maskspacing = maskspacing,
         extractfn = extractfn,
         version = desc$Version,
         starttime = starttime,
         proctime = (proc.time() - ptm)[3])
}
