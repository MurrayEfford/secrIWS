## functions to return :
## Telem.Dhat -- IWS Telem estimate of density on inner zone
## estimatePG -- proportion overlap for each of n detected animals

estimatePG <- function (CH, NTelem, pTelem) {
    
    ## CH     -- capthist object with individual covariate 'truep' that should 
    ##           be non-missing for all animals
    ## NTelem -- number of telemetry fixes per animal
    ## pTelem -- proportion of animals telemetered

    n <- nrow(CH)
    
    ## number of telemetered animals
    ntel <- round(pTelem * n)
    ## select sample of animals to telemeter
    stel <- sample.int(n, ntel)
    ## binomial sample of in/out locations for telemetered animals only
    nG <- matrix(nrow = n, ncol = 2)
    nG[stel,1] <- rbinom(ntel, size = NTelem, prob = covariates(CH)$truep[stel])
    nG[,2] <- NTelem - nG[,1]    # number outside
    if (all(nG[stel,1] == 0)) {
        warning ("No telemetry fixes on grid")
    }
    else if (all(nG[stel,1] == 1)) {
        warning ("All telemetry fixes on grid")
    }
    else {
        ## logistic regression: proportion inside vs distance to edge
        CH <- addDistanceToEdge(CH)
        DTE <- covariates(CH)$distancetoedge
        model <- try(glm(nG ~ DTE, family = binomial, na.action = na.omit), silent = TRUE)
        PG <- rep(NA,n)
        if (inherits(model, 'try-error')) {
            ## attempt to trap error "Argument eta must be a nonempty numeric vector"
            warning ("glm failed")
        }
        else {
            ## proportion inside for each of n detected animals
            ## missing except for telemetered sample (stel)
            tmp <- try(predict(model, newdata = data.frame(DTE = DTE[-stel]), type = 'response'), silent = TRUE)
            if (!inherits(tmp, 'try-error')) {
                ## attempt to trap error "Argument eta must be a nonempty numeric vector"
                PG <- nG[,1] / NTelem   
                ## replace missing PG with logistic regression estimate
                PG[-stel] <- tmp
            }
        }
    }
    PG   ## vector of length n; does not remember which were telemetered
}

TelemDhat <- function (CH, PG, innerarea) {
    ## IWS estimator, not using DTE as covariate of p
    
    ## CH        -- capthist object 
    ## PG        -- vector of (estimated) proportion on grid, length n
    ## innerarea -- area of inner zone
    
    Telem.Dhat <- NA
    if (!any(is.na(PG))) {
        ## IWS estimator, not using DTE as covariate of p
        p <- hugginsp(CH, dte = FALSE)
        Telem.Dhat <- sum(PG/p) / innerarea
    }
    Telem.Dhat
}

