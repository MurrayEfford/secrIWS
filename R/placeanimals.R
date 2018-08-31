placeanimals <- function (N, trapsites, buffer, template, Ndist = c('fixed', 'poisson'),
                          allinside = FALSE, gridaligned = FALSE, IWSlowerleft = FALSE,
                          recentre = FALSE, rotateHR = FALSE) {
    
    vertinside <- function (HR, arena) {
        if (BVN) {   ## 95% contour
            centre <- unlist(HR)
            vertices <- ellipse(centre, shape = cov, radius = 2.45)
        }
        else {
            corners <- cbind(c(-1,-1,+1,+1), c(-1,+1,+1,-1)) * spacing(HR)/2 * (1-tol)
            vertices <- apply(HR, 1, function(x) sweep(corners, FUN='+', STATS=x, MARGIN=2))
            vertices <- unique(vertices)
        }
        all(pointsInPolygon(vertices, arena))
    }
    tol <- 1e-6
    Ndist <- match.arg(Ndist)
    if (Ndist == 'poisson')
        N <- rpois (1, N)
    #     else
    #         N <- secr:::discreteN(N)
    if (N < 1)
        stop ("attempt to simulate zero animals")
    spc <- spacing(trapsites)
    if (gridaligned) {
        if (IWSlowerleft)
            cells <- spc * as.matrix(expand.grid(x = -3:9, y=-3:9))
        else
            cells <- as.matrix(make.mask(trapsites, buffer = buffer, spacing = spc))
    }
    HRlist <- vector('list', N)
    pop <- matrix(nrow = N, ncol = 2)
    maxtries <- 1000
    
    xrange <- range(trapsites$x) + c(-1,+1)*buffer
    yrange <- range(trapsites$y) + c(-1,+1)*buffer
    arena <- expand.grid(x = xrange, y = yrange)[c(1,3,4,2,1),]  ## before we mess with ranges
    
    ## 2015-06-01 adjust default BVN placement to 0-150 from 0-160 on each axis
    ## i.e. block placement in 1-cell strip on right and top sides
    ## change should not affect any other placements
#     if (IWSlowerleft) {  
#         xrange[2] <- xrange[2] - spc
#         yrange[2] <- yrange[2] - spc
#     }
    if (IWSlowerleft) {  
        xrange <- xrange +  c(spc/2, -spc/2)
        yrange <- yrange +  c(spc/2, -spc/2)
    }
    dx <- diff(xrange)
    dy <- diff(yrange)
    origin <- c(xrange[1], yrange[1])
    adjust <- c(0,0)

    ## prepare for BVN home ranges
    BVN <- !is.null(attr(template, 'sigmaX'))
    if (BVN) {
        sigmaX <- attr(template, 'sigmaX')
        sigmaY <- attr(template, 'sigmaY')
        if (is.null(sigmaY)) sigmaY <- sigmaX
        cov <- diag(c(sigmaX, sigmaY)^2)
    }
    
    ## loop over all animals
    for (i in 1:N) {
        if (is.function(template))
            temp <- template()
        else {
            temp <- template
            if (gridaligned & !IWSlowerleft) {
                if (rotateHR) {
                    temp[,] <- secr::rotate(temp, sample.int(4,1) * 90)
                }                
                adjustx <- ((diff(range(temp[,1]))/spc) %% 2) * spc/2 * sample(c(-1,1), 1)
                adjusty <- ((diff(range(temp[,2]))/spc) %% 2) * spc/2 * sample(c(-1,1), 1)
                adjust <- c(adjustx, adjusty)
            }
            else {
                if (rotateHR)
                    temp[,] <- secr::rotate(temp, runif(1)*360)
                adjust <- c(0,0)
            }
        }
        tries <- 0
        repeat {
            if (gridaligned) {
                centre <- cells[sample.int(nrow(cells),1),] + adjust
            }
            else
                centre <- runif(2)* c(dx,dy) + origin
            HR <- temp
            p <- runif(nrow(HR))
            p <- p/sum(p)
            covariates(HR) <- data.frame(p = p)
            if (recentre) {
                ## centre on centroid
                if (gridaligned) stop("recentre option incompatible with gridaligned option")
                centre <- centre - c(sum(HR$x*p), sum(HR$y*p))
            }
            if (IWSlowerleft & !(is.function(template))) {   ## zero at lower left
                HR[,] <- sweep(HR, STATS = apply(HR,2,min), MARGIN=2, FUN = '-' )
            }                
            HR[,] <- round(sweep(HR, STATS = centre , MARGIN=2, FUN = '+' ),6)
            OK <- if(allinside) vertinside(HR, arena) else TRUE
            if (OK) {
                HRlist[[i]] <- HR         
                pop[i,] <- centre
                break
            }
            tries <- tries+1
            if (tries > maxtries) break

        }
        if (tries > maxtries) stop("failed to place HR in ", maxtries, " attempts")
    }
    ## return this as a 'popn' object with extra attribute HRlist
    pop <- as.data.frame(pop)
    attr(pop, 'HRlist') <- HRlist
    attr(pop, 'boundingbox') <- arena[1:4,]
    attr(pop, 'Ndist') <- Ndist
    attr(pop, 'model2D') <- 'user'
    attr(pop, 'buffertype') <- 'rect'
    dimnames(pop) <- list(1:N, c('x','y'))
    class(pop) <- c('popn', 'data.frame')
    pop
}