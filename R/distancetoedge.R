## from secrlinear
make.sldf <- function (coord, f) {
    ## Form SpatialLinesDataFrame from coordinates data
    ## coord is dataframe of coordinates - must include columns 'x','y'
    ## f is vector of values by which to split coord rows
    if (missing(f) & ('LineID' %in% names(coord)))
        f <-  coord$LineID
    if (missing(f))
        coordlist <- list(coord)
    else
        coordlist <- split(coord[,c('x','y')], f)
    S0 <- lapply(coordlist, Line)
    S1 <- Lines(S0, ID = '1')
    S2 <- SpatialLines(list(S1))
    ldf <-  data.frame(ID = 1:length(S2), rownames=1)
    SpatialLinesDataFrame(S2, data = ldf)
}

addDistanceToEdge <- function (CH) {
    ## Only for traps, not proximity detectors yet 
    distancetoedge <- function (CH1, trapsites, trappoly) {
        xy <- trapsites[abs(CH1),]   ## drops zeros automatically
        xy <- matrix(apply(xy,2,mean), ncol = 2)
        SP <- sp::SpatialPoints(xy)
        d <- rgeos::gDistance(SP, Sldf)
        d
    }
    if (nrow(CH) > 0) {
        trapsites <- traps(CH)
        trappoly <- as.data.frame(trapsites[chull(trapsites),])
        trappoly <- rbind(trappoly, trappoly[1,])
        Sldf <- make.sldf(trappoly, rep(1,nrow(trappoly)))
        d <- apply(CH, 1, distancetoedge, trapsites, trappoly = Sldf)
        if (is.null(covariates(CH)))
            covariates(CH) <- data.frame(distancetoedge = d)
        else {
            if (nrow(covariates(CH))==0)
                covariates(CH) <- data.frame(distancetoedge = d)
            else
                covariates(CH)$distancetoedge <- d
        }
    }
    CH
}
