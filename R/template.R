build.irregular <- function (nc = 16, spacing = 10) {
    xy <- matrix(0, nrow=1, ncol=2)
    kernel <- matrix(c(-1,1,0,0, 0,0,-1,1), nrow = 4, ncol = 2)
    expand <- function(xy) sweep(kernel, FUN='+', STATS=xy, MARGIN=2)
    i <- 1
    while (i<nc) {
        tmp <- matrix(0, nrow=0, ncol=2)
        connected <- apply(xy,1,expand)
        connected <- matrix(t(connected), ncol=2)         
        #            connected <- do.call(rbind,connected)
        tmp <- unique(rbind(xy,connected))
        tmp <- tmp[-(1:i),,drop=FALSE]
        ni <- nrow(tmp)
        tmp <- tmp[sample.int(ni,1),,drop=FALSE]
        xy <- rbind(xy,tmp)
        i <- i+1
    }
    xy <- xy * spacing
    dimnames(xy) <- list(1:nc,c('x','y'))    
    read.mask(data = data.frame(xy), spacing = spacing)    
}

templatefn <- function (nx = 4, ny = 4, spacing = 10, sigmaX = NULL, sigmaY = NULL) {
    BVN <- !is.null(sigmaX)
    if (BVN & is.null(sigmaY)) sigmaY <- sigmaX
    if (BVN) {
        if ((nx!=1) | (ny!=1))
            warning("BVN template setting nx = 1, ny = 1 (centre point)")
        nx <- ny <- 1
    }
    template <- expand.grid(x = seq(-(nx-1)/2, (nx-1)/2, 1)*spacing, 
                            y = seq(-(ny-1)/2, (ny-1)/2, 1)*spacing)
    temp <- read.mask(data = template, spacing = spacing) 
    if (BVN) {
        attr(temp, 'sigmaX') <- sigmaX
        attr(temp, 'sigmaY') <- sigmaY
    }
    temp
}

#########################################
