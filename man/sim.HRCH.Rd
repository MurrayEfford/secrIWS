\name{sim.HRCH}
\alias{sim.HRCH}

\title{Generate Spatial Capture--Recapture Data for Simulated Population}

\description{

This function takes a simulated distribution of individuals and a trap
array and generates a set of capture histories. It is used internally by 
\code{\link{runsim}} and will not usually be called by the user.

}
\usage{

sim.HRCH(popn, trapsites, p1, noccasions, plt = TRUE, pltHR = TRUE, plttr = TRUE,
    pltfirst = FALSE)

}

\arguments{
  \item{popn}{secr popn object}
  \item{trapsites}{secr traps object}
  \item{p1}{numeric overall probability of capture}
  \item{noccasions}{integer number of sampling occasions}
  \item{plt}{logical; if TRUE the population is plotted}
  \item{pltHR}{logical; if TRUE home ranges are added to the plot}
  \item{plttr}{logical; if TRUE traps are added to the plot}
  \item{pltfirst}{logical; if TRUE a mark is placed on the first cell of each plotted home range}
}

\details{
  
This function combines some plotting capability with simulation of
capture-recapture sampling. Unlike the `secr' function sim.capthist,
there is no need to specify a detection function: the cellwise detection
probabilities are governed by probabilities already computed and stored
in an attribute of \code{popn}.

}

\value{

A 3-D secr capthist object (n x S x K).
  
}

\seealso{

  \code{\link{runsim}}
  
}

\examples{

opar <- par(xpd = TRUE, mar = c(5,5,5,5))
on.exit(par(opar))
grid10 <- make.grid(nx = 10, ny = 10, spacing = 10)
template44 <- templatefn(4, 4)
pop <- placeanimals (20, grid10, 35, template44)
CH <- sim.HRCH(pop, grid10, p1 = 0.2, noccasions = 7)
plot(CH, add = TRUE, varycol = FALSE)

}

\keyword{ datagen }
