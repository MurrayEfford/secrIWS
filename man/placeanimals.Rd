\name{placeanimals}
\alias{placeanimals}

\title{Simulate Population and Home Ranges}

\description{
  
  This function places home ranges and assigns cell-specific detection
  probabilities. Home ranges may lie entirely within the rectangular
  `arena' or be allowed to straddle the edge. It is used internally by 
  \code{\link{runsim}} and will not usually be called by the user.
  
}
\usage{

placeanimals(N, trapsites, buffer, template, Ndist = c("fixed", "poisson"),
allinside = FALSE, gridaligned = FALSE, IWSlowerleft = FALSE, recentre = FALSE, 
rotateHR = FALSE)

}

\arguments{
  \item{N}{integer number of animals to simulate}
  \item{trapsites}{secr traps object}
  \item{buffer}{numeric width of buffer zone between outermost trap and
  edge of arena}
  \item{template}{secr mask object describing one home range, or
  function generating same}
  \item{Ndist}{character distribution of N}
  \item{allinside}{logical; if TRUE home ranges must lie entirely inside
    the arena}
  \item{gridaligned}{logical; if TRUE each simulated home range comprises
  only cells}
  \item{IWSlowerleft}{logical; if TRUE the IWS algorithm is used for placing grid-cell home ranges} 
  \item{recentre}{logical; if TRUE each irregular home range is positioned by its centroid}
  \item{rotateHR}{ logical: if TRUE then template is randomly rotated (see Details)}
}

\details{

The default placement of irregular home ranges centres on the first
cell, the one on which the remaining home range is built.
  
Each home range is randomly rotated when \code{rotateHR = TRUE} and
\code{template} is not a function. When \code{gridaligned = TRUE}
rotation is a multiple of 90 degrees.  }

\value{

  `popn' object with attribute `HRlist'. Each component of HRlist is a
  small mask object representing the cells of a home range, with
  attribute `p' for cell-specific capture probabilities.

}

\seealso{
  
 \code{\link{templatefn}}, \code{\link{build.irregular}}, \code{\link{runsim}}
  
}

\examples{

template44 <- templatefn(nx = 4, ny = 4, spacing = 10)
grid10 <- make.grid(nx = 10, ny = 10, spacing = 10)
pop <- placeanimals(N = 20, trapsites = grid10, buffer = 35, template44, 
    allinside = TRUE)
plot(pop)
plot(grid10, add = TRUE)

}

\keyword{ datagen }
