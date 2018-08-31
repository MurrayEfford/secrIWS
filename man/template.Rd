\name{templatefn}
\alias{templatefn}
\alias{build.irregular}

\title{Construct Home Range Template or Instance}

\description{

The shape of home ranges in Ivan et al. (2013) was either constant or
`irregular' (constructed at random according to some rules). These cases
are handled by the functions \code{templatefn} and
\code{build.irregular} respectively. 

}

\usage{
templatefn(nx, ny, spacing = 10, sigmaX = NULL, sigmaY = NULL)
build.irregular(nc = 16, spacing = 10)

}
\arguments{
  \item{nx}{number of cells in x direction}
  \item{ny}{number of cells in y direction}
  \item{spacing}{cell width}
  \item{sigmaX}{numeric dispersion of BVN home range in x-direction}
  \item{sigmaY}{numeric dispersion of BVN home range in y-direction}  
  \item{nc}{number of cells}
}

\details{
  
  \code{build.irregular} generates one random instance; usually each
  individual is expected to have a uniquely generated HR, and for this
  build.irregular must be passed to runsim.
  
  If \code{sigmaX} is not NULL a bivariate normal (BVN) range is assumed. Only the centre point is used (nx = 1, ny = 1). 
  
}

\value{

  A mask object compatible with `secr', particularly the plot method for
  mask objects. 
  
  For BVN ranges the result is a token 1-cell mask; the values of sigmaX and sigmaY are saved as attributes.
  
}

\examples{
plot(templatefn(3,3), dots = FALSE, meshcol='blue')
plot(build.irregular(), dots = FALSE, meshcol='blue')

}
\keyword{ datagen }
