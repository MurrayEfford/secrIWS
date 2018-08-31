\name{secrIWS-package}
\alias{secrIWS}

\docType{package}

\title{Evaluation of Spatially Explicit Capture-Recapture with Irregular Data}

\description{

A set of functions used to reproduce and extend the simulations of Ivan et al. (2013b) for evaluating the estimator of Ivan et al. (2013a) in relation to spatially explicit capture--recapture estimators.
  
}

\details{
\tabular{ll}{
Package: \tab secrIWS\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2018-09-01\cr
License: \tab GNU General Public License Version 2 or later\cr
}

The important functions in \pkg{secrIWS} are:
\tabular{ll}{
  \code{\link{runsim}} \tab run replicate simulations of a particular scenario from Ivan et al. (2013b)\cr
  \code{\link{templatefn}} \tab function to generate home ranges with simple geometrical shapes \cr
  \code{\link{build.irregular}} \tab function to generate irregular home ranges on-the-fly \cr
  \code{\link{plotsim}} \tab generate quantile plots \cr
}

Other useful exported functions are:
\tabular{ll}{
  \code{\link{plotmaskedge}} \tab plots the perimeter of a \pkg{secr} habitat mask \cr
  \code{\link{addDistanceToEdge}} \tab add covariate to a capthist object \cr
  \code{\link{hugginsN}} \tab Higgins conditional likelihood estimator of N \cr
  \code{\link{predictPG}} \tab overlap of bivariate circular normal distribution with square region \cr
}

Other internal functions are:
\tabular{ll}{
  \code{\link{sim.HRCH}} \tab simulate IWS home ranges and capture histories \cr
}

Documentation will be provided in a vignette
\url{../doc/secrIWS-vignette.pdf} and in the pdf version of the help
pages \url{../doc/secrIWS-manual.pdf}

}

\seealso{

  \code{\link{runsim}}
  
}

\references{

Ivan, J. S., White, G. C. and Shenk, T. M. (2013a) Using auxiliary telemetry information to 
estimate animal density from capture--recapture data. \emph{Ecology} \bold{94}, 809--816.

Ivan, J. S., White, G. C. and Shenk, T. M. (2013b) Using simulation to compare methods for 
estimating density from capture--recapture data. \emph{Ecology} \bold{94}, 817--826.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\examples{

template44 <- templatefn(nx = 4, ny = 4, spacing = 10)
grid10 <- make.grid(10,10, spacing = 10)
out1 <- runsim(nrepl = 20, trapsites = grid10, template = template44, 
               CL = TRUE, start = list(g0 = 0.1, sigma = 10))

}

\keyword{ package }
