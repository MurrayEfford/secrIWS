\name{runsim}
\alias{runsim}

\title{Simulated Trapping and Estimation}

\description{
  
This function runs a set of simulations. In each simulation a new
population is generated, a capture--recapture sample is taken, and
population density is estimated.

}

\usage{

runsim (nrepl, trapsites, template, N = 20, popnbuffer = 35, allinside = FALSE,
gridaligned = FALSE, IWSlowerleft = FALSE, recentre = FALSE, rotateHR = FALSE, 
p1 = 0.4, noccasions = 7, pTelem = 0.75, NTelem = 10, maskbuffer, maskspacing,
extractfn = derived, plt = FALSE, SECR = TRUE, seed = 123, ncores = 1, ... )

}

\arguments{
  \item{nrepl}{integer number of replicates}
  \item{trapsites}{traps object from `secr'}
  \item{template}{either a (small) secr mask object or a function to
  generate such an object}
  \item{N}{integer number of individuals to place in the arena}
  \item{popnbuffer}{numeric width of buffer between outermost traps and edge
  of the arena}
  \item{allinside}{logical; if TRUE each simulated home range must lie entirely inside the arena}
  \item{gridaligned}{logical; if TRUE each simulated home range comprises
  only cells}
  \item{IWSlowerleft}{logical; if TRUE the IWS algorithm is used for placing home ranges}
  \item{recentre}{logical; if TRUE each irregular home range is positioned by its centroid rather than the first cell}
  \item{rotateHR}{ logical: if TRUE then template is randomly rotated (see Details in 
    \code{\link{placeanimals}})}
  \item{p1}{numeric overall probability of capture}
  \item{noccasions}{integer number of sampling occasions}
  \item{pTelem}{proportion of animals telemetered}
  \item{NTelem}{number of fixes per telemetered animal}
  \item{maskbuffer}{numeric buffer used in \code{\link{secr.fit}}}
  \item{maskspacing}{numeric spacing of mask used by secr.fit}
  \item{extractfn}{ function used to summarise each simulated model fit}
  \item{plt}{logical; if TRUE then each replicate will be plotted as it
    is generated}
  \item{SECR}{logical; if TRUE then a spatially explicit capture-recapture model will be fitted}
  \item{seed}{integer random number seed}
  \item{ncores}{integer number of cores to use (ncores > 1 for parallel processing)}
  \item{\dots}{other arguments passed to \code{\link{secr.fit}}}
}

\details{

\code{trapsites} and \code{popnbuffer} define an `arena' in which simulated animals are placed. By default, home ranges may extend outside the arena (\code{allinside = FALSE}).

\code{maskbuffer} defaults to the same value as \code{popnbuffer}.
  
\code{maskspacing} defaults to \code{spacing(trapsites) / 3}.

Areas are specified in hectares and linear measurements in metres.

Bivariate normal home range probabilities are computed if a value is
supplied for \code{sigmaX}. \code{sigmaY} defaults to the same value as
\code{sigmaX} (i.e., a circular home range).

Plotting is suppressed if processing is distributed across multiple cores (\code{ncores} > 1). Plotting uses function \code{\link{sim.HRCH}}; the arguments \code{pltHR}, \code{plttr} and \code{pltfirst} of that function are all set to the value of \code{plt} in \code{runsim}.

\code{extractfn} should take a fitted secr model as its first argument and return a single dataframe with a row labelled `D' and a column labelled `estimate'.

A spatially explicit capture-recapture model is fitted by maximizing the likelihood in \code{\link{secr.fit}}. The actual model may be controlled by providing named arguments in the `dots' argument of \code{runsim}. Care must be taken to ensure \code{extractfn} matches the output expected. The default value of \code{extractfn} is `derived' which works with either CL = TRUE or CL = FALSE (the former is faster). 

The IWS algorithm for placing home ranges is selected by setting \code{IWSlowerleft = TRUE}. Flaws in this method are identified in the vignette.

}

\value{
  A list with components
  \item{fit}{list of output from extractfn (default: derived(fit))}
  \item{Nr}{vector length nrepl (total number of individuals simulated (realised))}
  \item{n}{vector length nrepl (number of individuals caught)}
  \item{mask.area}{area of simulation arena (ha)}
  \item{N.inner}{number of simulated individuals whose centroids lay
  within `inner' zone (rectangular buffer extending half spacing beyond
  outermost traps) }
  \item{Telem.N}{true sum of home-range fractions on inner zone}
  \item{Telem.n}{sum of estimated home-range fractions on inner zone for detected individuals }
  \item{Telem.Dhat}{density of inner zone estimated by TELEM method}
  
  and these echoed inputs --  
  \item{N}{}
  \item{inner.area}{}
  \item{pTelem}{}
  \item{NTelem}{}
  \item{popnbuffer}{}
  \item{allinside}{}
  \item{gridaligned}{}
  \item{recentre}{}
  \item{p1}{} 
  \item{noccasions}{}
  \item{maskbuffer}{}
  \item{maskspacing}{}
  \item{extractfn}{}
  
  and runtime information --
  \item{version}{}
  \item{starttime}{}
  \item{proctime}{}

}

\references{

 Ivan, J. S., White, G. C. and Shenk, T. M. (2013) Using simulation to
  compare methods for estimating density from capture--recapture
  data. \emph{Ecology} \bold{94}, 817--826.

}

\seealso{
  
  \code{\link{sim.HRCH}},
  \code{\link{secr.fit}}
  
}

\section{Warning}{

The code can fail unpredictably with the message "finite coordinates are needed". The source of this error has not been traced; the solution for now is to re-run the simulations.

}

\examples{

\dontrun{
template44 <- templatefn(x = seq(-15,15,10), y = seq(-15,15,10), spacing = 10)
grid10 <- make.grid(10,10, spacing = 10)
out1 <- runsim(nrepl = 20, trapsites = grid10, template = template44, 
               CL = TRUE, start = list(g0 = 0.1, sigma = 10))
}
}
\keyword{ datagen }
