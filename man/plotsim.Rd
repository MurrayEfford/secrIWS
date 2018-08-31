\name{plotsim}
\alias{plotsim}
\alias{plotTelem}

\title{Plot Simulation Output}

\description{

  These functions reproduce reasonably closely the plots of Ivan et al. (2013b) Figs 2 and 3.

}

\usage{

plotsim(simout, add = FALSE, estimator = c("SECR", "TELEM", "MMDM2", "MMDM", "SECRT"), 
    type = c("inner","total", "telem"), p = seq(0, 1, 0.01), ylim = c(-100,100),
    xlab = 'Percentile', ylab = 'quantile(PE)', label = 'xy', ...)

}

\arguments{
  \item{simout}{list object output by \code{\link{runsim}}}
  \item{add}{logical; if TRUE a curve is added to an existing plot}
  \item{estimator}{character value coding ML SECR vs TELEM (Ivan et al 2013a)}
  \item{type}{character value coding the reference density}
  \item{p}{numeric vector of quantile probabilities}
  \item{ylim}{limits on y-axis}
  \item{xlab}{character label for x-axis}
  \item{ylab}{character label for y-axis}
  \item{label}{character string indicating which axes to be labelled}
  \item{\dots}{other arguments passed to plot() or lines()}
}
\details{
 
The default estimator is "SECR".

\code{plotsim} plots quantiles of percentage error of density estimates computed relative to the number of home range centres in the arena ("total"), the number of centres on the trapping grid itself (`inner'), or the summed home-range utilization of the trapping grid. The area of the trapping grid includes a half-trap-spacing marginal strip. The default is "inner".
  
}

\value{
  
Invisibly returns a dataframe with the estimated density (Dhat), total density (arenaD), and inner density (innerD) for each replicate.
 
}

\references{

Ivan, J. S., White, G. C. and Shenk, T. M. (2013a) Using auxiliary telemetry information to 
estimate animal density from capture--recapture data. \emph{Ecology} \bold{94}, 809--816.

Ivan, J. S., White, G. C. and Shenk, T. M. (2013b) Using simulation to compare methods for 
estimating density from capture--recapture data. \emph{Ecology} \bold{94}, 817--826.

}

\seealso{
  
  \code{\link{runsim}}
  
}

\examples{

\dontrun{
template44 <- templatefn(nx=4, ny=4, spacing = 10)
grid10 <- make.grid(10,10, spacing = 10)
out1 <- runsim(nrepl = 100, trapsites = grid10, template = template44, 
               CL = TRUE, start = list(g0 = 0.1, sigma = 10))
plotsim(out1)
}

}

\keyword{ hplot }
