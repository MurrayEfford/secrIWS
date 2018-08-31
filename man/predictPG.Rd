\name{predictPG}
\alias{predictPG}

\title{Predicted Activity on Grid}

\description{For each mask point of a fitted `secr' model, what proportion of a circular bivariate normal density lies within the buffered square bounding the detectors?
}

\usage{
predictPG(object, buffer = 0, plt = FALSE, ...)
}

\arguments{
  \item{object}{single-session 'secr' model from \code{\link{secr.fit}}}
  \item{buffer}{numeric distance from outermost detector to buffered edge}
  \item{plt}{logical; if TRUE a plot is generated}
  \item{\dots}{other arguments passed to \code{\link{plot.mask}}}
}

\details{

A square detector array is assumed. It must be aligned with the axes. The buffered corners are not rounded.

The model object may have been fitted by maximizing the full likelihood (\code{CL = TRUE}) or the conditional likelihood (\code{CL = TRUE}).

The present algorithm is restricted to models fitted with half-normal detection functions because it uses the multivariate normal density function \code{\link{pmvnorm}} in package \pkg{mvtnorm} to perform the necessary integrations. The cumulative hazard detection function \code{detectfn = 'HHN'} is preferred because it is interpretable as the distribution of activity.

}

\value{

A numeric vector of length equal to the number of rows in \code{object$mask}.

}

\seealso{

\code{\link{PG}}

}

\examples{

fitCL <- secr.fit(captdata, model = lambda0~1, detectfn = 'HHN', CL = TRUE)
pPG <- predictPG(fitCL, plt = TRUE, dots = FALSE, buffer = 30)
plot(traps(captdata), add = TRUE)

}

\keyword{ model }
