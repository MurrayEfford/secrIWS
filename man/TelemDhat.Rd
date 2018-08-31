\name{TelemDhat}
\alias{TelemDhat}
\alias{estimatePG}

\title{Telemetry-based Density Estimate}

\description{

This function simulates the telemetry-based density estimator of Ivan et
al. 2013). It is used internally by \code{\link{runsim}} and will not
usually be called by the user.

The function \code{estimatePG} simulates telemetry sampling to determine the proportion of each detected animal's activity that falls in the inner zone. \code{estimatedPG} is called by \code{TelemDhat}.

}

\usage{
estimatePG(CH, NTelem, pTelem)

TelemDhat(CH, PG, innerarea)
}

\arguments{
  \item{CH}{capthist object as defined in secr }
  \item{NTelem}{integer number of telemetry fixes per animal}
  \item{pTelem}{numeric proportion of animals for which telemetry data obtained}
  \item{PG}{numeric vector of length n output by \code{estimatePG}}
  \item{innerarea}{numeric area of trapping grid for telemetry determinations}
}

\details{

For \code{estimatePG} the input capthist object must have an individual covariate `truep'
giving for each detected animal the true proportion of activity on the trapping grid 
i.e. inner area.

The individuals for which telemetry data are simulated in \code{estimatePG} are a binomial
random sample of the rows in CH. The estimator fails when this hits zero by chance.

}

\value{

For \code{estimatePG}, a vector of PG (overlap) values, one per detected animal.

For \code{TelemDhat}, scalar estimate of density (animals per hectare), or NA if density cannot be estimated.

}

\references{

Ivan, J. S., White, G. C. and Shenk, T. M. (2013a) Using auxiliary telemetry information to estimate animal density from capture--recapture data. \emph{Ecology} \bold{94}, 809--816.

}

\seealso{

\code{\link{runsim}}

}

\keyword{ datagen }

