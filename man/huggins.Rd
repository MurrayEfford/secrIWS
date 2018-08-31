\name{hugginsN}
\alias{hugginsN}
\alias{addDistanceToEdge}

\title{Huggins Closed Population Estimate}

\description{
  
A partial implementation of Huggins' conditional likelihood method.

}

\usage{

hugginsN(CH)
addDistanceToEdge(CH)

}

\arguments{
  \item{CH}{secr capthist object}
}

\details{

  If the capthist object has a covariates attribute with the column
  'distancetoedge' this will be used as an individual covariate.

  \code{addDistanceToEdge} adds the required covariate.

  [A more complete interface will follow].
  
}

\value{
  
  For hugginsN, a naked scalar estimate of N.

  For addDistanceToEdge, a capthist object with covariates attribute
  including the column `distancetoedge'.
  
}

\references{

Huggins, R. M. (1989) On the statistical analysis of capture
experiments. \emph{Biometrika} \bold{76}, 133--140.

Huggins, R. M. (1991) Some practical aspects of a conditional likelihood
  approach to capture experiments. \emph{Biometrics} \bold{47}, 725--732.

}

\examples{

hugginsN(captdata)
captdata2 <- addDistanceToEdge(captdata)
hugginsN(captdata2)

}

\keyword{ model }
