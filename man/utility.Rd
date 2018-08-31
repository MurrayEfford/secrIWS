\name{utility}
\alias{rCH}
\alias{Px}
\alias{HRcentres}

\title{Miscellaneous Functions}

\description{Functions used internally and not usually called by the user.}

\usage{

Px(trapsites, template)
rCH(px, p1, noccasions, hazard = FALSE)
HRcentres(pop)

}

\arguments{
  \item{trapsites}{secr traps object}
  \item{template}{secr mask object comprising the cell centres of one home range}

  \item{px}{numeric vector of capture probabilities as returned by \code{Px}}
  \item{p1}{numeric overall detection probability}
  \item{noccasions}{integer number of sampling occasions}
  \item{hazard}{logical; if TRUE the model is expressed in terms of
    hazard rather than probability (NOT YET USED) }
  \item{pop}{secr popn object to summarise}
}
\details{

  \code{Px} returns the probability of capture at each trap.

  \code{rCH} generates one capture history.

  \code{HRcentres} computes home-range centres (centroids) using the
  weightings in the covariate `p'.

}

\value{
  
\code{Px} returns a vector of trap-specific capture probabilities.
  
\code{rCH} returns a capture history, a vector of integer values either
0 (no capture) or the trap in which this animal was caught, on each occasion.

\code{HRcentres} returns a 2-column matrix, the x-y coordinates of the
centroids (one row per animal).

}

\keyword{ datagen }
