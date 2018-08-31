\name{plotmaskedge}
\alias{plotmaskedge}

\title{Outline Around Mask Cells}

\description{
  
Plots the outer edge of a secr mask.

}
\usage{

plotmaskedge(mask, add = TRUE, ...)

}

\arguments{
  \item{mask}{secr habitat mask object}
  \item{add}{logical; if TRUE the line is added to an existing plot}
  \item{\dots}{other line plotting arguments passed to \code{\link{segments}}}
}

\details{
  
  May be slow. 
  
}

\value{
  None
}

\examples{

HR <- build.irregular()
plotmaskedge (HR, add = FALSE)

}

\keyword{ hplot }
