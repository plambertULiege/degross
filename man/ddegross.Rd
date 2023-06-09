% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddegross.R
\name{ddegross}
\alias{ddegross}
\title{Density function based on an object resulting from the estimation procedure in \link{degross}.}
\usage{
ddegross(x, degross.fit, phi)
}
\arguments{
\item{x}{Scalar or vector where the fitted density must be evaluated.}

\item{degross.fit}{A \link{degross.object} generated using \link{degross} and containing the density estimation results.}

\item{phi}{(Optional) vector of spline parameters for the log density (default: \code{degross.fit$phi} if missing).}
}
\value{
A scalar or vector of the same length as \code{x} containing the value of the fitted density at \code{x}.
}
\description{
Density function based on an object resulting from the estimation procedure in \link{degross}.
}
\examples{
## Generate grouped data
sim = simDegrossData(n=1500, plotting=TRUE, choice=2)

## Create a degrossData object
obj.data = degrossData(Big.bins=sim$Big.bins, freq.j=sim$freq.j, m.j=sim$m.j)
print(obj.data)

## Estimate the density
obj.fit = degross(obj.data)

## Superpose the fitted density using the <ddegross> function
curve(ddegross(x,obj.fit),add=TRUE,lty="dashed")
legend("topright",lty="dashed",lwd=2,legend="Estimated",box.lty=0, inset=.04)

}
\references{
Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
}
\seealso{
\code{\link{degross.object}}, \code{\link{pdegross}}, \code{\link{qdegross}}, \code{\link{degross}}.
}
\author{
Philippe Lambert \email{p.lambert@uliege.be}
}
