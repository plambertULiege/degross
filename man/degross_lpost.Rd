% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/degross_lpost.R
\name{degross_lpost}
\alias{degross_lpost}
\title{Log-posterior (with gradient and Fisher information) for given spline parameters, small bin frequencies, tabulated sample moments and roughness penalty parameter.
This function is maximized during the M-step of the EM algorithm to estimate the B-spline parameters entering the density specification.}
\usage{
degross_lpost(phi, tau, n.i, degross.data,
                     use.moments = rep(TRUE,4), freq.min = 20, diag.only=FALSE,
                     penalize = TRUE, aa = 2, bb = 1e-6, pen.order = 3)
}
\arguments{
\item{phi}{Vector of K B-spline parameters \eqn{\phi} to specify the log-density.}

\item{tau}{Roughness penalty parameter.}

\item{n.i}{Small bin frequencies.}

\item{degross.data}{A \link{degrossData.object} created using the \link{degrossData} function.}

\item{use.moments}{Vector with 4 logicals indicating which tabulated sample moments to use as soft constraints. Defaults: \code{rep(TRUE,4)}.}

\item{freq.min}{Minimal big bin frequency required to use the corresponding observed moments as soft constraints. Default: \code{20}.}

\item{diag.only}{Logical indicating whether to ignore the off-diagonal elements of the variance-covariance matrix of the sample central moments. Default: FALSE.}

\item{penalize}{Logical indicating whether a roughness penalty of order \code{pen.order} is required (with \eqn{\tau \sim G(aa,bb)}). Default: \code{TRUE}.}

\item{aa}{Positive real giving the first parameter in the Gamma prior for \code{tau}. Default: \code{2}.}

\item{bb}{Positive real giving the second parameter in the Gamma prior for \code{tau}. Default: \code{1e-6}.}

\item{pen.order}{Integer giving the order of the roughness penalty. Default: \code{3}.}
}
\value{
A list containing :
\itemize{
\item{\code{lpost}, \code{lpost.ni} : \verb{ }}{value of the log-posterior based on the given small bin frequencies \code{n.i} and the tabulated sample moments.}
\item{\code{lpost.mj} : \verb{ }}{value of the log-posterior based on the big bin frequencies \code{degross.data$freq.j} and the tabulated sample moments.}
\item{\code{llik.ni} : \verb{ }}{multinomial log-likelihood based on the given small bin frequencies \code{n.i}.}
\item{\code{llik.mj} : \verb{ }}{multinomial log-likelihood based on the big bin frequencies \code{degross.data$freq.j}.}
\item{\code{moments.penalty} : \verb{ }}{log of the joint (asymptotic) density for the observed sample moments.}
\item{\code{penalty} : \verb{ }}{\eqn{\log p(\phi|\tau) + \log p(\tau)}.}
\item{\code{Score}, \code{Score.ni} : \verb{ }}{score (w.r.t. \eqn{\phi}) of \code{lpost.ni}.}
\item{\code{Score.mj} : \verb{ }}{score (w.r.t. \eqn{\phi}) of \code{lpost.mj}.}
\item{\code{Fisher} & \code{Fisher.ni}: \verb{ }}{information matrix (w.r.t. \eqn{\phi}) of \code{lpost.ni}.}
\item{\code{Fisher.mj} : \verb{ }}{information matrix (w.r.t. \eqn{\phi}) of \code{lpost.mj}.}
\item{\code{M.j} : \verb{ }}{theoretical moments of the density (resulting from \eqn{\phi}) within a big bin.}
\item{\code{pi.i} : \verb{ }}{small bin probabilities.}
\item{\code{ui} : \verb{ }}{small bin midpoints.}
\item{\code{delta} : \verb{ }}{width of the small bins.}
\item{\code{gamma.j} : \verb{ }}{Big bin probabilities.}
\item{\code{tau} : \verb{ }}{reminder of the value of the roughness penalty parameter \eqn{\tau}.}
\item{\code{phi} : \verb{ }}{reminder of the vector of spline parameters (defining the density).}
\item{\code{n.i} : \verb{ }}{reminder of the small bin frequencies given as input.}
}
}
\description{
Log-posterior (with gradient and Fisher information) for given spline parameters, small bin frequencies, tabulated sample moments and roughness penalty parameter.
This function is maximized during the M-step of the EM algorithm to estimate the B-spline parameters entering the density specification.
}
\examples{
sim = simDegrossData(n=3500, plotting=TRUE,choice=2) ## Generate grouped data
obj.data = degrossData(Big.bins=sim$Big.bins, freq.j=sim$freq.j, m.j=sim$m.j)
print(obj.data)
obj.fit = degross(obj.data) ## Estimate the underlying density
## Evaluate the log-posterior at convergence
res = with(obj.fit, degross_lpost(phi, tau, n.i, obj.data, diag.only=diag.only))
print(res$Score) ## Score of the log posterior at convergence

}
\references{
Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
}
\seealso{
\code{\link{degross_lpostBasic}}, \code{\link{degross}}, \code{\link{degross.object}}.
}
\author{
Philippe Lambert \email{p.lambert@uliege.be}
}
