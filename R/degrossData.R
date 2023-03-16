#' Creates a \link{degrossData.object} from the observed tabulated frequencies and central moments.
#' @usage degrossData(Big.bins, freq.j, m.j, I=300, K=25)
#' @importFrom cubicBsplines Bsplines
#' @param Big.bins Vector of length \code{J+1} with the limits of the \code{J} big bins containing the data used to produce the tabulated statistics.
#' @param freq.j The number of data observed within each big bin.
#' @param m.j A matrix of dim \code{J} by 4 giving the first 4 sample central moments within each of the \code{J} big bins.
#' @param I The number of small bins used for quadrature during the normalization of the density during its estimation. Default: \code{300}.
#' @param K The desired number of B-splines in the basis used for density estimation. Default= \code{25}.
#'
#' @return A \link{degrossData.object}, i.e. a list containing:
#' \itemize{
#'   \item{\code{small.bins} : \verb{ }}{a vector of length \code{I+1} with the small bin limits.}
#'   \item{\code{ui} : \verb{ }}{the \code{I} midpoints of the small bins.}
#'   \item{\code{delta} : \verb{ }}{width of the small bins.}
#'   \item{\code{I} : \verb{ }}{the number of small bins.}
#'   \item{\code{B.i} : \verb{ }}{a matrix of dim \code{I} by \code{K} with the B-spline basis evaluated at the small bin midpoints.}
#'   \item{\code{K} : \verb{ }}{number of B-splines in the basis.}
#'   \item{\code{knots} : \verb{ }}{equidistant knots supporting the B-splines basis.}
#'   \item{\code{Big.bins} : \verb{ }}{vector of length \code{J+1} with the limits of the \code{J} big bins containing the data used to produce the tabulated statistics.}
#'   \item{\code{freq.j} : \verb{ }}{the number of data observed within each big bin.}
#'   \item{\code{m.j} : \verb{ }}{a matrix of dim \code{J} by 4 giving the first 4 sample central moments within each big bin.}
#'   \item{\code{J} : \verb{ }}{the number of big bins.}
#'   \item{\code{small.to.big} : \verb{ }}{a vector of length \code{I} indicating to what big bin each element of \code{ui} belongs.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#' @seealso \code{\link{print.degrossData}}, \code{\link{degross}}.

#' @examples
#' sim = simDegrossData(n=3500, plotting=TRUE)
#' obj.data = degrossData(Big.bins=sim$Big.bins, freq.j=sim$freq.j, m.j=sim$m.j)
#' print(obj.data)
#'
#' @export
degrossData = function(Big.bins, freq.j, m.j,
                   I=300, K=25){
  ymin = min(Big.bins) ; ymax = max(Big.bins)
  small.bins = seq(ymin,ymax,length=I+1) ## Small bins
  delta = diff(small.bins[1:2])  ## Width of small bins
  ui = small.bins[-1] - .5*delta ## Midpoints of the small bins
  small.to.big = as.numeric(cut(ui,Big.bins))  ## Id of the big bin to which a small bin belongs
  knots = seq(ymin,ymax,length=K-2) ## knots
  B.i = Bsplines(ui,knots) ## B-spline matrix at midpoints of small bins
  ##
  degross.data = list(small.bins=small.bins, ui=ui, delta=delta, I=length(ui),
                      B.i=B.i, K=K, knots=knots,
                      Big.bins=Big.bins, freq.j=freq.j, m.j=m.j, J=length(freq.j),
                      small.to.big=small.to.big
  )
  ##
  return(structure(degross.data,class="degrossData"))
}
