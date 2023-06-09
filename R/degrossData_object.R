#' Object generated from grouped summary statistics, including tabulated frequencies and central moments of order 1 up to 4, to estimate the underlying density using \code{\link{degross}}.
#'
#' An object returned by the \code{\link{degrossData}} function from tabulated frequencies and central moments of order 1 up to 4. It is used in a second step by \code{\link{degross}} to estimate the underlying density.
#'
#' @return A list containing :
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
#' @seealso \code{\link{degrossData}}, \code{\link{print.degrossData}}
#'
#' @name degrossData.object
NULL

