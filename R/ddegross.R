#' Density function based on an object resulting from the estimation procedure in \link{degross}.
#'
#' @usage ddegross(x, degross.fit, phi)
#'
#' @param x Scalar or vector where the fitted density must be evaluated.
#' @param degross.fit A \link{degross.object} generated using \link{degross} and containing the density estimation results.
#' @param phi (Optional) vector of spline parameters for the log density (default: \code{degross.fit$phi} if missing).
#'
#' @return A scalar or vector of the same length as \code{x} containing the value of the fitted density at \code{x}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#'
#' @seealso \code{\link{degross.object}}, \code{\link{pdegross}}, \code{\link{qdegross}}, \code{\link{degross}}.
#'
#' @examples
#' ## Generate grouped data
#' sim = simDegrossData(n=1500, plotting=TRUE, choice=2)
#'
#' ## Create a degrossData object
#' obj.data = degrossData(Big.bins=sim$Big.bins, freq.j=sim$freq.j, m.j=sim$m.j)
#' print(obj.data)
#'
#' ## Estimate the density
#' obj.fit = degross(obj.data)
#'
#' ## Superpose the fitted density using the <ddegross> function
#' curve(ddegross(x,obj.fit),add=TRUE,lty="dashed")
#' legend("topright",lty="dashed",lwd=2,legend="Estimated",box.lty=0, inset=.04)
#'
#' @export
ddegross = function(x,degross.fit,phi){
  f.x = 0*x
  obj = degross.fit
  ymin = min(obj$degross.data$Big.bins) ; ymax = max(obj$degross.data$Big.bins) ## Support of the density
  idx = which((x >= ymin) & (x <= ymax))
  if (length(idx) >0){
    y = x[idx]
    K = obj$degross.data$K ## Number of B-splines in the basis
    if (!missing(phi)){  ## Recompute fitted density if new <phi> vector provided
        if (length(phi)!=K){
            cat("<phi> should be of length ",K,"\n")
            return(NULL)
        }
        eta.i = c(obj$degross.data$B.i %*% phi)
        temp = exp(eta.i)
        logNormCst = -log(sum(temp)) - log(obj$degross.data$delta)
    } else { ## ... otherwise, use the estimated value for <phi> in <degross.fit>
        phi = obj$phi ## Estimated spline parameters
        logNormCst = obj$logNormCst
    }
    ##
    B.y = cubicBsplines::Bsplines(y,seq(ymin,ymax,length=K-2)) ## B-spline basis at <y>
    eta.y = c(B.y %*% phi)
    f.y = exp(eta.y + logNormCst)
    f.x[idx] = f.y
  }
  return(f.x)
}
