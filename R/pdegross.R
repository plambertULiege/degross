#' Cumulative distribution function (cdf) based on an object resulting from the estimation procedure in \link{degross}.
#'
#' @usage pdegross(x, degross.fit, phi)
#' @param x Scalar or vector where the fitted cdf must be evaluated.
#' @param degross.fit A \code{\link{degross.object}} generated using \link{degross} and containing the density estimation results.
#' @param phi (Optional) vector of spline parameters for the log density (default: \code{degross.fit$phi} if missing).
#'
#' @return a scalar or vector of the same length as \code{x} containing the value of the fitted cdf at \code{x}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#' @seealso \code{\link{degross.object}}, \code{\link{ddegross}}, \code{\link{qdegross}}, \code{\link{degross}}.
#'
#' @examples
#' ## Generate grouped data
#' sim = simDegrossData(n=3500, plotting=TRUE, choice=2)
#'
#' ## Create a degrossData object
#' obj.data = degrossData(Big.bins=sim$Big.bins, freq.j=sim$freq.j, m.j=sim$m.j)
#' print(obj.data)
#'
#' ## Estimate the density
#' obj.fit = degross(obj.data)
#'
#' ## Superpose the fitted cdf using the <pdegross> function
#' with(sim, curve(true.cdf(x),min(Big.bins),max(Big.bins),
#'      col="red",lwd=2, ylab="F(x)"))
#' curve(pdegross(x,obj.fit),add=TRUE,lty="dashed")
#' legend("topleft", legend=c("Target cdf","Estimated cdf"), lwd=2,
#'        lty=c("solid","dashed"), col=c("red","black"), box.lty=0, inset=.04)
#'
#' @export
pdegross = function(x,degross.fit,phi){
    F.x = 0*x
    obj = degross.fit
    ymin = min(obj$degross.data$Big.bins) ; ymax = max(obj$degross.data$Big.bins) ## Support of the density
    idx = which((x >= ymin) & (x <= ymax))
    idx0 = which(x < ymin) ; idx1 = which(x > ymax)
    if (length(idx0 > 0)) F.x[idx0] = 0.0
    if (length(idx1 > 0)) F.x[idx1] = 1.0
    if (length(idx) > 0){
        xx = c(ymin,obj$ui + .5*obj$delta)
        if (!missing(phi)){  ## Recompute fitted small bin probabilities if new <phi> vector provided
            K = obj$degross.data$K ## Number of B-splines in the basis
            if (length(phi)!=K){
                cat("<phi> should be of length ",K,"\n")
                return(NULL)
            }
            eta.i = c(obj$degross.data$B.i %*% phi)
            temp = exp(eta.i)
            pi.i = temp/sum(temp)
        } else {  ## ... otherwise, use the estimated values for <pi.i>
            pi.i = obj$pi.i
        }
        ##
        yy = c(0,cumsum(pi.i))
        F.x[idx] = stats::splinefun(xx,yy)(x[idx])
    }
    return(F.x)
}
