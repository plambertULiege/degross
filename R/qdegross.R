#' Quantile function based on an object resulting from the estimation procedure in \link{degross}.
#'
#' @usage qdegross(p, degross.fit, phi, get.se=FALSE, cred.level=.95, eps=1e-4)
#' @importFrom stats integrate qnorm
#' @param p Scalar or vector of probabilities in (0,1) indicating the requested fitted quantiles Q(p) based on the density estimation results in \code{degross.fit}.
#' @param degross.fit A \code{\link{degross.object}} generated using \link{degross} and containing the density estimation results.
#' @param phi (Optional) vector of spline parameters for the log density (default: \code{degross.fit$phi} if missing).
#' @param get.se Logical indicating if standard errors for Q(p) are requested (default: FALSE).
#' @param cred.level Level of credible intervals for Q(p).
#' @param eps Precision with which each quantile should be computed (default: 1e-4).
#'
#' @return A scalar or vector \code{x} of the same length as \code{p} containing the values Q(p) at which the cdf \code{pdegross(x,degross.fit)} is equal to \code{p}.
#'         When \code{get.se} is TRUE, a vector or a matrix containing the quantile estimate(s), standard errors and credible interval limits for Q(p) is provided.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#' @seealso \code{\link{degross.object}}, \code{\link{ddegross}}, \code{\link{pdegross}}, \code{\link{degross}}.
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
#' ## Corresponding fitted quantiles
#' p = c(.01,.05,seq(.1,.9,by=.1),.95,.99) ## Desired probabilities
#' Q.p = qdegross(p,obj.fit) ## Compute the desired quantiles
#' print(Q.p) ## Estimated quantiles
#'
#' ## Compute the standard error and a 90% credible interval for the 60% quantile
#' Q.60 = qdegross(.60,obj.fit,get.se=TRUE,cred.level=.90) ## Compute the desired quantile
#' print(Q.60) ## Estimated quantile, standard error and credible interval
#'
#' @export
qdegross = function(p,degross.fit,phi,get.se=FALSE,cred.level=.95,eps=1e-4){
    obj = degross.fit
    x = 0*p
    x[(p<0)|(p>1)] = NA
    ymin = min(obj$degross.data$Big.bins) ; ymax = max(obj$degross.data$Big.bins) ## Support of the density
    idx0 = which(p==0) ; if (length(idx0)>0) x[idx0] = ymin
    idx1 = which(p==1) ; if (length(idx1)>0) x[idx1] = ymax
    idx = which((p>0)&(p<1))
    if (length(idx) > 0){
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
            phi = obj$phi
            pi.i = obj$pi.i
        }
        ##
        F.i = cumsum(pi.i)
        ui = with(obj, ui)
        gg = function(p) ui[which.min(abs(F.i-p))]
        x[idx] = sapply(p[idx],gg) ## Good starting guesses for the desired quantiles
        ok = FALSE
        while(!ok){
            dx = (p[idx]-pdegross(x[idx],obj,phi)) / ddegross(x[idx],obj,phi)
            x[idx] = x[idx] + dx
            ok = all(abs(dx)<eps)
        }
    }
    names(x) = p
    ##
    if (get.se){ ## Computation of s.e. and credible interval for Q(p)
        Q.p = x
        ## (a) Partial derivative of the estimated quantile function wrt <phi.k>
        ##      Careful: p & Q.p should be scalars (i.e. length 1)
        D1Q = function(p,Q.p,k){
            knots = obj$degross.data$knots
            bf = function(x,k) cubicBsplines::Bsplines(x,knots)[,k] * ddegross(x,obj,phi)
            ## Q.p = qdegross(p,obj,phi)
            const = (1/ddegross(Q.p,obj,phi))
            temp = integrate(bf,min(knots),Q.p,k=k)$val
            temp = temp - p*integrate(bf,min(knots),max(knots),k=k)$val
            ##
            ans = const*temp
            return(ans)
        }
        ## (b) Gradient of the estimated quantile function wrt <phi>
        ##      Careful: p & Q.p should be scalars (i.e. length 1)
        gradQ = function(p,Q.p){
            K = obj$degross.data$K
            ans = rep(0,K)
            for (k in 1:K) ans[k] = D1Q(p,Q.p,k)
            return(ans)
        }
        ##
        se.Q = 0*x
        idx = which.max(obj$phi)
        V.mj = with(obj, solve(Fisher.mj[-idx,-idx]))
        ## VQ.ni = rep(0,length(p))
        VQ.mj = rep(0,length(p))
        for (k in 1:length(p)){ ## Loop over quantiles
            gr = gradQ(p[k],Q.p[k])[-idx] ## \partial Q(p[k]|theta) \partial theta
            VQ.mj[k] = sum(gr * c(V.mj%*%gr))   ## Variance of the estimated quantiles
        }
        se.Q = sqrt(VQ.mj)
        names(se.Q) = names(x)
        ##
        z.alpha = qnorm(1-.5*(1-cred.level))
        ci.low = Q.p - z.alpha * se.Q
        ci.up  = Q.p + z.alpha * se.Q
        cred.int = rbind(c(x),se.Q,ci.low,ci.up)
        rownames(cred.int) = c("Estimate","s.e.","ci.low","ci.up")
        x = t(cred.int)
        attr(x,"cred.level") = cred.level
    }
    return(x)
}
