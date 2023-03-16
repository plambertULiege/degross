#' Log-posterior for given spline parameters, big bin (and optional: small bin) frequencies, tabulated sample moments and roughness penalty parameter.
#' Compared to \link{degross_lpost}, no Fisher information matrix is computed and the gradient evaluation is optional, with a resulting computational gain.
#' @usage degross_lpostBasic(phi, tau, n.i, degross.data,
#'                           use.moments = rep(TRUE,4), freq.min = 20, diag.only=FALSE,
#'                           gradient=FALSE,
#'                           penalize = TRUE, aa = 2, bb = 1e-6, pen.order = 3)
#' @param phi Vector of K B-spline parameters \eqn{\phi} to specify the log-density.
#' @param tau Roughness penalty parameter.
#' @param n.i Small bin frequencies.
#' @param degross.data A \link{degrossData.object} created using the \link{degrossData} function.
#' @param use.moments Vector with 4 logicals indicating which tabulated sample moments to use as soft constraints. Defaults: \code{rep(TRUE,4)}.
#' @param freq.min Minimal big bin frequency required to use the corresponding observed moments as soft constraints. Default: \code{20}.
#' @param diag.only Logical indicating whether to ignore the off-diagonal elements of the variance-covariance matrix of the sample central moments. Default: FALSE.
#' @param gradient Logical indicating if the gradient (Score) of the \eqn{\log p(\phi|\tau,data)} should be computed (default: FALSE).
#' @param penalize Logical indicating whether a roughness penalty of order \code{pen.order} is required (with \eqn{tau \sim G(aa,bb)}). Default: \code{TRUE}.
#' @param aa Real giving the first parameter in the Gamma prior for \code{tau}. Default: \code{2}.
#' @param bb Real giving the second parameter in the Gamma prior for \code{tau}. Default: \code{1e-6}.
#' @param pen.order Integer giving the order of the roughness penalty. Default: \code{3}.
#'
#' @return A list containing :
#' \itemize{
#' \item{\code{lpost.ni} : \verb{ }}{value of the log-posterior based on the given small bin frequencies \code{n.i} and the tabulated sample moments.}
#' \item{\code{lpost.mj} : \verb{ }}{value of the log-posterior based on the big bin frequencies \code{degross.data$freq.j} and the tabulated sample moments.}
#' \item{\code{llik.ni} : \verb{ }}{multinomial log-likelihood based on the given small bin frequencies \code{n.i}.}
#' \item{\code{llik.mj} : \verb{ }}{multinomial log-likelihood based on the big bin frequencies \code{degross.data$freq.j} resulting from \code{n.i}.}
#' \item{\code{moments.penalty} : \verb{ }}{log of the joint (asymptotic) density for the observed sample moments.}
#' \item{\code{penalty} : \verb{ }}{\eqn{\log p(\phi|\tau) + \log p(\tau)}.}
#' \item{\code{M.j} : \verb{ }}{theoretical moments of the density (resulting from \eqn{\phi}) within a big bin.}
#' \item{\code{pi.i} : \verb{ }}{small bin probabilities.}
#' \item{\code{ui} : \verb{ }}{small bin midpoints.}
#' \item{\code{delta} : \verb{ }}{width of the small bins.}
#' \item{\code{gamma.j} : \verb{ }}{big bin probabilities.}
#' \item{\code{tau} : \verb{ }}{reminder of the value of the roughness penalty parameter \eqn{\tau}.}
#' \item{\code{phi} : \verb{ }}{reminder of the vector of spline parameters (defining the density).}
#' \item{\code{n.i} : \verb{ }}{reminder of the small bin frequencies given as input.}
#' \item{\code{freq.j} : \verb{ }}{reminder of the big bin frequencies in \code{degross.data$freq.j}.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#' @seealso \code{\link{degross_lpost}}, \code{\link{degross}}, \code{\link{degross.object}}.
#'
#' @examples
#' sim = simDegrossData(n=3500, plotting=TRUE,choice=2) ## Generate grouped data
#' obj.data = degrossData(Big.bins=sim$Big.bins, freq.j=sim$freq.j, m.j=sim$m.j)
#' print(obj.data)
#' obj.fit = degross(obj.data) ## Estimate the underlying density
#' phi.hat = obj.fit$phi ; tau.hat = obj.fit$tau
#' ## Evaluate the log-posterior at convergence
#' res = degross_lpostBasic(phi=phi.hat, tau=tau.hat, degross.data=obj.data,
#'                          gradient=TRUE)
#' print(res)
#'
#' @export
degross_lpostBasic = function(phi,tau,n.i,degross.data,
                              use.moments=rep(TRUE,4), freq.min=20, diag.only=FALSE,
                              gradient=FALSE,
                              penalize=TRUE,aa=2,bb=1e-6,pen.order=3){
    if (missing(n.i)) n.i = NULL
    ## Spline parameters
    K = length(phi)
    ##
    ## Penalty matrix
    Dd = diff(diag(K),diff=pen.order)
    Pd = t(Dd) %*% Dd
    ##
    ## freq.min = 20 ## Minimal big bin frequency for using the observed moments as soft constraints
    ## Data
    freq.j = degross.data$freq.j # Frequencies for big bins
    m.tot = sum(freq.j)
    J = length(freq.j) # Nbr of big bins
    m.j = degross.data$m.j # Sample central moments
    small.to.big = degross.data$small.to.big # Big bin to which each small bin belongs
    ui = degross.data$ui ## Midpoints of the small bins
    delta = diff(ui[1:2]) ## Width of a small bin
    ##
    ## Fitted mean for small & big bins + small bin probs
    B.i = degross.data$B.i
    eta.i = c(B.i %*% phi)
    temp = exp(eta.i) ## Small bin means
    pi.i = temp / sum(temp) ## Small bin probs
    gamma.j = c(rowsum(pi.i,small.to.big)) # Big bin probs
    ## gamma.j = c(tapply(pi.i,small.to.big,sum)) # Big bin probs
    ##
    ## Theoretical moments of the fitted density within a big bin
    M.j = matrix(nrow=J,ncol=8)
    rownames(M.j) = paste("Bin",1:J,sep="")
    colnames(M.j) = paste("mu",1:8,sep="")
    M.j[,1] = rowsum(ui*pi.i,small.to.big) / gamma.j # Theoretical mean within big bin j
    ## M.j[,1] = tapply(ui*pi.i,small.to.big,sum) / gamma.j # Theoretical mean within big bin j
    for (r in 2:8){
        ## kth theoretical central moment within each of the big bins
        M.j[,r] = rowsum((ui-M.j[small.to.big,1])^r*pi.i,small.to.big) / gamma.j
        ## M.j[,r] = tapply((ui-M.j[small.to.big,1])^r*pi.i,small.to.big,sum) / gamma.j
    }
    ## Gradient
    if (gradient){
        Bi.tilde = t(t(B.i) - c(apply(pi.i*B.i,2,sum)))
        PiBi.tilde = pi.i*Bi.tilde
        ## Score.ni
        if (!is.null(n.i)) Score.ni = c(t(B.i)%*%(n.i-m.tot*pi.i))
        else Score.ni = NULL
        ## Score.mj
        Score.mj = Score.lprior = rep(0,ncol(B.i))
        for (j in 1:J){
            idx = which(small.to.big == j)
            Score.mj = Score.mj + freq.j[j]/gamma.j[j]*c(apply(PiBi.tilde[idx,],2,sum))
        }
    }
    R = 4 ## Up to moments of order 4 of interest
    dmujr.k = array(dim=c(J,R,K)) #  \partial\mu_{rj} / \partial\theta_k
    ##
    ## Contribution due to sample central moments within each big bin
    if (sum(use.moments) > 0){
        invSigma.j = array(dim=c(J,4,4))
        for (j in 1:J){
            idx = which(small.to.big == j)
            if (gradient){
                ## \partial\mu_{1j} / \partial\theta_k
                dmujr.k[j,1,] = c(apply((ui[idx]-M.j[j,1])*pi.i[idx]*Bi.tilde[idx,],2,sum)/gamma.j[j])
                ## \partial\mu_{rj} / \partial\theta_k  (r > 1)
                for (r in 2:4){
                    dmujr.k[j,r,] = c(apply((ui[idx]-M.j[j,1])^r*pi.i[idx]*Bi.tilde[idx,],2,sum)) / gamma.j[j]
                    temp = sum((ui[idx]-M.j[j,1])^(r-1)*pi.i[idx])
                    dmujr.k[j,r,] = dmujr.k[j,r,] - r * temp * dmujr.k[j,1,] / gamma.j[j]
                    temp2 = c(apply(pi.i[idx]*Bi.tilde[idx,],2,sum)/gamma.j[j])
                    dmujr.k[j,r,] = dmujr.k[j,r,] - M.j[j,r] * temp2
                }
            }
            ##
            ## Sample moments of interest with non-NA values AND big bin sample size >= freq.min
            idx2 = which((freq.j[j] >= freq.min) & (!is.na(m.j[j,])) & use.moments)
            if (length(idx2) > 0){   # i.e. if at least one observed sample moment of interest
                Mat = Sigma_fun(M.j[j,1:8]) # Theoretical covariance matrix for moments
                if (diag.only) Mat = diag(diag(Mat)) # Just keep variances (as covariances rely on moments of very high orders !)
                invSigma.j[j,,] = freq.j[j] * solve(Mat) # Theoretical version corrected for frequencies
                if (gradient){
                    Mat = matrix(invSigma.j[j,idx2,idx2],nrow=length(idx2)) # Make sure that sub-matrix is a matrix...
                    temp2 = matrix(dmujr.k[j,idx2,],nrow=length(idx2))
                    temp = Mat %*% temp2
                    Score.lprior = Score.lprior + matrix(m.j[j,idx2]-M.j[j,idx2],nrow=1) %*% temp
                }
            }
        }
    }
    ## Contribution due to roughness penalty
    if (gradient){
        Score.lprior = c(Score.lprior)
        if (penalize){
            Score.lprior = Score.lprior - tau * c(Pd%*%phi)
        }
        ##
        if (!is.null(n.i)) Score.ni = Score.ni + Score.lprior
        Score.mj = Score.mj + Score.lprior
    }
    ## lpost
    llik.ni = NULL
    if (!is.null(n.i)) llik.ni = sum(n.i*log(pi.i)) ## log complete likelihood
    llik.mj = sum(freq.j * log(gamma.j)) ## log observed likelihood
    ##
    moments.penalty = 0
    if (sum(use.moments) > 0){
        for (j in 1:J){
            ## Sample moments of interest with non-NA values AND big bin sample size >= 20
            idx2 = which((freq.j[j] >= freq.min) & (!is.na(m.j[j,])) & use.moments)
            if (sum(idx2) > 0){
                temp = matrix(invSigma.j[j,idx2,idx2],nrow=length(idx2))
                eig.vals = svd(temp)$d ; eig.vals = eig.vals[eig.vals > 1e-4]
                moments.penalty = moments.penalty + .5*sum(log(eig.vals))
                moments.penalty = moments.penalty -.5 * sum(c(m.j[j,idx2]-M.j[j,idx2]) * c(temp %*% c(m.j[j,idx2]-M.j[j,idx2])))
            }
        }
    }
    ## Roughness penalty
    penalty = 0
    if (penalize){
        penalty = (aa+.5*nrow(Dd))*log(tau) -tau*(bb+.5*sum(phi*c(Pd%*%phi)))
    }
    ##
    lpost.ni = llik.ni + moments.penalty + penalty    ## log complete posterior
    lpost.mj = llik.mj + moments.penalty + penalty ## log observed posterior
    ##
    ## Output
    res = list()
    res$lpost.ni = lpost.ni
    res$lpost.mj = lpost.mj
    res$llik.ni = llik.ni
    res$llik.mj = llik.mj
    res$moments.penalty = moments.penalty
    res$penalty = penalty
    ## res$Score = Score
    if (gradient){
        if (!is.null(n.i)) res$Score.ni = Score.ni
        res$Score.mj = Score.mj
        res$Score.lprior = Score.lprior
    }
    res$M.j = M.j
    res$pi.i = pi.i
    res$ui = ui ; res$delta = delta
    res$gamma.j = gamma.j
    res$tau = tau
    res$phi = phi
    res$n.i = n.i
    res$freq.j = freq.j
    ##
    return(res)
}
