#' Log-posterior (with gradient and Fisher information) for given spline parameters, small bin frequencies, tabulated sample moments and roughness penalty parameter.
#' This function is maximized during the M-step of the EM algorithm to estimate the B-spline parameters entering the density specification.
#' @usage degross_lpost(phi, tau, n.i, degross.data,
#'                      use.moments = rep(TRUE,4), freq.min = 20, diag.only=FALSE,
#'                      penalize = TRUE, aa = 2, bb = 1e-6, pen.order = 3)
#' @param phi Vector of K B-spline parameters \eqn{\phi} to specify the log-density.
#' @param tau Roughness penalty parameter.
#' @param n.i Small bin frequencies.
#' @param degross.data A \link{degrossData.object} created using the \link{degrossData} function.
#' @param use.moments Vector with 4 logicals indicating which tabulated sample moments to use as soft constraints. Defaults: \code{rep(TRUE,4)}.
#' @param freq.min Minimal big bin frequency required to use the corresponding observed moments as soft constraints. Default: \code{20}.
#' @param diag.only Logical indicating whether to ignore the off-diagonal elements of the variance-covariance matrix of the sample central moments. Default: FALSE.
#' @param penalize Logical indicating whether a roughness penalty of order \code{pen.order} is required (with \eqn{\tau \sim G(aa,bb)}). Default: \code{TRUE}.
#' @param aa Positive real giving the first parameter in the Gamma prior for \code{tau}. Default: \code{2}.
#' @param bb Positive real giving the second parameter in the Gamma prior for \code{tau}. Default: \code{1e-6}.
#' @param pen.order Integer giving the order of the roughness penalty. Default: \code{3}.
#'
#' @return A list containing :
#' \itemize{
#' \item{\code{lpost}, \code{lpost.ni} : \verb{ }}{value of the log-posterior based on the given small bin frequencies \code{n.i} and the tabulated sample moments.}
#' \item{\code{lpost.mj} : \verb{ }}{value of the log-posterior based on the big bin frequencies \code{degross.data$freq.j} and the tabulated sample moments.}
#' \item{\code{llik.ni} : \verb{ }}{multinomial log-likelihood based on the given small bin frequencies \code{n.i}.}
#' \item{\code{llik.mj} : \verb{ }}{multinomial log-likelihood based on the big bin frequencies \code{degross.data$freq.j}.}
#' \item{\code{moments.penalty} : \verb{ }}{log of the joint (asymptotic) density for the observed sample moments.}
#' \item{\code{penalty} : \verb{ }}{\eqn{\log p(\phi|\tau) + \log p(\tau)}.}
#' \item{\code{Score}, \code{Score.ni} : \verb{ }}{score (w.r.t. \eqn{\phi}) of \code{lpost.ni}.}
#' \item{\code{Score.mj} : \verb{ }}{score (w.r.t. \eqn{\phi}) of \code{lpost.mj}.}
#' \item{\code{Fisher} & \code{Fisher.ni}: \verb{ }}{information matrix (w.r.t. \eqn{\phi}) of \code{lpost.ni}.}
#' \item{\code{Fisher.mj} : \verb{ }}{information matrix (w.r.t. \eqn{\phi}) of \code{lpost.mj}.}
#' \item{\code{M.j} : \verb{ }}{theoretical moments of the density (resulting from \eqn{\phi}) within a big bin.}
#' \item{\code{pi.i} : \verb{ }}{small bin probabilities.}
#' \item{\code{ui} : \verb{ }}{small bin midpoints.}
#' \item{\code{delta} : \verb{ }}{width of the small bins.}
#' \item{\code{gamma.j} : \verb{ }}{Big bin probabilities.}
#' \item{\code{tau} : \verb{ }}{reminder of the value of the roughness penalty parameter \eqn{\tau}.}
#' \item{\code{phi} : \verb{ }}{reminder of the vector of spline parameters (defining the density).}
#' \item{\code{n.i} : \verb{ }}{reminder of the small bin frequencies given as input.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#' @seealso \code{\link{degross_lpostBasic}}, \code{\link{degross}}, \code{\link{degross.object}}.
#'
#' @examples
#' sim = simDegrossData(n=3500, plotting=TRUE,choice=2) ## Generate grouped data
#' obj.data = degrossData(Big.bins=sim$Big.bins, freq.j=sim$freq.j, m.j=sim$m.j)
#' print(obj.data)
#' obj.fit = degross(obj.data) ## Estimate the underlying density
#' ## Evaluate the log-posterior at convergence
#' res = with(obj.fit, degross_lpost(phi, tau, n.i, obj.data, diag.only=diag.only))
#' print(res$Score) ## Score of the log posterior at convergence
#'
#' @export
degross_lpost = function(phi,tau,n.i,degross.data,
                         use.moments=rep(TRUE,4), freq.min=20, diag.only=FALSE,
                         penalize=TRUE,aa=2,bb=1e-6,pen.order=3){
    ## Spline parameters
    K = length(phi)
    ##
    ## Penalty matrix
    Dd = diff(diag(K),diff=pen.order)
    Pd = t(Dd) %*% Dd
    ##
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
    temp = exp(eta.i-max(eta.i)) ## Small bin means
    pi.i = temp / sum(temp) ## Small bin probs
    ## gamma.j = c(tapply(pi.i,small.to.big,sum)) # Big bin probs
    gamma.j = c(rowsum(pi.i,small.to.big)) # Big bin probs
    ##
    ## Theoretical moments of the fitted density within a big bin
    M.j = matrix(nrow=J,ncol=8)
    rownames(M.j) = paste("Bin",1:J,sep="")
    colnames(M.j) = paste("mu",1:8,sep="")
    M.j[,1] = rowsum(ui*pi.i,small.to.big) / gamma.j # Theoretical mean within big bin j
    ## M.j[,1] = tapply(ui*pi.i,small.to.big,sum) / gamma.j # Theoretical mean within big bin j
    for (r in 2:8){
        ## kth theoretical central moment within each of the big bins
        M.j[,r] = rowsum((ui-M.j[small.to.big,1])^r*pi.i,small.to.big) / gamma.j ## <rowsum> much faster than <tapply> !!
        ## M.j[,r] = tapply((ui-M.j[small.to.big,1])^r*pi.i,small.to.big,sum) / gamma.j
    }
    ## Gradient
    Bi.tilde = t(t(B.i) - c(apply(pi.i*B.i,2,sum)))
    PiBi.tilde = pi.i*Bi.tilde
    ## Contribution due to bin frequencies
    Score.ni = Score.mj = Score.lprior = rep(0,ncol(B.i))
    ##
    Score.ni = c(t(B.i)%*%(n.i-m.tot*pi.i))
    for (j in 1:J){
        idx = which(small.to.big == j)
        Score.mj = Score.mj + freq.j[j]/gamma.j[j]*c(apply(PiBi.tilde[idx,],2,sum))
    }
    ## dmu1j.k = matrix(nrow=J,ncol=K) #  \partial\mu_{1j} / \partial\theta_k
    R = 4 ## Up to moments of order 4 of interest
    dmujr.k = array(dim=c(J,R,K)) #  \partial\mu_{rj} / \partial\theta_k
    ##
    ## Contribution due to sample central moments within each big bin
    if (sum(use.moments) > 0){
        invSigma.j = array(dim=c(J,4,4))
        for (j in 1:J){
            idx = which(small.to.big == j)
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
            ##
            ## Sample moments of interest with non-NA values AND big bin sample size >= freq.min
            idx2 = which((freq.j[j] >= freq.min) & (!is.na(m.j[j,])) & use.moments)
            ## if (!is.na(m.j[j,1])){
            if (length(idx2) > 0){   # i.e. if at least one observed sample moment of interest
                Mat = Sigma_fun(M.j[j,1:8]) # Theoretical covariance matrix for moments
                ## diag(Mat) = diag(Mat) + 1e-6
                if (diag.only) Mat = diag(diag(Mat)) # If desired, just keep variances (as covariances rely on moments of very high orders !)
                invSigma.j[j,,] = freq.j[j] * solve(Mat) # Theoretical version corrected for frequencies
                Mat = matrix(invSigma.j[j,idx2,idx2],nrow=length(idx2)) # Make sure that sub-matrix is a matrix...
                temp2 = matrix(dmujr.k[j,idx2,],nrow=length(idx2))
                temp = Mat %*% temp2
                Score.lprior = Score.lprior + matrix(m.j[j,idx2]-M.j[j,idx2],nrow=1) %*% temp
            }
        }
    }
    Score.lprior = c(Score.lprior)
    ## Contribution due to roughness penalty
    if (penalize){
        Score.lprior = Score.lprior - tau * c(Pd%*%phi)
    }
    ##
    Score.ni = Score.ni + Score.lprior
    Score.mj = Score.mj + Score.lprior
    ## Fisher (=expected) information matrix
    Fisher = Fisher.ni = Fisher.mj = Fisher.lprior = matrix(0,nrow=K,ncol=K)
    ##
    ## Contribution due to bin frequencies
    Bpi = c(t(B.i) %*% pi.i)
    Fisher = Fisher.ni = m.tot*t(B.i)%*%(pi.i*B.i) - m.tot*(Bpi %o% Bpi)
    Fisher.mj = Fisher.ni
    for (j in 1:J){
        idx = which(small.to.big == j)
        piBtilde.j = pi.i[idx]*Bi.tilde[idx,]
        cs = colSums(piBtilde.j)
        ## Fisher.mj = Fisher.mj + freq.j[j]/gamma.j[j]^2 * t(piBtilde.j)%*%piBtilde.j
        Fisher.mj = Fisher.mj + freq.j[j]/gamma.j[j]^2 * (cs %o% cs)
        Fisher.mj = Fisher.mj - freq.j[j]/gamma.j[j] * t(piBtilde.j)%*%Bi.tilde[idx,]
    }
    ## Fisher = m.tot*t(B.i) %*% (diag(pi.i)-pi.i%o%pi.i) %*% B.i ## Equivalent !!
    ##
    ## Contribution due to observed sample central moments within each big bin
    if (sum(use.moments) > 0){
        for (j in 1:J){
            ## Sample moments of interest with non-NA values AND big bin sample size >= 20
            idx2 = which((freq.j[j] >= freq.min) & (!is.na(m.j[j,])) & use.moments)
            if (length(idx2) > 0){   # i.e. if at least one observed sample moment of interest
                temp = matrix(invSigma.j[j,idx2,idx2],nrow=length(idx2))
                temp2 = matrix(dmujr.k[j,idx2,],nrow=length(idx2))
                Fisher.lprior = Fisher.lprior + t(temp2) %*% temp %*% temp2
                Fisher = Fisher + t(temp2) %*% temp %*% temp2
                ## Fisher = Fisher + (m.tot*gamma.j[j] / M.j[j,2]) * (c(dmu1j.k[j,]) %o% c(dmu1j.k[j,]))
            }
        }
    }
    ##
    ## Contribution due to roughness penalty
    if (penalize){
        Fisher = Fisher + tau*Pd
        Fisher.lprior = Fisher.lprior + tau*Pd
    }
    ##
    Fisher.ni = Fisher.ni + Fisher.lprior
    Fisher.mj = Fisher.mj + Fisher.lprior
    ## lpost
    llik.ni = sum(n.i*log(pi.i)) ## log complete likelihood
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
    lpost = llik.ni + moments.penalty + penalty    ## log complete posterior
    lpost.mj = llik.mj + moments.penalty + penalty ## log observed posterior
    ##
    Score = Score.ni
    ## Output
    res = list()
    res$lpost = lpost
    res$lpost.ni = lpost
    res$lpost.mj = lpost.mj
    res$llik.ni = llik.ni
    res$llik.mj = llik.mj
    res$moments.penalty = moments.penalty
    res$penalty = penalty
    res$Score = Score
    res$Score.ni = Score.ni
    res$Score.mj = Score.mj
    res$Score.lprior = Score.lprior
    res$Fisher = Fisher
    res$Fisher.ni = Fisher.ni
    res$Fisher.mj = Fisher.mj
    res$M.j = M.j
    res$pi.i = pi.i
    res$ui = ui ; res$delta = delta
    res$gamma.j = gamma.j
    res$tau = tau
    res$phi = phi
    res$n.i = n.i
    ##
    return(res)
}
