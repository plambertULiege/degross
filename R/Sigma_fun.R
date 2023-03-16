#' Variance-covariance of sample central moments (root-n approximation)
#'  given the vector mu with the theoretical moments of order 1 to 8.
#'  CAREFUL: the result must be divided by n (= sample size)!
#'
#' @param mu Vector of length 8 with the first 8 theoretical central moments.
#'
#' @return Variance-covariance matrix of the first four sample central moments (CAREFUL: a division by the sample size is further required !)
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#' @examples
#' mu = numeric(8)
#' dfun = function(x) dgamma(x,10,5)
#' mu[1] = integrate(function(x) x*dfun(x),0,Inf)$val
#' for (j in 2:8) mu[j] = integrate(function(x) (x-mu[1])^j*dfun(x),0,Inf)$val
#' Sigma_fun(mu)
#'
#' @export
Sigma_fun= function(mu){
  m = mu
  #
  Lambda0 = matrix(c(m[2:5], m[3],m[4]-m[2]^2,m[5]-m[2]*m[3],m[6]-m[2]*m[4], m[4],m[5]-m[3]*m[2],m[6]-m[3]^2,m[7]-m[3]*m[4], m[5],m[6]-m[2]*m[4],m[7]-m[3]*m[4],m[8]-m[4]^2),byrow=4,nrow=4)
  #
  ## G0 = matrix(c(-1,0,0,0, 0,-1,0,0, -3*m[2],0,-1,0, -4*m[3],0,0,-1),nrow=4,byrow=T)
  G0 = matrix(c(-1,0,0,0, 0,-1,0,0, -3*m[2],0,-1,0, -4*m[3],0,0,-1),nrow=4,byrow=T)
  ##
  invG0 = solve(G0)
  Sigma = invG0%*%Lambda0%*%t(invG0)
  Sigma = .5*(Sigma + t(Sigma))
  attr(Sigma,"mu") = mu
  return(Sigma)
}
