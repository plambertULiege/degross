#' Simulation of grouped data and their sample moments to illustrate the \link{degross} density estimation procedure
#' @usage simDegrossData(n, plotting=TRUE, choice=2, J=3)
#' @importFrom graphics axis curve legend
#' @importFrom stats dnorm pnorm rnorm dbeta pbeta rbeta dgamma pgamma rgamma rbinom
#'
#' @param n Desired sample size
#' @param plotting Logical indicating whether the histogram of the simulated data should be plotted. Default: FALSE
#' @param choice Integer in 1:3 indicating from which mixture of distributions to generate the data
#' @param J Number of big bins
#'
#' @return A list containing tabulated frequencies and central moments of degrees 1 to 4 for data generated using a mixture density. This list contains :
#' \itemize{
#' \item{\code{n} : \verb{ }}{total sample size.}
#' \item{\code{J} : \verb{ }}{number of big bins.}
#' \item{\code{Big.bins} : \verb{ }}{vector of length \code{J+1} with the big bin limits.}
#' \item{\code{freq.j} : \verb{ }}{vector of length \code{J} with the observed big bin frequencies.}
#' \item{\code{m.j} : \verb{ }}{\code{J} by \code{4} matrix with on each row the observed first four sample central moments within a given big bin.}
#' \item{\code{true.density} : \verb{ }}{density of the raw data generating mechanism (to be estimated from the observed grouped data).}
#' \item{\code{true.cdf} : \verb{ }}{cdf of the raw data generating mechanism (to be estimated from the observed grouped data).}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references
#' Lambert, P. (2023) Nonparametric density estimation and risk quantification from tabulated sample moments. Insurance: Mathematics and Economics, 108: 177-189. doi:10.1016/j.insmatheco.2022.12.004
#' Lambert, P. (2021) Moment-based density and risk estimation from grouped summary statistics. arXiv:2107.03883.
#'
#' @seealso \code{\link{degrossData}}.
#'
#' @examples
#' ## Generate data
#' sim = simDegrossData(n=3500, plotting=TRUE, choice=2, J=3)
#' print(sim$true.density) ## Display density of the data generating mechanism
#'
#' # Create a degrossData object
#' obj.data = with(sim, degrossData(Big.bins=Big.bins, freq.j=freq.j, m.j=m.j))
#' print(obj.data)
#'
#' @export
simDegrossData = function(n, plotting=TRUE, choice=2, J=3){
  if (!(choice %in% 1:3)){
    message("Distribution mixture <choice> should be in 1:3\n")
    return(NULL)
  }
  Big.bins = c(-1,round(seq(1,6,length=J),1))
  w1 = .20 ; w2 = 1-w1
  w = rbinom(n,1,w1) ## Indicator of the selected mixture component for a given unit
  n1 = round(w1*n) ; n2 = n - n1 ## Percentage of data per mixture component
  ## true.density: density corresponding to the data generating mechanism
  ## y: generated (ungrouped) raw data
  if (choice ==1){
    ## Normal mixture
      true.density = function(x) w1*dnorm(x,1,1/3) + w2*dnorm(x,3.5,.5)
      true.cdf = function(x) w1*pnorm(x,1,1/3) + w2*pnorm(x,3.5,.5)
      y = w*rnorm(n,1,1/3) + (1-w)*rnorm(n,3.5,.5)
  }
  if (choice ==2){
    ## Normal-Reverse Gamma mixture
      true.density = function(x) w1*dnorm(x,1,1/3) + w2*dgamma(5.6-x,11,6)
      true.cdf = function(x) w1*pnorm(x,1,1/3) + w2*(1-pgamma(5.6-x,11,6))
      y = w*rnorm(n,1,1/3) + (1-w)*(5.6-rgamma(n,11,6))
  }
  if (choice ==3){
      ## Normal-Beta mixture
      true.density = function(x) w1*dnorm(x,1,1/3) + w2*(1/4)*dbeta((x-2)/4,4,8)
      true.cdf = function(x) w1*pnorm(x,1,1/3) + w2*pbeta((x-2)/4,4,8)
      y = w*rnorm(n,1,1/3) + (1-w)*(2+4*rbeta(n,4,8))
  }
  ## Observed grouped data
  ## Big.bins = c(-1,6) ## Big bin limits
  ## Big.bins = c(-1,0,6) ## Big bin limits
  ## Big.bins = c(-1,0,4,6) ## Big bin limits
  ## Big.bins = c(-1,0,2,3,4,6) ## Big bin limits
  ## Big.bins = c(-1,0,2,3.5,4.5,6) ## Big bin limits
  ## Big.bins = c(-1,0,2,2.7,3.5,4.5,6) ## Big bin limits
  ## Big.bins = c(-1,0,2,2.7,3.5,4,4.5,6) ## Big bin limits
  # Big.bins = c(-1,2,4,6) ## Big bin limits <-----------
  # J = length(Big.bins) - 1 ## Number of Big bins
  temp = hist(y,breaks=Big.bins,plot=FALSE)
  if (plotting) {
    tab = hist(y,breaks=Big.bins,
               col="grey85",border="white",
               ylim=c(0,1.1*max(temp$density,.95*dnorm(3.5,3.5,.5))),
               plot=TRUE,xaxt='n',main="",xlab="") ## Observed histogram
    axis(1,at=Big.bins,labels=Big.bins)
    curve(true.density,add=TRUE,col="red",lwd=2)
    legend("topleft",legend=c("Observed freq.","Target density"),
           col=c("grey85","red"),lwd=c(10,2),lty="solid",
           box.lty=0, inset=.02)
  } else {
    tab = temp
  }
  freq.j = tab$counts ## Observed big bin frequencies
  ## Function to compute the observed sample central moments
  sample.moments= function(x,...){
    m = rep(NA,4)
    m[1] = mean(x,...)
    z = (x - m[1])
    if (length(x) > 1){
      m[2] = mean(z**2,...)
      m[3] = mean(z**3,...)
      m[4] = mean(z**4,...)
    }
    return(m)
  }
  ## Observed (first four) sample central moments
  m.j = matrix(NA,nrow=J,ncol=4)
  colnames(m.j) = paste("m",1:4,sep="")
  rownames(m.j) = paste("Bin",1:J,sep="")
  ##
  for (j in 1:J){ ## Loop over big bins
    xx = y[(y > Big.bins[j]) & (y <= Big.bins[j+1])] ## Data within jth big bin
    m.j[j,] = sample.moments(xx) ## Sample central moments within jth big bin
  }
  ##
  ans = list(n=n, J=J, Big.bins=Big.bins, freq.j=freq.j, m.j=m.j,
             true.density=true.density, true.cdf=true.cdf)
  return(ans)
}
