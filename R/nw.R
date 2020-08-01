#' Nadaraya-Watson estimator
#'
#' Estimate a regression function using the Nadaraya-Watson estimator, using a user-specified bandwidth \eqn{h}.
#'
#' @inheritParams local_average
#' @param kernel a kernel function. The package supplies \code{uniform}, \code{gauss}, \code{epanechnikov} and \code{biweight} (the default).
#'
#' @inherit local_average return params
#' @export
#'
#' @examples
#'   # simulate and plot some data
#'   m <- function(x) (x^2+1)*sin(2*pi*x*((1-x) + 4*x))
#'   x <- sort(runif(100))
#'   y <- m(x) + rnorm(length(x), sd=0.1)
#'   simdata <- data.frame(x=x,y=y)
#'   plot(simdata)
#'
#'   # calculate the estimator at x=0.1, with bandwidth 0.02
#'   nw(simdata,h=0.02,t=0.1)
#'
#'   # a specialised print method has been provided to make life easier
#'   # however, we can still access the underlying numbers e.g.
#'
#'   fit <- nw(simdata,h=0.02,t=0.1)
#'   fit$mhat
#'   print(fit) # the same output as before
#'
#'   # plot the estimator with bandwidth 0.02 using default biweight kernel
#'   plot(nw(simdata,h=0.02))
#'
#'    # add a line for the estimator with bandwidth 0.4
#'   lines(nw(simdata,h=0.4), col=2)
#'
#'   # add a line for the estimator using Gaussian kernel
#'   lines(nw(simdata,h=0.02,kernel=gauss), col=4)
#'
#'   # NB the first plot is equivalent to the following:
#'   fit <- nw(simdata,h=0.02)
#'   plot(fit$data)
#'   lines(fit$t,fit$mhat)
#'
#'   # get smoother matrix
#'   fit$A
nw <- function(data,
               h,
               t=NULL,
               kernel=biweight,
               empty_nhood=NaN) {
  if (is.null(t)) {
    t <- seq(from=min(data$x),to=max(data$x),length=200)
    t <- sort(c(t, data$x-h, data$x+h, data$x+h-1e-6, data$x-h+1e-6))
    # NB the purpose of this elaborate evaluation grid is so that any discontinuities
    # are shown accurately in plots
    # relies on the support of the kernel being [-1,1]
  }
  m <- length(t)
  n <- length(data$x)
  A <- matrix(0,nrow=m,ncol=n)
  A <- outer(t, data$x, function(s,t) kernel((s-t)/h))
  den <- rowSums(A)
  A <- t(apply(A,1,function(x) x/sum(x)))
  out <- list(t=t, h=h, mhat=ifelse(den==0, empty_nhood, A%*%data$y), data=data, A=A)
  class(out) <- "npfit"
  out
}


#' Local average estimator
#'
#' Estimate a regression function m(x) by local averaging. To calculate the estimate, we average the response values for design points within distance \eqn{h} of \eqn{x}. The quantity \eqn{h} is known as the bandwidth.
#'
#' The function calls \code{nw} using the uniform kernel.
#'
#' @param data the data used to fit the estimator. Must be a data frame with columns \code{x} and \code{y}, where \code{x} contains the design points \eqn{x_1,\ldots,x_n}
#'  and \code{y} contains the response values \eqn{Y_1,\ldots,Y_n}
#' @param h a scalar giving the user-specified bandwidth (N.B. the cross-validation bandwidth can be computed using \code{find_hcv})
#' @param t (optional) a vector of points at which the estimator is evaluated. If unspecified, a sequence of 200 points is created that spans the range of the x-values in the data.
#' @param empty_nhood a scalar specfying a custom value to be returned at locations where the estimator is undefined (as occurs when there are no nearby data points to average).
#'   Default is \code{NaN}.
#' @return An object of class \code{npfit}, which is a list with 5 items:
#'   \item{t}{the vector of evaluation points}
#'   \item{h}{the bandwidth used}
#'   \item{mhat}{evaluations of the estimator \eqn{\hat{m}(t_1),\ldots, \hat{m}(t_n)}}
#'   \item{data}{the data used to fit the estimator}
#'   \item{A}{the smoother matrix, such that \eqn{\hat{m}=AY}.}
#'   Specialised \code{print}, \code{plot}, and \code{lines} methods are available for these objects, to facilitate analysis. See examples below.
#' @export
#'
#' @examples
#'   # simulate and plot some data
#'   m <- function(x) (x^2+1)*sin(2*pi*x*((1-x) + 4*x))
#'   x <- sort(runif(100))
#'   y <- m(x) + rnorm(length(x), sd=0.1)
#'   simdata <- data.frame(x=x,y=y)
#'   plot(simdata)
#'
#'   # calculate the estimator at x=0.1, with bandwidth 0.02
#'   local_average(simdata,h=0.02,t=0.1)
#'
#'   # a specialised print method has been provided to make life easier
#'   # however, we can still access the underlying numbers e.g.
#'
#'   fit <- local_average(simdata,h=0.02,t=0.1)
#'   fit$mhat
#'   print(fit) # the same output as before
#'
#'   # plot the estimator with bandwidth 0.02
#'   plot(local_average(simdata,h=0.02))
#'
#'    # add a line for the estimator with bandwidth 0.4
#'   lines(local_average(simdata,h=0.4), col=2)
#'
#'   # NB the first plot is equivalent to the following:
#'   fit <- local_average(simdata,h=0.02)
#'   plot(fit$data)
#'   lines(fit$t,fit$mhat)
#'
#'   # get smoother matrix
#'   fit$A
local_average <- function(data,
                          h,
                          t=NULL,
                          empty_nhood=NaN) { nw(data,h,t,uniform,empty_nhood) }
