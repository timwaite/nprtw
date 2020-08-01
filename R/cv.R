#' Least Squares Cross-Validation criterion function
#'
#' This function computes the leave-one-out cross-validation criterion \eqn{CV(h)} for the trial bandwidths supplied in the vector \code{h}, where
#' \deqn{ CV(h) = \frac{1}{n}\sum_{i=1}^{n} (Y_i-\hat{m}^{(-i)}_h(x_i))^2 }
#' Used in \code{find_hcv} to find an optimal data-driven bandwidth.
#'
#' @seealso find_hcv
#'
#' @param h a vector of bandwidths at which to evaluate \eqn{CV(h)}
#' @param data the data, specified as a data frame containing two components \code{x} and \code{y}
#' @param estimator an R function for the estimator, with inputs and return type similar to \code{local_average} and \code{nw}
#' @param ... additional arguments to be passed to the estimator (e.g. \code{kernel}, \code{empty_nhood})
#'
#' @return The function returns a list with components
#'    \item{h}{the vector of bandwidths used}
#'    \item{cv}{a vector containing the computed values of \eqn{CV(h)}}
#' @export
#'
#' @examples
#'   # simulate and plot some data
#'   m <- function(x) (x^2+1)*sin(2*pi*x*((1-x) + 4*x))
#'   x <- sort(runif(100))
#'   y <- m(x) + rnorm(length(x), sd=0.1)
#'   simdata <- data.frame(x=x,y=y)
#'
#'   # create sequence of bandwidths, compute and plot CV function
#'   h <- seq(from=0,to=0.1,length=100)
#'   info <- lscv(h, simdata, nw, kernel=gauss)
#'   plot(info$h,info$cv, type="l", xlab="h", ylab="CV(h)")
#'
lscv <- function(h, data, estimator,  ...) {
  cv <- rep(0,length(h))
  n <- length(data$x)
  loo_data <- list()
  for (k in 1:length(h)) {
    for (i in 1:n) {
      loo_data$x <- data$x[-i]
      loo_data$y <- data$y[-i]
      mhat_xi_loo <- estimator(loo_data, h[k], data$x[i], ...)$mhat
      cv[k] <- cv[k] + (data$y[i] - mhat_xi_loo)^2 /n
    }
  }
  list(h=h,cv=cv)
}


#' Find optimal bandwidth using cross-validation
#'
#' The function finds a data-driven optimal bandwidth using cross-validation.
#' The output is an approximation to \eqn{\hat{h}_{CV}}, which is the minimum of the function
#' \deqn{ CV(h) = \frac{1}{n}\sum_{i=1}^{n} (Y_i-\hat{m}^{(-i)}_h(x_i))^2 }
#' where \eqn{\hat{m}^{(-i)}_h} denotes the leave-one-out estimator.
#'
#' The function is minimized approximately via a grid search.
#' The user specifies an interval over which to search for \eqn{h}.
#' The function then constructs a sequence of evenly-spaced trial bandwidths in that interval, of length \code{num_bws}, and computes \eqn{CV(h)} for each.
#' The best of these trial bandwidths is reported.
#'
#' By default, the function also produces a plot of the function \eqn{CV(h)}. This enables the user to check that a suitable interval has been specified.
#' The interval should be wide enough that it is clear that the identified point is a minimum. However, if it is too wide then the discretization error from the grid search may be substantial.
#'
#' @param data the data used to fit the estimator, a dataframe with columns \code{x} and \code{y}
#' @param estimator the estimator (\code{nw}, \code{local_average}, or another function with the same input and return types)
#' @param hrange a vector of length 2 specifying the range of \eqn{h}-values to try
#' @param num_bws number of different h-values to try in the range (default 100)
#' @param plot if set to \code{TRUE}, produces a plot of the CV function. Useful for checking that \code{hrange} is appropriate
#' @param ... additional arguments to pass to the estimator (e.g. \code{kernel}, \code{empty_nhood})
#'
#' @return A list with 4 components:
#' \item{hcv}{the identified optimal bandwidth}
#' \item{mincv}{the minimal value of \eqn{CV(h)}}
#' \item{h}{the vector of bandwidths that have been tried}
#' \item{cvs}{the values of \eqn{CV(h)} for the trial bandwidths}
#' @export
#'
#' @examples
#'   # simulate and plot some data
#'   m <- function(x) (x^2+1)*sin(2*pi*x*((1-x) + 4*x))
#'   x <- sort(runif(100))
#'   y <- m(x) + rnorm(length(x), sd=0.1)
#'   simdata <- data.frame(x=x,y=y)
#'
#'   find_hcv(simdata, nw, c(0,0.4))
find_hcv <- function(data, estimator, hrange=c(0,1), num_bws=100, plot=T, ...) {
  h <- seq(from=hrange[1],to=hrange[2],length=num_bws)
  cvs <- lscv(h,data,estimator,...)$cv
  hcv <- h[which.min(cvs)]
  mincv <- cvs[which.min(cvs)]
  if (plot) {
    plot(h, cvs, ylab="CV(h)", type="l")
    abline(v=hcv,lty=2)
  }
  list(hcv=hcv, mincv=mincv, h=h, lscv=cvs)
}
