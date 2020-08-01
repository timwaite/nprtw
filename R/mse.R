#' Compute exact MSE of a kernel estimator
#'
#' This function computes the exact conditional bias, variance, and MSE of a kernel estimator. 'Conditional' means conditional given the design.
#' The estimator must be a linear smoother from this package.
#' The true regression function and response variance must both be known, therefore this function is only useful for theoretical examples.
#'
#' @param x_data a vector giving the \eqn{x}-values in the data set
#' @param h a scalar, giving the bandwidth
#' @param m an R function, corresponding to the true regression function
#' @param s2 a scalar giving the true response variance
#' @param estimator the estimator used (\code{nw}, \code{local_average}, or another estimator with the same input and output type)
#' @param t (optional) a vector specifying the points at which to evaluate the estimator. The default is a sequence of length 100 spanning the range of the \eqn{x}-data.
#' @param ... additional arguments to be passed to \code{estimator} (e.g. \code{kernel} or \code{empty_nhood}).
#'
#' @return A list with 4 components:
#' \item{t}{a vector containing the evaluation points used}
#' \item{bias}{a vector containing the bias of the estimator at the evaluation points}
#' \item{var}{a vector containing the variance of the estimator at the evaluation points}
#' \item{mse}{a vector containing the mean-squared error of the estimator at the evaluation points}
#' @export
#'
#' @examples
#' m <- function(x) { (x^2+1)*sin(2*pi*x*((1-x) + 4*x)) }
#' x <- sort(runif(80))
#' res <- mse(x, 0.02, m, 0.1^2, nw)
#' par(mfrow=c(1,2))
#' plot(res$t, res$mse, type="l", xlab="x", ylab="MSE")
#' plot(res$t, res$var, type="l", xlab="x", ylab="Variance")
#' rug(x)
#'
mse <- function(x_data, h, m, s2, estimator, t=seq(from=min(x_data),to=max(x_data),length=100), ...) {
  m_data <- sapply(x_data,m)
  m_test <- sapply(t,m)
  fit_using_m_data <- estimator(data.frame(x=x_data,y=m_data),h,t,...)
  bias <- fit_using_m_data$mhat - m_test
  sum_W2 <- rowSums(fit_using_m_data$A^2)
  var <-  ifelse(sum_W2>0, s2*sum_W2, Inf)
  list(t=t,bias=bias,var=var,mse=bias^2+var)
}

