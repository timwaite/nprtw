#' @rdname npfit-plot
#' @aliases  plot.npfit lines.npfit
#'
#' @title Plotting functions for npfit objects
#'
#' @description
#' Plotting methods for fitted nonparametric regression objects.
#' The \code{plot} command creates a new plot, with the observed data as points and the fitted nonparametric estimator as a line.
#' The \code{lines} command adds the nonparametric estimate to an existing plot, as a line.
#'
#' @param x  an \code{npfit} object, the result of running \code{nw} or \code{local_average}
#' @param ... additional arguments to be passed to the plotting function, e.g. \code{xlim}, \code{col}, \code{lty} etc.
#'
#' @examples
#'  m <- function(x) (x^2+1)*sin(2*pi*x*((1-x) + 4*x))
#'  x <- sort(runif(100))
#'  y <- m(x) + rnorm(length(x), sd=0.1)
#'  simdata <- data.frame(x=x,y=y)
#'  plot(local_average(simdata,h=0.02))
#'  lines(local_average(simdata,h=0.04), col=4)
#'

plot.npfit <- function(x,  ...) {
  plot(x$data, ...);
  ord <- order(x$t)
  lines(x$t[ord],x$mhat[ord])
}

#' @export
#' @rdname npfit-plot
lines.npfit <- function(x, ...) {
  ord <- order(x$t)
  lines(x$t[ord],x$mhat[ord], ...)
}


#' Print npfit object
#'
#' Pretty printing for nonparametric fits from \code{nw} or \code{local_average} etc.
#'
#' @param x an object of class \code{npfit}
#' @param ... additional arguments (not needed)
#' @export
#'
print.npfit <- function(x, ...) {
  cat("Nonparametric estimate (n=",length(x$data$x),", bandwidth=",x$h,")\n",sep="");
  print(data.frame(t=x$t, mhat=x$mhat))
}
