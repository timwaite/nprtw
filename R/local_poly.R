
#' Local polynomial estimation
#'
#' Estimate a regression function using local polynomial estimation, essentially estimating a local Taylor series using locally weighted least squares.
#' The function requires the user to specify a bandwidth, \eqn{h}.

#' @inheritParams nw
#' @param degree the degree \eqn{p} of local polynomial to use. Defaults to \eqn{p=1} for local linear estimation.
#' @inherit nw return params
#' @return
#' @export
#'
#' @examples
#'   #  simulate and plot some data
#'   m <- function(x) (x^2+1)*sin(2*pi*x*((1-x) + 4*x))
#'   x <- sort(runif(100))
#'   y <- m(x) + rnorm(length(x), sd=0.1)
#'   simdata <- data.frame(x=x,y=y)
#'   plot(simdata)
#'
#'   # calculate the estimator at x=0.1, with bandwidth 0.02
#'   local_poly(simdata,h=0.02,t=0.1)
#'
#'   # a specialised print method has been provided to make life easier
#'   # however, we can still access the underlying numbers e.g.
#'
#'   fit <- local_poly(simdata,h=0.02,t=0.1)
#'   fit$mhat
#'   print(fit) # the same output as before
#'
#'   # plot the estimator with bandwidth 0.02 using default biweight kernel
#'   plot(local_poly(simdata,h=0.02))
#'
#'    # add a line for the estimator with bandwidth 0.4
#'   lines(local_poly(simdata,h=0.4), col=2)
#'
#'   # add a line for the estimator using Gaussian kernel
#'   lines(local_poly(simdata,h=0.02,kernel=gauss), col=4)
local_poly <- function(data,
                       h,
                       t=NULL,
                       kernel=biweight,
                       degree=1,
                       empty_nhood=NaN) {
  # if evaluation points not provided, form a plotting grid
  if (is.null(t)) {
    t <- seq(from=min(data$x),to=max(data$x),length=200)
    t <- sort(c(t, data$x-h, data$x+h, data$x+h-1e-6, data$x-h+1e-6))
  }
  m <- length(t)
  n <- length(data$x)
  A <- matrix(0,nrow=m,ncol=n)
  mhat <- rep(0,m)

  for (k in 1:m) {
    X <- outer(data$x, 0:degree, function(x,j) (x-t[k])^j)
    W <- diag( (1/h)*kernel((data$x-t[k])/h) )
    info <- t(X)%*%W%*%X
    tmp <- tryCatch(solve(info, t(X)%*%W), error=function(e) e)
    if (!identical(class(tmp),"matrix")) {
      mhat[k] <- empty_nhood
      A[k,] <- NaN
    } else {
      A[k,] <- tmp[1,]
    }
  }
  out <- list(t=t,h=h,mhat=A%*%data$y,A=A, data=data)
  class(out) <- "npfit"
  return(out)
}
