

#' Uniform kernel
#'
#' Evaluates the uniform kernel, \deqn{K(u)=\frac{1}{2} I[|u|\leq 1].}
#'
#' @family kernel functions
#' @param u a vector of points at which to evaluate the kernel.
#'
#' @return a vector giving the evaluations.
#' @export
#'
uniform = function(u) { dunif(u,-1,1) }

#' Biweight kernel
#'
#' Evaluates the biweight kernel, \deqn{K(u)=\frac{15}{16} (1-u^2)^2 I[|u|\leq 1].} This is a highly efficient kernel which is continuous and has a continuous first derivative.
#'
#' @family kernel functions
#'
#' @inheritParams uniform
#' @inherit uniform return
#' @export
#'
biweight = function(u) { ifelse(abs(u)<=1, (15/16)*(1-u^2)^2, 0) }



#' Epanenchnikov kernel
#'
#' Evaluates the Epanechnikov kernel, \deqn{K(u)=\frac{3}{4} (1-u^2) I[|u|\leq 1].} This is the asymptotically optimal second order kernel. It is continuous, but the first derivative is discontinuous at \eqn{u=-1,1}.
#
#' @family kernel functions
#'
#' @inheritParams uniform
#' @inherit uniform return
#' @export
#'
epanechnikov = function(u) { ifelse(abs(u)<=1, (3/4)*(1-u^2), 0) }

#' Gaussian kernel
#'
#' Evaluates the Gaussian kernel, equal to the standard normal density. This is an alias for \code{dnorm}.
#'
#' @family kernel functions
#'
#' @inheritParams uniform
#' @inherit uniform return
#' @export
#'
gauss = function(u) { dnorm(u) }
