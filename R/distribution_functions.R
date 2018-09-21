#### distribution_functions.R ####


#' Inverse Gamma Density Function
#'
#' Computes the (log) density of the inverse gamma distribution using
#' either the scale or rate parametrization. The scale parametrization has
#' expected value of scale/(shape - 1)
#'
#' @param x vector of nonnegative values.
#' @param shape,scale shape and scale parameters. Must be stricly positive.
#' @param rate alternative way to specify scale.
#' @param log logical; if TRUE, density calculations are computed on the log scale.
#' @export

dinvgamma <- function(x, shape, rate, scale = 1 / rate, log = FALSE) {
  if (shape <= 0 || scale <= 0 || x < 0) {
    stop("Shape, rate, scale, and x must all be positive")
  }
  if (!is.logical(log)) {
    stop("log needs to be logical")
  }

  a <- shape
  b <- scale
  if(x == 0){
    out <- 0
  }
  else{
    out <- a * log(b) - lgamma(a) + (-a - 1) * log(x) - b / x
  }
  if (log == FALSE) out <- exp(out)

  return(out)
}


#' Random Multivariate Normal Generator
#'
#' Generates a normally distributed vector given a mean vector and
#' either a covariance matrix or a precision matrix. Computationally, this method
#' has an advantage since an inverse is never taken.
#'
#' @param mu mean vector.
#' @param cov covariance matrix.
#' @param prec precision matrix.
#'
#' @return A randomly generated vector with the same length as mu.
#' @export

rmnorm <- function(mu, cov, prec) {
  if (missing(mu)) stop("Provide a mean vector")
  if (!xor(missing(prec), missing(cov))) {
    stop("Provide either Precision or Covariance, but not both")
  }

  if (missing(prec)) {
    out <- mu + t(chol(cov)) %*% rnorm(length(mu))
  } else {
    out <- mu + backsolve(chol(prec), rnorm(length(mu)))
  }
  return(as.numeric(out))
}


