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
  if (any(shape <= 0) || any(scale <= 0)) {
    stop("Shape and rate/scale must be positive")
  }
  if (!is.logical(log)) {
    stop("log needs to be logical")
  }

  a <- shape
  b <- scale
  out <- a * log(b) - lgamma(a) + (-a - 1) * log(x) - b / x
  out[is.nan(out)] <- -Inf

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
  if(missing(prec)) {
    cov <- format_Matrix(cov, sparse = TRUE, symmetric = TRUE)
    return(as.numeric(mu + crossprod(chol(cov), rnorm(ncol(cov)))))
  } else if(missing(cov)) {
    prec <- format_Matrix(prec, sparse = TRUE, symmetric = TRUE)
    return(as.numeric(mu + solve(chol(prec), rnorm(ncol(prec)))))
  }
  stop("Provide either Precision or Covariance, but not both")
}


# #' Multivariate Normal Density Function
# #'
# #' Computes the (log) density of the multivariate normal distribution
# #' using either the covariance or precision parametrization.
# #'
# #' @param y vector of values.
# #' @param mu mean vector.
# #' @param cov covariance matrix.
# #' @param prec precision matrix.
# #' @param log logical; if TRUE, density calculations are computed on the log scale.
# #' @param unnorm logical; if TRUE then only density terms dependent on y are calculated
# #' @export
#
# dmnorm <- function(y, mu, cov, prec, log = FALSE, unnorm = FALSE) {
#
#   if(missing(prec)) {
#     cov <- format_Matrix(cov, sparse = TRUE, symmetric = TRUE)
#     prec <- chol2inv(chol(cov))
#   } else if(missing(cov)) {
#     prec <- format_Matrix(prec, sparse = TRUE, symmetric = TRUE)
#   } else {
#     stop("Provide either cov or prec, but not both")
#   }
#
#   n <- ncol(prec)
#   out <- - .5 * t(y - mu) %*% prec %*% (y - mu)
#   if(unnorm == FALSE){
#     if(!missing(cov)){
#       out <- out + -n/2 * log(2*pi) - .5 * det_spd(cov, log = TRUE)
#     } else {
#       out <- out + -n/2 * log(2*pi) + .5 * det_spd(prec, log = TRUE)
#     }
#   }
#
#   if (log == FALSE) out <- exp(out)
#   return(as.numeric(out))
# }
