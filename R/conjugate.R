#### conjugate.R ####


#' Conjugate Normal - Inverse Gamma Update
#'
#' Generates a value from the posterior distribution in the case where there
#' is a multivariate normal likelihood and an inverse gamma prior.\cr \cr
#' Argument model:\cr
#' y ~ Nn(mu, sig2*R)\cr
#' sig2 ~ IG(a,b)\cr
#' @param y vector of values at the likelihood level.
#' @param a prior shape value for inverse gamma.
#' @param b prior SCALE value for inverse gamma.
#' @param mu mean vector for multivariate normal.
#' @param R correlation matrix for multivariate normal.
#' @param R_inv scaled precision matrix, alternative specification to R.
#' @export

update_normal_invgamma <- function(y, a, b, mu, R, R_inv) {
  if (!xor(missing(R), missing(R_inv))) stop("Provide either R or R_inv, but not both")
  if (missing(R_inv)) R_inv <- chol2inv(chol(R))
  return(1 / rgamma(1, .5 * length(y) + a,
                    scale = .5 * t(y - mu) %*% R_inv %*% (y - mu) + b))
}


