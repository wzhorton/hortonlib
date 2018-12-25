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
  if(missing(R)){
    R_inv <- format_Matrix(R_inv, sparse = TRUE, symmetric = TRUE)
    qf <- qform(y - mu, R_inv)
  } else if(missing(R_inv)) {
    R <- format_Matrix(R, sparse = TRUE, symmetric = TRUE)
    qf <- crossprod(solve(chol(R), y - mu))
  } else {
    stop("Provide either R or R_inv, but not both")
  }
  return(1 / rgamma(1, .5 * length(y) + a, rate = as.numeric(.5 * qf + b)))
}


#' Conjugate Multivariate Normal - Multivariate Normal Update
#'
#' Generates a value from the posterior distribution in the case where
#' there is a multivariate normal likelihood and a multivariate normal prior.\cr \cr
#' Argument Model:\cr
#' y ~ Nn(X*beta, Sig)\cr
#' beta ~ Np(mu, V)\cr
#' @param y vector of values at the likelihood level.
#' @param X fixed design matrix in likelihood.
#' @param mu prior mean vector.
#' @param Sig,Sig_inv likelihood covariance/precision matrix.
#' @param V,V_inv prior covariance/precision matrix.
#' @export

update_normal_normal <- function(y, X, mu, Sig, V, Sig_inv, V_inv) {
  if(missing(Sig_inv)){
    Sig <- format_Matrix(Sig, sparse = TRUE, symmetric = TRUE)
    Sig_inv <- chol2inv(chol(Sig))
  } else if(missing(Sig)){
    Sig_inv <- format_Matrix(Sig_inv, sparse = TRUE, symmetric = TRUE)
  } else {
    stop("Provide either Sig or Sig_inv, but not both")
  }

  if(missing(V_inv)){
    V <- format_Matrix(V, sparse = TRUE, symmetric = TRUE)
    V_inv <- chol2inv(chol(V))
  } else if(missing(V)){
    V_inv <- format_Matrix(V_inv, sparse = TRUE, symmetric = TRUE)
  } else {
    stop("Provide either V or V_inv, but not both")
  }

  X <- format_Matrix(X, sparse = TRUE)

  vv <- qform(X, Sig_inv) + V_inv
  vterm <- chol2inv(chol(vv))
  return(rmnorm(vterm %*% (t(X) %*% (Sig_inv %*% y) + V_inv %*% mu), prec = vv))
}


#' Gaussian Process Update
#'
#' Updates and returns the mean vector for a gaussian process given evaluation points,
#' observed data, a mean function, and a covariance function.\cr \cr
#' Argument model:\cr
#' y(time) ~ GP(mnfun(x),covfun(x1-x2))
#'
#' @param x,y coordinate vectors for the observed data points.
#' @param time locations to evaluate the curve at.
#' @param mnfun mean function that takes a single argument.
#'   If using a constant like 0 use function(z)\{0\}.
#' @param covfun covariance function that takes a single argument that represents
#'   an absolute distance. Distances are internally calculated using fields::rdist.
#' @param random logical; FALSE will return the updated mean function.
#'   TRUE will generate a random curve.
#' @return a vector corresponding to the time variable that represents the updated
#'   mean vector.
#' @export

update_gaussian_process <- function(x, y, time, mnfun, covfun, random = FALSE) {
  R <- format_Matrix(covfun(fields::rdist(c(time,x))), sparse = TRUE, symmetric = TRUE)
  m <- length(time)
  f <- length(x)

  R11 <- R[1:m,1:m]
  R22 <- R[(m+1):(m+f),(m+1):(m+f)]
  R22i <- chol2inv(chol(R22))
  R12 <- R[1:m,(m+1):(m+f)]

  mu1 <- mnfun(time)
  mu2 <- mnfun(x)
  up_mean <- as.numeric(mu1 + R12 %*% (R22i %*% (y - mu2)))
  up_var <- R11 - R12 %*% tcrossprod(R22i, R12)
  if(random == FALSE){
    return(list(up_mean = up_mean, up_var = as.matrix(up_var)))
  } else {
    return(rmnorm(up_mean, cov = up_var))
  }
}



