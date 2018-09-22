#### utility.R ####


#' Determinant of Symmetric Positive Definite Matrix
#'
#' Computes the (log) determinant of a symmetric positive definite matrix using the cholesky
#' factorization.
#'
#' @param x square matrix
#' @param log logical; defaults to FALSE
#' @export

det_spd <- function(x, log = FALSE){
  out <- 2 * sum(log(diag(chol(x))))
  if(log == FALSE) out <- exp(out)
  return(out)
}


#' Spline Curve Interpolation
#'
#' Interpolates out equally spaced values by fitting a cubic spline on given points and
#' matching the endpoints. Often used to redefine a curve with warped time axis on an
#' even time basis. Note that quality of interpolation decreases with inadequate density
#' of defining points.
#'
#' @param x,y coordinate vectors on which to fit a cubic spline.
#' @param nout number of points to interpolate out
#' @export

interp_spline <- function(x, y, nout = length(y)) {
  ind_out <- seq(min(x), max(x), len = nout)
  spfit <- splinefun(x, y)
  return(spfit(ind_out))
}


#' Second Order Penalty Matrix
#'
#' Creates a second order penalty matrix used in fitting penalized b-splines.
#' Courtesy of CITATION NEEDED. Note that this function produces singular matrices.
#'
#' @param dim dimension of output matrix. Must be at least four.
#' @export

K2 <- function(dim) {
  tmp <- list(c(-2, rep(-4, dim - 3), -2), rep(1, dim))
  K <- Matrix::bandSparse(dim, k = -c(1:2), diag = tmp, symmetric = TRUE)
  diag(K) <- c(1, 5, rep(6, dim - 4), 5, 1)
  K <- matrix(as.numeric(K), nrow = dim, byrow = TRUE)
  K
}


#' First Order Penalty Matrix
#'
#' Creates a first order penalty matrix used in fitting penalized b-splines.
#' Courtesy of CITATION NEEDED. Note that this function produces singular matrices.
#'
#' @param dim dimension of output matrix. Must be at least three.
#' @export

K1 <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(1,rep(2,dim-2),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}


#' Acceptance Rate Calculator
#'
#' Calculates the acceptance rate of an MCMC chain by looking at the number of repeats.
#'
#' @param chain vector of mcmc values
#' @export

acc_rate <- function(chain) {
  n <- length(chain)
  accepts <- sapply(1:(n - 1), function(i) chain[i] != chain[i + 1])
  return(mean(accepts))
}


#' Stack a List of Vectors or Matrices
#'
#' Takes a list of vectors or matrices and returns one vector or matrix along with
#' a vector that specifies the first indices corresponding to the original elements.
#'
#' @param x list of vectors or matrices with the same number of columns
#' @return a list containing the stacked object and a vector of indices the specify where
#'   each of the original elements begins.
#' @export

stack <- function(x){
  n <- length(x)
  x <- lapply(x, as.matrix)

  lens <- sapply(x, nrow)
  inds <- c(1, sapply(1:n, function(i) 1 + sum(lens[1:i]))[-n])
  out <- list()
  out$stack <- do.call(rbind, x)
  out$inds <- inds
  return(out)
}


#' Unstack a Vector or Matrix
#'
#' Takes a vector or matrix and an index vector and returns a list containing the pieces.
#'
#' @param x vector or matrix.
#' @param inds vector indicating how to split x. If not provided then it will
#'   automatically attempt to split into equal parts.
#' @param n indicates the number of elements to split into.
#'   Only needed when inds is not provided.
#' @return list containing the split elements.
#' @export

unstack <- function(x, inds, n){
  x <- as.matrix(x)
  nr <- nrow(x)

  if(missing(inds)){
    if((nper_sub <- nr/n) %% 1 != 0) stop("the length of x is not evenly divided by n")
    inds <- ((1:n)-1)*nper_sub + 1
  }

  inds <- c(inds, nr + 1)

  return(lapply(1:(length(inds)-1), function(i) x[inds[i]:(inds[i+1] - 1),]))
}


