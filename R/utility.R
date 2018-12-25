#### utility.R ####

#' Imbue Matrix With Property Classes
#'
#' Gives a matrix classes such as sparseness and symmetry from the Matrix package. Note that
#' a matrix cannot be both symmetric and triangular unless it is diagonal.
#'
#' @param x a matrix with (or without) zero values
#' @param sparse logical; gives class sparseMatrix and is almost always true.
#' @param symmetric logical; gives class symmetricMatrix
#' @param triangular logical; gives class triangularMatrix
#' @param diagonal logical; gives class diagonalMatrix
#' @export

format_Matrix <- function(x, sparse = TRUE, symmetric = FALSE, triangular = FALSE, diagonal = FALSE){
  if(sparse && !is(x, "sparseMatrix")){
    x <- as(x, "sparseMatrix")
  }
  if(symmetric && !is(x, "symmetricMatrix")){
    x <- as(x, "symmetricMatrix")
  }
  if(triangular && !is(x, "triangularMatrix")){
    x <- as(x, "triangularMatrix")
  }
  if(diagonal && !is(x, "diagonalMatrix")){
    x <- as(x, "diagonalMatrix")
  }
  return(x)
}

#' Determinant of Symmetric Positive Definite Matrix
#'
#' Computes the (log) determinant of a symmetric positive definite matrix using the cholesky
#' factorization.
#'
#' @param x square matrix
#' @param log logical; defaults to FALSE
#' @export

det_spd <- function(x, log = FALSE){
  x <- format_Matrix(x, sparse = TRUE, symmetric = TRUE)
  out <- 2 * sum(log(diag(chol(x))))
  if(log == FALSE) out <- exp(out)
  return(out)
}

#' Efficient Computation of Quadratic Form
#'
#' Evaluates the quadratic form x'Ax using optimized cross products.
#' An upper triangular matrix with U'U = A can also be given to dramatically increase speed.
#'
#' @param x a vector
#' @param A a matrix. Only one of A and U needs to be provided.
#' @param U an upper triangular matrix. Note that plugging in chol(A) is slower.
#' @export

qform <- function(x, A, U){
  if(missing(A)){
    U <- format_Matrix(U, sparse = TRUE, triangular = TRUE)
    return(crossprod(U %*% x))
  } else if(missing (U)) {
    A <- format_Matrix(A, sparse = TRUE)
    return(crossprod(crossprod(A,x),x))
  }
  stop("Provide either A or U but not both")
}

#' Acceptance Rate Calculator
#'
#' Calculates the acceptance rate of an MCMC chain by looking at the number of repeats.
#'
#' @param chain vector of mcmc values
#' @export

acc_rate <- function(chain) {
  return(1 - mean(diff(chain) == 0))
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

#' Stack A List of Matrices
#'
#' Stack a list of matrices together into one tall matrix. This is a simplified version of the
#' stack function. It quickly stacks matrices and gives the result sparseMatrix class.
#'
#'  @param x list of matrices.
#'  @export

stack_Matrix <- function(x){
  out <- do.call(rbind, x)
  return(format_Matrix(out, sparse = TRUE))
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

#' Square Root Matrix using Eigen Values/Vectors
#'
#' Compute a square root matrix using the eigenvalue/vector spectral decomposition.
#'
#' @param x semi-positive definite matrix
#' @param symmetric logical; indicates to eigen if x is symmetric
#' @export

sqrt_eigen <- function(x, symmetric = FALSE) {
  eig <- eigen(x, symmetric)
  C <- eig$vectors
  D <- diag(eig$values)
  return(C %*% D^(.5) %*% t(C))
}

#' Check is vector is monotone
#'
#' Determines if the sequence of vector elements are increasing or decreasing with the option to
#' specify strictness.
#'
#' @param x numeric vector
#' @param strict logical; if true then strict monotonicity is determined
#' @export

is_monotone <- function(x, strict){
  ds <- diff(x)
  if(all(ds > 0) || all(ds < 0)){
    return(TRUE)
  }
  else {
    if(strict == FALSE) {
      if(all(ds >= 0) || all(ds <= 0)){
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    }
    else {
      return(FALSE)
    }
  }
}

#' Project a vector into monotone space
#'
#' Makes a vector monotone. (REFERENCE) proposed a method to project a vector onto monotone space.
#'
#' @param x non-empty numeric vector
#' @param type character indicating if the function is increasing or decreasing
#' @param forced numeric vector of indeces specifying points that cannot be moved
#' @export

monotonize <- function(x, type = "increasing", forced = NULL){
  n <- length(x)

  if(type == "increasing"){
    x <- x
  } else if(type == "decreasing"){
    x <- rev(x)
  } else {
    stop("Type must either be 'increasing' or 'decreasing'")
  }

  for(i in 2:n){
    if(x[i] < x[i-1]){
      j <- 1
      lagmean <- lagsum <- x[i]
      while(x[i-j] > lagmean){
        lagsum <- lagsum + x[i-j]
        j <- j + 1
        lagmean <- lagsum/j
        forced_check <- na.omit(match(((i-j+1):i), forced))
        if(length(forced_check) != 0){
          if(length(forced_check) != 1) {stop("SERIOUS ISSUES")}
          lagmean <- x[forced[forced_check]]
          lagsum <- j*lagmean
        }
        if(i == j) break
      }
      x[(i-j+1):i] <- lagmean
    }
  }

  if(type == "increasing") return(x)
  if(type == "decreasing") return(rev(x))
}

#' Fast version of Matrix :: .bdiag()
#'
#' Fast version of Matrix :: .bdiag() -- for the case of *many*  (k x k) matrices. Comes from example
#' section of Matrix::bdiag documentation.
#'
#' @param lmat list(<mat1>, <mat2>, ....., <mat_N>)  where each mat_j is a  k x k 'matrix'
#' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
#' @export

bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

