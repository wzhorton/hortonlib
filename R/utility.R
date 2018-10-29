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
