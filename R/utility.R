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

