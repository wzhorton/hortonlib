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
  if(x == 0){
    if(log == TRUE){
      stop("x cannot be 0 if log is TRUE (negative infinity)")
    }
    else{
      return(0)
    }
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



