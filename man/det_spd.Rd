% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{det_spd}
\alias{det_spd}
\title{Determinant of Symmetric Positive Definite Matrix}
\usage{
det_spd(x, log = FALSE)
}
\arguments{
\item{x}{square matrix}

\item{log}{logical; defaults to FALSE}
}
\description{
Computes the (log) determinant of a symmetric positive definite matrix using the cholesky
factorization.
}
