#' Apply no filtering
#' 
#' A \code{filter_method}.
#' @export
no_filter <- function(X, Y, B, ...) rep(TRUE, dim(Y)[2])

#' Filters out columns with estimated 0 signal size
#' 
#' A \code{filter_method}.
#' @param X design matrix
#' @param Y response matrix
#' @return \code{filt}, logical vector
#' @export
filter_eigenprism <- function(X, Y, ... ) {
  res_EP <- eigenprisms(X, Y)
  T2 <- res_EP$T2
  T2 > 0
}
