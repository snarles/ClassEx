#####################
## Sigma_e methods ## (Correlation between responses)
#####################

#' Covariance estimation using residuals
#' 
#' A \code{Sigma_e_method}.
#' Estimates within-row covariances (or \code{Sigma_e}) 
#' using sample covariance of residuals,
#' applying off-diagonal shrinkage.
#' @param X design matrix
#' @param Y response matrix
#' @param B fitted coefficient matrix
#' @param filt Logical vector: which responses to use
#' @param shrink shrinkage factor
#' @param mc.cores parallelization
#' @return Estimate for covariance matrix
#' @export
residual_offdiag <- function(X, Y, B, filt = rep(TRUE, dim(Y)[2]),
                             shrink = 0.5, mc.cores = 1, ...) {
  Y <- Y[, filt]
  Yh <- paramultiply(X, B, mc.cores)
  resids <- Y - Yh
  (1- shrink) * cov(resids) + shrink * diag(diag(cov(resids)))
}

#####################
## Sigma_t methods ## (Autocorrelation)
#####################

#' Produces an identity matrix
#' 
#' A \code{Sigma_t_method}.
#' @param X design matrix
#' @return Estimate for covariance matrix
#' @import pracma
#' @export
assume_iid <- function(X, ...) {
  eye(dim(X)[1])
}

#####################
## Sigma_b methods ## (Prior covariance of B)
#####################

#' Produces an identity matrix
#' 
#' A \code{Sigma_b_method}.
#' @param Y response matrix
#' @param sigma2 Scaling of matrix
#' @return Estimate for Sigma_b, a diagonal matrix
#' @import pracma
#' @export
assume_iid_B <- function(Y, sigma2 = 1.0, ...) {
  sigma2 * eye(dim(Y)[1])
}

#' Uses norms of a previously fitted B to set prior
#' 
#' A \code{Sigma_b_method}.
#' @param B fitted coefficient matrix
#' @return Estimate for Sigma_b, a diagonal matrix
#' @import pracma
#' @export
use_norms_of_Bhat <- function(B, ... ) {
  diag(colsums(B^2)/dim(B)[1])
}

#' Uses eigenprism to estimate prior covariances
#' 
#' A \code{Sigma_b_method}.
#' It is recommended to use filter_eigenprism in conjunction with this method.
#' @param X design matrix
#' @param Y response matrix
#' @param filt Logical vector: which responses to use
#' @return Estimate for Sigma_b, a diagonal matrix
#' @export
use_eigenprism <- function(X, Y, filt = rep(TRUE, dim(Y)[2]),... ) {
  res_EP <- eigenprisms(X, Y[, filt])
  T2 <- res_EP$T2
  lambdas_EP <-  T2/dim(X)[1]
  diag(lambdas_EP)
}






