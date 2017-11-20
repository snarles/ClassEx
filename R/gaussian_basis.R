
#' Log sum exponential
#' 
#' @export
logsumexp <- function(v) log(sum(exp(v - max(v)))) + max(v)

#' Log mean exponential
#' 
#' @export
logmeanexp <- function(v) logsumexp(v - log(length(v)))


#' Formula for Gaussian-kernel based extrapolation basis functions
#' 
#' @param mu the mean of the basis function.
#' @param sigma2 the variance parameter of the basis function.
#' @param ks the number of classes for which we want to compute the basis function.
#' @param n.sample number of sample points for computing Gaussian integrals
gaussian_basis <- function(mu, sigma2, ks, n.sample = 1e4) {
  if (mu == Inf) return(rep(1, length(ks)))
  samp_qs <- (1:n.sample - 0.5)/n.sample
  samp <- qnorm(samp_qs, mu, sqrt(sigma2))
  logps <- pnorm(samp, log.p = TRUE)
  Logps <- repmat(logps, length(ks), 1)
  Ks <- repmat(t(t(ks - 1)), 1, n.sample)
  ans <- apply(Logps * Ks, 1, logmeanexp)
  exp(ans)
}

#' Get a Gaussian-kernel based basis matrix for ClassExReg.
#' 
#' @param max.mu maximum mean to use in the basis
#' @param kernel_sigma the sd. dev parameter of the basis function.
#' @param Ktrain the number of classes k for which we have computed average test accuracies ATA_k.
#' @param Ktarg The number of classes k for which we want to extrpolate the accuracy AGA_k (can be a vector).
#' @param n.sample number of sample points for computing Gaussian integrals
#' @export
get_gaussian_basis_mat <- function(max.mu, kernel_sigma, Ktrain, Ktarg, n.sample = 1e4) {
  (n_half <- ceiling(max.mu/kernel_sigma))
  (seq_half <- seq(0, max.mu, length.out = n_half))
  (kernel_mus <- c(rev(seq_half[-1]), seq_half, Inf))
  Xmat <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, Ktrain, n.sample))
  Xtarg <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, Ktarg, n.sample))
  list(Xmat = Xmat, Xtarg = Xtarg)
}

#' Get predictions for a list of basis sets
#' 
#' @param accs_sub subsampled accuracies ATA_k
#' @param basis_sets a list of bases obtained e.g. using lapply and get_gaussian_basis_mat
#' @param nonnegative constrain coefficients of basis to be nonnegative? default TRUE.
#' @export
bdwid_all_preds <- function(accs_sub, basis_sets, nonnegative = TRUE) {
  ans <- matrix(NA, length(basis_sets), nrow(basis_sets[[1]]$Xtarg))
  ntr <- length(accs_sub)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    if (nonnegative) {
      ls_fit <- nnls(Xmat, accs_sub)
    } else {
      ls_fit <- solve(t(Xmat) %*% Xmat, t(Xmat) %*% accs_sub)
    }
    ans[j, ] <- set1$Xtarg %*% ls_fit$x
  }
  ans
}