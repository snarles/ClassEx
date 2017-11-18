#' Computing posterior likelihoods from predictive moments
#' 
#' Computes gaussian likelihoods for every (response, class) combo.
#' @param X_te Test class covariates
#' @param y_star Responses to be classified
#' @param pre_moments Predictive means and covariances for test classes
#' @param filt Logical vector: which responses to use
#' @param mc.cores Number of parallel cores
#' @return Matrix of log-likelihoods, rows = responses, columns = classes
#' @export
post_likes <- function(X_te, y_star, pre_moments,
                       filt = rep(TRUE, dim(y_star)[2]),
                       mc.cores = 0, ...) {
  y_filt <- y_star[, filt]
  L <- dim(X_te)[1]; n_te <- dim(y_filt)[1]; pY <- dim(y_filt)[2]
  colf <- function(i) {
    Mu <- pre_moments[[i]]$Mu
    Cov <- pre_moments[[i]]$Cov
    ld <- log(det(Cov))
    resid <- t(t(y_filt) - Mu)
    ss <- solve(Cov, t(resid))
    ips <- rowSums(resid * t(ss))
    -ld - ips
  }
  res <- mclapply0(1:L, colf, mc.cores = mc.cores)
  plikes <- do.call(cbind, res)
  plikes
}



#' Predictive distribution corresponding to MLE rule
#' 
#' A \code{backward_method}.
#' MLE rule is equivalent to using the same covariance for each response
#' @param X_te Test class covariates
#' @param B fitted coefficients
#' @param Sigma_e Fitted response covariance
#' @param filt Logical vector: which responses to use
#' @return A list with moments
#' @export
#' @examples
#' pars <- gen_params(n = 50, pY = 20, pX = 30, n_te = 100)
#' dat <- do.call(gen_data, pars)
#' zattach(dat)
#' res <- identification_pipeline1(X, Y, X_te, y_star, i_chosen,
#'   filter_method = no_filter, forward_method = fit_elnet_CV,
#'   Sigma_e_method = residual_offdiag, Sigma_t_method = assume_iid,
#'   backward_method = pre_mle, scoring_method = topk_score,
#'   mc.cores = 3)
#' res$score
pre_mle <- function(X_te, B, Sigma_e, filt, ...) {
  L <- dim(X_te)[1]
  ans <- as.list(numeric(L))
  Mus <- X_te %*% B
  for (i in 1:L) ans[[i]] <- list(Mu = Mus[i, ], Cov = Sigma_e)
  ans
}

#' Predictive distribution given true or estimated Sigma_b
#' 
#' A \code{backward_method}.
#' Computes a different covariance for each response due to posterior variance of B.
#' NOTE: To obtain Bayes rule (oracle) using pipeline, include
#' \code{backward_params = list(Sigma_b = Sigma_b, Sigma_e = Sigma_e, Sigma_t = Sigma_t)}
#' @param X design matrix
#' @param Y response matrix
#' @param X_te Test class covariates
#' @param Sigma_e Response covariance
#' @param Sigma_t Autocorrelation over time
#' @param Sigma_b Prior covariance of each row of B (rows assumed independent).
#' Constrained to be diagonal!
#' @param filt Logical vector: which responses to use
#' @return A list with moments
#' @export
#' @examples
#' pars <- gen_params(n = 50, pY = 20, pX = 30, n_te = 100)
#' dat <- do.call(gen_data, pars)
#' zattach(dat)
#' "Empirical Bayes"
#' ep_res <- identification_pipeline1(X, Y, X_te, y_star, i_chosen,
#'   filter_method = filter_eigenprism, forward_method = fit_elnet_CV,
#'   Sigma_e_method = residual_offdiag, Sigma_t_method = assume_iid,
#'   Sigma_b_method = use_eigenprism,
#'   backward_method = pre_Bayes,
#'   scoring_method = topk_score,
#'   mc.cores = 3)
#' ep_res$score
#' "Oracle Bayes rule"
#' bayesrule <- identification_pipeline1(X, Y, X_te, y_star, i_chosen,
#'   filter_method = no_filter, forward_method = fit_elnet_CV,
#'   Sigma_e_method = residual_offdiag, Sigma_t_method = assume_iid,
#'   backward_method = pre_Bayes,
#'   backward_params = list(Sigma_b = Sigma_b, Sigma_e = Sigma_e,
#'     Sigma_t = Sigma_t),
#'   scoring_method = topk_score,
#'   mc.cores = 3)
#' bayesrule$score
pre_Bayes <- function(X, Y, X_te, Sigma_e, Sigma_t, Sigma_b, filt, mc.cores = 0,...) {
  Y <- Y[, filt]
  pX <- dim(X)[2]; pY <- dim(Y)[2]
  B <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, computeCov = FALSE)
  res <- post_predictive(X, Y, X_te, Sigma_e, Sigma_b, Sigma_t, mc.cores = mc.cores)
  res
}

