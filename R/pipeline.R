#' Automates fitting forward/backward models and scoring
#' 
#' Pipeline 1 involves a single variable-selection (filter) pass over the data
#' and a single forward model fitting step, covariance fitting step,
#' and backwards modeling step.
#' @param X design matrix
#' @param Y response matrix
#' @param X_te Test class covariates
#' @param y_star Responses to be classified
#' @param i_chosen True classes
#' @param filter_method (function) Which variable selection method 
#' @param filter_params (list) Additional arguments for variable selection
#' @param forward_method (function) Which forward model 
#' @param forward_params (list) Additional arguments for forward model
#' @param Sigma_e_method (function) Which covariance estimation method
#' @param Sigma_e_params (list) Additional arguments for covariance estimation
#' @param Sigma_t_method (function) Which covariance estimation method
#' @param Sigma_t_params (list) Additional arguments for covariance estimation
#' @param Sigma_b_method (function) Empirical bayes heuristic for choosing B prior
#' @param Sigma_b_params (list) Additional arguments for Sigma_b estimation
#' @param backward_method (function) Which forward model 
#' @param backward_params (list) Additional arguments for forward model
#' @param scoring_method (function) Which forward model 
#' @param scoring_params (list) Additional arguments for forward model
#' @param mc.cores Number of parallel cores to use
#' @return A list containing model fits and classification results
#' @export
#' @examples
#' pars <- gen_params(n = 50, pY = 20, pX = 30, n_te = 100)
#' dat <- do.call(gen_data, pars)
#' zattach(dat)
#' "ML approach"
#' ml_res <- identification_pipeline1(X, Y, X_te, y_star, i_chosen,
#'   filter_method = no_filter, forward_method = fit_elnet_CV,
#'   Sigma_e_method = residual_offdiag, Sigma_t_method = assume_iid,
#'   backward_method = pre_mle, scoring_method = topk_score,
#'   mc.cores = 3)
#' ml_res$score
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
identification_pipeline1 <- function(X, Y, X_te, y_star, i_chosen,
                                    filter_method = no_filter,
                                    filter_params = list(),
                                    forward_method = fit_ridge_kernel,
                                    forward_params = list(),
                                    Sigma_e_method = residual_offdiag,
                                    Sigma_e_params = list(),
                                    Sigma_t_method = assume_iid,
                                    Sigma_t_params = list(),
                                    Sigma_b_method = assume_iid_B,
                                    Sigma_b_params = list(),
                                    backward_method = pre_mle,
                                    backward_params = list(),
                                    scoring_method = topk_score,
                                    scoring_params = list(),
                                    mc.cores = 1, ...) {
  res <- list(X = X, Y = Y, X_te = X_te, y_star = y_star,
                i_chosen = i_chosen, mc.cores = mc.cores)
  res$filt <- do.call(filter_method, modifyList(res, filter_params))
  res$B <- do.call(forward_method, modifyList(res, forward_params))
  res$Sigma_e <- do.call(Sigma_e_method, modifyList(res, Sigma_e_params))
  res$Sigma_t <- do.call(Sigma_t_method, modifyList(res, Sigma_t_params))
  res$Sigma_b <- do.call(Sigma_b_method, modifyList(res, Sigma_b_params))
  res$pre_moments <- do.call(backward_method, modifyList(res, backward_params))
  res$plikes <- do.call(post_likes, res)
  res$score <- do.call(scoring_method, modifyList(res, scoring_params))
  res
}