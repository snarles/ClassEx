#' Classification Extrapolation with Regression (and bootstrap model selection)
#' 
#' @import nnls
#' @param pmat Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.
#' @param i_chosen True labels.
#' @param Ktrain the number of classes k for which we have computed average test accuracies ATA_k.
#' @param basis_sets a list of bases obtained e.g. using lapply and get_gaussian_basis_mat.
#' If there is only one element in the list, then model selection is not applied.
#' @param nboot Number of bootstrap samples to use in model selection.
#' @param nonnegative constrain coefficients of basis to be nonnegative? default TRUE.
#' @export
ClassExReg <- function(pmat, i_chosen, basis_sets, Ktrain=2:ncol(pmat), nboot = 25, nonnegative = TRUE) {
  accs_sub <- sub_accuracies(pmat, i_chosen, Ktrain)
  ksub <- max(Ktrain)
  if (length(basis_sets) > 1) {
    lsub2 <- floor(length(Ktrain)/2)
    sub_basis_sets <- lapply(basis_sets, function(set1) {
      list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(Ktrain), , drop = FALSE])
    })
    boot_accs <- matrix(NA, nboot, lsub2)
    for (ii in 1:nboot) {
      subinds <- sort(sample(ksub, Ktrain[lsub2], replace = FALSE))
      i_chosen2 <- match(i_chosen[i_chosen %in% subinds], subinds)
      boot_accs[ii, ] <- sub_accuracies(pmat[i_chosen %in% subinds, subinds], i_chosen2, Ktrain[1:lsub2])
    }
    all_sub_preds <- t(apply(boot_accs, 1, bdwid_all_preds, basis_sets = sub_basis_sets))
    cv_curve <- sqrt(colMeans((all_sub_preds - accs_sub[length(Ktrain)])^2))
    sel_ind <- which.min(cv_curve)
  } else {
    sel_ind <- 1
  }
  Xmat <- basis_sets[[sel_ind]]$Xmat
  Xpred <- basis_sets[[sel_ind]]$Xtarg
  if (nonnegative) {
    bt <- nnls(Xmat, accs_sub)$x
  } else {
    bt <- solve(t(Xmat) %*% Xmat, t(Xmat) %*% accs_sub)
  }
  return(Xpred %*% bt)
}



#' Classification Extrapolation with Regression (and bootstrap model selection) for 1 nearest neighbor
#' 
#' @import nnls
#' @param Xm Matrix of centroids, rows are centroids and columns are dimensions.
#' @param Ym Data matrix, rows are observations and columns are dimensions.
#' @param i_chosen True labels.
#' @param Ktrain the number of classes k for which we have computed average test accuracies ATA_k.
#' @param basis_sets a list of bases obtained e.g. using lapply and get_gaussian_basis_mat.
#' If there is only one element in the list, then model selection is not applied.
#' @param nboot Number of bootstrap samples to use in model selection.
#' @param nonnegative constrain coefficients of basis to be nonnegative? default TRUE.
#' @export
ClassExReg1nn <- function(Xm, Ym, i_chosen, basis_sets, Ktrain=2:nrow(Xm), nboot = 25, nonnegative = TRUE) {
  accs_sub <- sub_accuracies_1nn(Xm, Ym, i_chosen, Ktrain)
  ksub <- max(Ktrain)
  if (length(basis_sets) > 1) {
    lsub2 <- floor(length(Ktrain)/2)
    sub_basis_sets <- lapply(basis_sets, function(set1) {
      list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(Ktrain), , drop = FALSE])
    })
    boot_accs <- matrix(NA, nboot, lsub2)
    for (ii in 1:nboot) {
      subinds <- sort(sample(ksub, Ktrain[lsub2], replace = FALSE))
      i_chosen2 <- match(i_chosen[i_chosen %in% subinds], subinds)
      boot_accs[ii, ] <- sub_accuracies_1nn(Xm[subinds, ], Ym[i_chosen %in% subinds, ], 
                                            i_chosen2, Ktrain[1:lsub2])
    }
    all_sub_preds <- t(apply(boot_accs, 1, bdwid_all_preds, basis_sets = sub_basis_sets))
    cv_curve <- sqrt(colMeans((all_sub_preds - accs_sub[length(Ktrain)])^2))
    sel_ind <- which.min(cv_curve)
  } else {
    sel_ind <- 1
  }
  Xmat <- basis_sets[[sel_ind]]$Xmat
  Xpred <- basis_sets[[sel_ind]]$Xtarg
  if (nonnegative) {
    bt <- nnls(Xmat, accs_sub)$x
  } else {
    bt <- solve(t(Xmat) %*% Xmat, t(Xmat) %*% accs_sub)
  }
  return(Xpred %*% bt)
}