#' Classification Extrapolation with Regression (and bootstrap model selection)
#' 
#' @param pmat Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.
#' @param i_chosen True labels.
#' @param Ktrain the number of classes k for which we have computed average test accuracies ATA_k.
#' @param basis_sets a list of bases obtained e.g. using lapply and get_gaussian_basis_mat.
#' If there is only one element in the list, then model selection is not applied.
#' @param nboot Number of bootstrap samples to use in model selection.
#' @param nonnegative constrain coefficients of basis to be nonnegative? default TRUE.
ClassExReg <- function(pmat, i_chosen, basis_sets, Ktrain=2:ncol(pmat), nboot = 25, nonnegative = TRUE) {
  accs_sub <- sub_accuracies(pmat, i_chosen, Ktrain)
  if (length(basis_sets) > 0) {
    lsub2 <- floor(length(Ktrain)/2)
    sub_basis_sets <- lapply(basis_sets, function(set1) {
      list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(Ktrain), , drop = FALSE])
    })
    boot_accs <- matrix(NA, nboot, lsub2)
    for (ii in 1:nboot) {
      subinds <- sample(ksub, ksub2, replace = FALSE)
      counts_subsub <- fastRank(pmat[i_chosen %in% subinds, subinds], i_chosen[i_chosen %in% subinds])
      boot_accs[ii, ] <- count_acc(counts_subsub, kref[1:lsub2])
    }
  } else {
    sel_ind <- 1
  }

  
}