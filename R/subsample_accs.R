#' Subroutine to compute average test accuracies ATA_k given correct-label ranks.
#' 
#' @param counts The ranks of the correct label for each tets instance.
#' @param K Total number of labels in test set.
#' @param subKs The number of classes k for which we want to calculate ATA_k.
count_acc <- function(ranks, K=length(counts), subKs = 1:length(counts)) {
  p_i <- 1 - ranks/K
  rowMeans(dhyper(zeros(length(subKs), length(p_i)), 
                  repmat(ranks, length(subKs), 1) - 1, 
                  K - repmat(ranks, length(subKs), 1), 
                  repmat(t(t(subKs)), 1, length(p_i)) - 
                    1))
}

#' Compute average test accuracies ATA_k.
#' 
#' @param pmat Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.
#' @param i_chosen True labels.
#' @param subKs The number of classes k for which we want to calculate ATA_k.
#' @export
sub_accuracies <- function(pmat, i_chosen, subKs = 1:ncol(pmat)) {
  count_acc(fastRank(pmat, i_chosen), ncol(pmat), subKs)
}

#' Compute average test accuracies ATA_k for 1-nearest neighbor.
#' @param Xm Matrix of centroids, rows are centroids and columns are dimensions.
#' @param Ym Data matrix, rows are observations and columns are dimensions.
#' @param i_chosen True labels.
#' @param subKs The number of classes k for which we want to calculate ATA_k.
#' @export
sub_accuracies_1nn <- function(Xm, Ym, i_chosen, subKs = 1:nrow(Xm)) {
  count_acc(fastRank1nn(Xm, Ym, i_chosen), nrow(Xm), subKs)
}