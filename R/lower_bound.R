
#' Lower bound for ABA as a function of information
#' 
#' @param k Number of classes
#' @param aba average Bayes accuracy
#' @export
aba_to_mi_lower <- function(k, aba) {
  if (aba < 1/k) return(0)
  if (aba == 1) return(Inf)
  abas <- iota_aba_table[iota_aba_table$k==k, "aba"]
  a1 <- max(abas[abas < aba])
  a2 <- min(abas[abas > aba])
  i1 <- iota_aba_table$iota[iota_aba_table$k==k & iota_aba_table$aba==a1]
  i2 <- iota_aba_table$iota[iota_aba_table$k==k & iota_aba_table$aba==a2]
  (aba - a1)/(a2 - a1) * i2 + (a2 - aba)/(a2 - a1) * i1
}
