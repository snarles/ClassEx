% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassExReg.R
\name{ClassExReg}
\alias{ClassExReg}
\title{Classification Extrapolation with Regression (and bootstrap model selection)}
\usage{
ClassExReg(pmat, i_chosen, basis_sets, Ktrain = 2:ncol(pmat), nboot = 25,
  nonnegative = TRUE)
}
\arguments{
\item{pmat}{Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.}

\item{i_chosen}{True labels.}

\item{basis_sets}{a list of bases obtained e.g. using lapply and get_gaussian_basis_mat.
If there is only one element in the list, then model selection is not applied.}

\item{Ktrain}{the number of classes k for which we have computed average test accuracies ATA_k.}

\item{nboot}{Number of bootstrap samples to use in model selection.}

\item{nonnegative}{constrain coefficients of basis to be nonnegative? default TRUE.}
}
\description{
Classification Extrapolation with Regression (and bootstrap model selection)
}
