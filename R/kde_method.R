
#' Subroutine to get the KDE estimate of Pr[X > x] using Gaussian Kernels
#' 
#' @param xs Observations X_i
#' @param x Query point x
#' @param bw Bandwidth estimation method: see help("bw.nrd")
gaussian_kernel_cdf <- function(xs, x, bw = "bcv") {
  dens <- density(xs, bw=bw)
  ps <- pnorm(x, mean = xs, sd = dens$bw)
  mean(ps)
}

#' Subroutine to get all of the KDE estimates of the favorability for each test instance
#' 
#' @param pmat Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.
#' @param i_chosen The true labels for the observations.
#' @param bw Bandwidth estimation method: see help("bw.nrd")
raccs <- function(pmat, i_chosen, bw) {
  sapply(1:nrow(pmat), function(ind) gaussian_kernel_cdf(pmat[ind, -i_chosen[ind]], pmat[ind, i_chosen[ind]], bw))
}

#' KDE-based extrapolator for multi-class classification accuracy based on Kay (2008)
#' 
#' @param pmat Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.
#' @param i_chosen The true labels for the observations.
#' @param Ks The number of classes k for which we want to extrpolate the accuracy AGA_k (can be a vector).
#' @param bw Bandwidth estimation method: see help("bw.nrd")
#' @export
kernel_extrap <- function(pmat, i_chosen, Ks, bw) {
  racs <- raccs(pmat, i_chosen, bw)
  sapply(Ks, function(k) mean(racs^(k-1)))
}