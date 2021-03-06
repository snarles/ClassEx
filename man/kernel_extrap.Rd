% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kde_method.R
\name{kernel_extrap}
\alias{kernel_extrap}
\title{KDE-based extrapolator for multi-class classification accuracy based on Kay (2008)}
\usage{
kernel_extrap(pmat, i_chosen, Ks, bw)
}
\arguments{
\item{pmat}{Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.}

\item{i_chosen}{The true labels for the observations.}

\item{Ks}{The number of classes k for which we want to extrpolate the accuracy AGA_k (can be a vector).}

\item{bw}{Bandwidth estimation method: see help("bw.nrd")}
}
\description{
KDE-based extrapolator for multi-class classification accuracy based on Kay (2008)
}
