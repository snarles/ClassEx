#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

####
##  Usage: run as a shell script
##  Rscript mc_gaussian_sim_batch.R n1 n2 fdir
##  where n1 is the starting seed, n2 is the ending seed,
##  and fdir is the directory for saving the data
##  Use `collect_batch.R` to assemble the simulation output into a list,
##  and then `plot_sim.R` to view the results.
####

library(ClassEx)
#source("approximation/gaussian_identity_finsam.R")
#source("approximation/gaussian_identity_finsam2.R")
#source("extrapolation/ku_source.R")
#source("extrapolation/kay_method.R")
#source("par2/objective_function.R")
#source("extrapolation/basis_source.R")

####
##  SIMULATION PARAMETERS
####

p <- 10
sigma2_seq <- 0.01 * 1:50
K <- 10000 ## multiple of 1000
Ktarg <- c(1000, 2000, 5000, 10000)
ksub <- 500 ## multiple of 250
ksub2 <- ksub/2

nboot <- 20 ## number of bootstraps for CV
mc.reps <- 100000
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))

kseq <- function(nr, ksub) {
  interv <- floor(ksub/nr)
  seq(interv, ksub, by = interv)
}

nspline <- 100
nrow <- 100
kref <- kseq(nrow, ksub)
kde_bdwids <- list("bcv", "ucv", 0.1, 0.15, 0.2, 0.25, 0.3)
(max.mu <- (qnorm(1- 1/(max(kref)^2))))
## fixed bandwidths
fixed.gb <- c(0.4, 0.45, 0.5, 0.55, 0.6)
lsub2 <- length(kref)/2
basis_vecs <- list()

for (gb in fixed.gb) {
  nm <- paste0("r.gauss", gb)
  basis_vecs[[nm]] <- get_basis_mat(max.mu, kernel_sigma = gb)
}

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_basis_mat(max.mu, bd))

sub_basis_sets <- lapply(basis_sets, function(set1) {
  list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(kref), , drop = FALSE])
})

(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(names(basis_vecs), kde.names, "r.cv.gauss", "par2"))

repno <- 25
subfun <- function (repno) {
  set.seed(repno)
  sigma2 <- sigma2s[repno]
  sigma2_tr <- sigma2s[repno]
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  pmat_sub <- -pdist2(ys[1:ksub, ], mu_hats[1:ksub, ])
  rSqs <- rowSums((ys - mu_hats)^2)
  counts <- countDistEx(mu_hats, ys, rSqs)
  accs <- count_acc(counts, Ktarg)    
  counts_sub <- countDistEx(mu_hats[1:ksub,], ys[1:ksub,], rSqs[1:ksub])
  accs_sub <- count_acc(counts_sub, kref)
  boot_accs <- matrix(NA, nboot, lsub2)
  for (ii in 1:nboot) {
    subinds <- sample(ksub, ksub2, replace = FALSE)
    counts_subsub <- countDistEx(mu_hats[subinds,], ys[subinds,], rSqs[subinds])
    boot_accs[ii, ] <- count_acc(counts_subsub, kref[1:lsub2])
  }
  preds <- matrix(NA, length(column_names), length(Ktarg))
  rownames(preds) <- column_names
  for (ind in 1:length(column_names)) {
    if (ind <= length(basis_vecs)) {
      Xmat <- basis_vecs[[ind]]$Xmat
      Xpred <- basis_vecs[[ind]]$Xtarg
      bt <- nnls::nnls(Xmat, accs_sub)
      preds[ind, ] <- (Xpred %*% bt$x)
    } else if(ind <= length(basis_vecs) + length(kde_bdwids)) {
      bw <- kde_bdwids[[ind - length(basis_vecs)]]
      preds[ind, ] <- kernel_extrap(pmat_sub, Ktarg, bw = bw)
    } else if (column_names[ind] == "par2") {
      preds[ind, ] <- par2_extrapolate(kref, accs_sub, Ktarg)
    } else if (column_names[ind] == "r.cv.gauss") {
      all_sub_preds <- t(apply(boot_accs, 1, bdwid_all_preds, basis_sets = sub_basis_sets))
      cv_curve <- sqrt(colMeans((all_sub_preds - accs_sub[length(kref)])^2))
      sel_ind <- which.min(cv_curve)
      Xmat <- basis_sets[[sel_ind]]$Xmat
      Xpred <- basis_sets[[sel_ind]]$Xtarg
      bt <- nnls::nnls(Xmat, accs_sub)
      preds[ind, ] <- (Xpred %*% bt$x)
    }
  }
  list(preds = preds, accs = accs)
}  

res <- lapply(as.numeric(args[1]):as.numeric(args[2]), subfun)
save(res, file = paste0(args[3], "_", args[1], "_", args[2], ".RData"))

