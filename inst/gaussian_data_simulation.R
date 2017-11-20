library(ClassEx)
library(pracma)
library(parallel)

####
##  SIMULATION PARAMETERS
####

p <- 10 ## dimension of Gaussian vectors
sigma2_seq <- 0.01 * 1:50 ## noise levels which determine SNR range of simulation
K <- 10000 ## number of classes to extrapolate to
ksub <- 500 ## number of classes in training set
nboot <- 20 ## number of bootstraps for model selection
mc.reps <- 500 ## total number of simulation instances (should be a multiple of the length of sigma2_seq)
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))
Ktrain <- 1:25 * 20 ## classes to subsample
kde_bdwids <- list("bcv", "ucv")
(max.mu <- (qnorm(1- 1/(ksub^2))))
(bdwids <- seq(0.1, 1, 0.1)) ## bandwidths for model selection
## get the ClassExReg Basis Sets for varying guassian bandwidths
#basis_sets <- lapply(bdwids, function(bd) get_gaussian_basis_mat(max.mu, bd, Ktrain, K))
#saveRDS(basis_sets, file = "inst/basis_sets.rds")
basis_sets <- readRDS("inst/basis_sets.rds")

## inspect basis matrices
for (i in 1:length(basis_sets)) {
  matplot(basis_sets[[i]]$Xmat, type = "l", main = paste("Xmat bdwid=",bdwids[i]))
  plot(basis_sets[[i]]$Xtarg, type = "l", main = paste("Xtarg bdwid=", bdwids[i]))
}

####
##  Simulation code
####

subfun <- function(repno) {
  ## Generate Data
  set.seed(repno)
  sigma2 <- sigma2s[repno]
  sigma2_tr <- sigma2s[repno]
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  pmat_sub <- -pdist2(ys[1:ksub, ], mu_hats[1:ksub, ])
  accs <- sub_accuracies_1nn(mu_hats, ys, 1:K, K)
  Xm <- mu_hats[1:ksub, ]
  Ym <- ys[1:ksub, ]
  i_chosen <- 1:ksub
  nonnegative <- TRUE
  preds <- c(true = accs,
           KDE_UCV = kernel_extrap(pmat_sub, 1:ksub, K, bw = "ucv"),
           KDE_BCV = kernel_extrap(pmat_sub, 1:ksub, K, bw = "bcv"),
           ClassExReg = ClassExReg1nn(mu_hats[1:ksub, ], ys[1:ksub, ], 1:ksub, basis_sets, Ktrain, nboot))
  preds
}

####
##  Run simulation
####

res <- mclapply(1:mc.reps, subfun, mc.cores = 8)
res <- do.call(rbind, res)

####
##  Analysis of Simulation results
####

true_accs <- numeric(length(sigma2_seq))
for (ii in 1:length(sigma2_seq)) {
  true_accs[ii] <- mean(res[sigma2s == sigma2_seq[ii], 1])
}
plot(sigma2_seq, true_accs, type = "l", ylim = c(0, 1))

resids <- res[, -1] - true_accs[match(sigma2s, sigma2_seq)]
rmses <- matrix(NA, length(sigma2_seq), ncol(resids))
biases <- matrix(NA, length(sigma2_seq), ncol(resids))
rmse_sd <- matrix(NA, length(sigma2_seq), ncol(resids))
bias_sd <- matrix(NA, length(sigma2_seq), ncol(resids))
nboot <- 100
for (ii in 1:length(sigma2_seq)) {
  rmses[ii, ] <- sqrt(colMeans(resids[sigma2s == sigma2_seq[ii], ]^2))
  biases[ii, ] <- colMeans(resids[sigma2s == sigma2_seq[ii], ])
  rmse_boots <- matrix(NA, nboot, ncol(rmses))
  bias_boots <- matrix(NA, nboot, ncol(rmses))
  orig_inds <- which(sigma2s == sigma2_seq[ii])
  for (jj in 1:nboot) {
    boot_inds <- sample(orig_inds, length(orig_inds), replace = TRUE)
    rmse_boots[jj, ] <- sqrt(colMeans(resids[boot_inds, ]^2))
    bias_boots[jj,] <- colMeans(resids[boot_inds, ])
  }
  rmse_sd[ii, ] <- apply(rmse_boots, 2, sd)
  bias_sd[ii,] <- apply(bias_boots,2,sd)
}

column_names <- colnames(resids)
colnames(rmses) <- column_names
colnames(rmse_sd) <- column_names
colnames(bias_sd) <- column_names
colnames(biases) <- column_names

####
##  Plot results
####

## RMSES

library(reshape2)
library(ggplot2)
sd_mult <- 2.95
temp <- data.frame(true_acc = true_accs, rmses)
temp2 <- melt(data = temp, id.vars = "true_acc")
temp_se <- data.frame(true_acc = true_accs, rmse_sd)
temp_se <- melt(data = temp_se, id.vars = "true_acc")
temp3 <- data.frame(temp2, rmse_low = temp2[, 3] - temp_se[, 3], rmse_high = temp2[, 3] + temp_se[, 3])
colnames(temp3)[3] <- "rmse"
ggplot(data = temp3, aes(x = true_acc, y = rmse, colour = variable, linetype=variable)) +
  geom_line() + coord_cartesian(xlim = c(0, 1)) + 
  geom_errorbar(aes(ymin = rmse_low, max = rmse_high)) + 
  scale_linetype_manual(values = c(1,4,6))+
  ggtitle(paste0("Predicting K=", K, " from k=", ksub)) + ylim(c(0,0.11))

## Biases

temp <- data.frame(true_acc = true_accs, biases)
temp2 <- melt(data = temp, id.vars = "true_acc")
temp_se <- data.frame(true_acc = true_accs, bias_sd)
temp_se <- melt(data = temp_se, id.vars = "true_acc")
temp3 <- data.frame(temp2, rmse_low = temp2[, 3] - temp_se[, 3], rmse_high = temp2[, 3] + temp_se[, 3])
colnames(temp3)[3] <- "bias"
ggplot(data = temp3, aes(x = true_acc, y = bias, colour = variable, linetype=variable)) +
  geom_line() + coord_cartesian(xlim = c(0, 1)) + 
  geom_errorbar(aes(ymin = rmse_low, max = rmse_high)) + 
  scale_linetype_manual(values = c(1,4,6))+
  ggtitle(paste0("Predicting K=", K, " from k=", ksub)) + ylim(c(-0.11,0.11))
