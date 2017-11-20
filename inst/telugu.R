library(class)
library(pracma)
library(ClassEx)
library(parallel)

load("inst/telugu.rda", verbose = TRUE)
set.seed(0)

accZ20 <- lapply(lp20, function(v) sub_accuracies(v, rep(1:20, each = 50), 1:20))
accZ100 <- lapply(lp100, function(v) sub_accuracies(v, rep(1:100, each = 50), 1:100))

matplot(1:20, do.call(cbind, accZ20), type = "l", ylim = c(0,1))
matplot(1:100, do.call(cbind, accZ100), type = "l", ylim = c(0,1))

####
## predict 100 -> 400
####

K <- 400
ksub <- 100
ksub2 <- ksub/2
nboot <- 20 ## number of bootstraps for CV
i_chosen <- rep(1:100, each = 50)

kde_bdwids <- list("bcv", "ucv")
(max.mu <- (qnorm(1- 1/(ksub^2))))
lsub2 <- ksub2

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_gaussian_basis_mat(max.mu, bd, 1:ksub, K))


(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(kde.names, "r.cv.gauss"))

subfun <- function(pmat) {
  res <- c(kernel_extrap(pmat, i_chosen, K, bw = "ucv"),
           kernel_extrap(pmat, i_chosen, K, bw = "bcv"),
           ClassExReg(pmat, i_chosen, basis_sets, 1:ksub, nboot))
  names(res) <- column_names
  res
}

pmat <- lp100[[1]]
Ktrain <- 1:ksub


(preds <- lapply(lp100, subfun))
(resids <- lapply(methods, function(m) preds[[m]] - accs400[m]))
do.call(rbind, resids)

(preds100 <- cbind(do.call(rbind, preds), accs400))

####
## predict 20 -> 400
####

K <- 400
ksub <- 20
ksub2 <- ksub/2
nboot <- 20 ## number of bootstraps for CV
i_chosen <- rep(1:20, each = 50)

kde_bdwids <- list("bcv", "ucv")
(max.mu <- (qnorm(1- 1/(ksub^2))))
lsub2 <- ksub2

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_gaussian_basis_mat(max.mu, bd, 1:ksub, K))


(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(kde.names, "r.cv.gauss"))

(preds <- lapply(lp20, subfun))
(resids <- lapply(methods, function(m) preds[[m]] - accs400[m]))
do.call(rbind, resids)

(preds20 <- cbind(do.call(rbind, preds), accs400))

####
## predict 20 -> 100
####

K <- 100
ksub <- 20
ksub2 <- ksub/2
nboot <- 20 ## number of bootstraps for CV
i_chosen <- rep(1:20, each = 50)

kde_bdwids <- list("bcv", "ucv")
(max.mu <- (qnorm(1- 1/(ksub^2))))
lsub2 <- ksub2

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_gaussian_basis_mat(max.mu, bd, 1:ksub, K))


(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(kde.names, "r.cv.gauss"))

(preds <- lapply(lp20, subfun))
(resids <- lapply(methods, function(m) preds[[m]] - accs100[m]))
do.call(rbind, resids)

(preds20_100 <- cbind(do.call(rbind, preds), accs100))