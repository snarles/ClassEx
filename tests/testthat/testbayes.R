library(lineId)
context("Testing Bayesian Regression")

# complex Frobenius norm
cf2 <- function(x, y = 0) Re(sum(Conj(x - y) * (x - y)))

# Normal cf
charfn <- function(tt, mu, Sigma) exp(1i * sum(mu * tt) - (t(tt) %*% Sigma %*% tt)/2)


# Tests Bayesian method
#
# Let (B, Y) have a joint distribution, where Y is observed
# Then for any functions f(B), g(Y) we have
# E[f(B)g(Y)] = E[h(Y)g(Y)]
# where h(Y) = E[f(B)|Y]


test_that("Posterior distribution is correct", {
  ntrial <- 1000
  n <- 100; pX <- 20; pY <- 30
  X <- randn(n, pX)
  Sigma_b <- diag(rexp(pY))
  Sigma_e <- cor(randn(3 * pY, pY))
  Sigma_t <- cor(randn(3 * n, n))
  sst <- sqrtm(Sigma_t)
  sse <- sqrtm(Sigma_e)
  Sigma_B <- (eye(pX) %x% Sigma_b)
  Sigma_Y <- (X %*% t(X)) %x% Sigma_b + Sigma_t %x% Sigma_e
  Sigma_BY <- rbind(cbind(Sigma_B, t(X) %x% Sigma_b),
                    cbind(X %x% Sigma_b, Sigma_Y))
  BYs <- lclapply(1:ntrial, function(i) {
    B <- t(sqrt(diag(Sigma_b)) * randn(pY, pX))
    E <- sst %*% randn(n, pY) %*% sse
    list(B = B, Y = X %*% B + E)
  }, mc.cores = 3)
  fmu <- rnorm(pX * pY)
  nfmu <- (t(fmu) %*% Sigma_B %*% fmu)[1]
  fmu <- fmu/sqrt(nfmu) * runif(1)
  gmu <- rnorm(n * pY)
  ngmu <- (t(gmu) %*% Sigma_Y %*% gmu)[1]
  gmu <- gmu/sqrt(ngmu) * runif(1)
  fBs <- sapply(BYs$B, function(B) exp(1i * sum(as.numeric(B) * fmu)))
  gYs <- sapply(BYs$Y, function(Y) exp(1i * sum(as.numeric(Y) * gmu)))
  fBgYs <- fBs * gYs
  fBs0 <- charfn(fmu, rep(0, pX * pY), Sigma_B)[1]
  gYs0 <- charfn(gmu, rep(0, n * pY), Sigma_Y)[1]
  fBgYs0 <- charfn(c(fmu, gmu), rep(0, (n + pX) * pY), Sigma_BY)[1]
  ## c(empirical value, theoretical value)
  c(mean(fBs), fBs0)
  c(mean(gYs), gYs0)
  c(mean(fBgYs), fBgYs0)
  ## computation using h
  Y <- BYs$Y[[1]]
  hY <- function(Y) {
    res <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, TRUE)
    charfn(fmu, as.numeric(res$Mu), res$Cov)[1]
  }
  hYs <- unlist(mclapply(BYs$Y, hY, mc.cores = 3))
  c(mean(fBs), mean(hYs), fBs0)
  c(mean(fBgYs), mean(hYs * gYs), fBgYs0)
  expect_less_than(cf2(mean(hYs * gYs), fBgYs0), 1e-2)  
})



test_that("Predictive distribution is correct", {
  ntrial <- 1000
  n <- 100; pX <- 20; pY <- 30; L <- 1
  X <- randn(n, pX)
  Xte <- randn(L, pX)
  Sigma_b <- diag(rexp(pY))
  Sigma_e <- cor(randn(3 * pY, pY))
  Sigma_t <- cor(randn(3 * n, n))
  sst <- sqrtm(Sigma_t)
  sse <- sqrtm(Sigma_e)
  Sigma_B <- (eye(pX) %x% Sigma_b)
  Sigma_Y <- (X %*% t(X)) %x% Sigma_b + Sigma_t %x% Sigma_e
  #Sigma_BY <- rbind(cbind(Sigma_B, t(X) %x% Sigma_b),
  #                  cbind(X %x% Sigma_b, Sigma_Y))
  Sigma_Yte <- (Xte %*% t(Xte)) %x% Sigma_b + eye(L) %x% Sigma_e
  Sigma_tt <- rbind(cbind(Sigma_t, zeros(n, L)),
                    cbind(zeros(L, n), eye(L)))
  XX <- rbind(X, Xte)
  Sigma_YY <- (XX %*% t(XX)) %x% Sigma_b + Sigma_tt %x% Sigma_e
  
  BYs <- lclapply(1:ntrial, function(i) {
    B <- t(sqrt(diag(Sigma_b)) * randn(pY, pX))
    E <- sst %*% randn(n, pY) %*% sse
    E_te <- randn(L, pY) %*% sse
    list(B = B, Y = X %*% B + E, Yte = Xte %*% B + E_te)
  }, mc.cores = 3)
  
  fmu <- rnorm(L * pY)
  nfmu <- (t(fmu) %*% Sigma_Yte %*% fmu)[1]
  fmu <- fmu/sqrt(nfmu) * runif(1)
  
  gmu <- rnorm(n * pY)
  ngmu <- (t(gmu) %*% Sigma_Y %*% gmu)[1]
  gmu <- gmu/sqrt(ngmu) * runif(1)
  
  fBs <- sapply(BYs$Yte, function(B) exp(1i * sum(as.numeric(B) * fmu)))
  gYs <- sapply(BYs$Y, function(Y) exp(1i * sum(as.numeric(Y) * gmu)))
  fBgYs <- fBs * gYs
  
  fBs0 <- charfn(fmu, rep(0, L * pY), Sigma_Yte)[1]
  gYs0 <- charfn(gmu, rep(0, n * pY), Sigma_Y)[1]
  fBgYs0 <- charfn(c(fmu, gmu), rep(0, (n + L) * pY), Sigma_YY)[1]
  
  ## c(empirical value, theoretical value)
  c(mean(fBs), fBs0)
  c(mean(gYs), gYs0)
  c(mean(fBgYs), fBgYs0)
  
  ## computation using h
  hY <- function(Y) {
    res <- post_predictive(X, Y, Xte, Sigma_e, Sigma_b, Sigma_t)[[1]]
    charfn(fmu, as.numeric(res$Mu), res$Cov)[1]
  }
  hYs <- unlist(mclapply(BYs$Y, hY, mc.cores = 3))
  c(mean(fBs), mean(hYs), fBs0)
  c(mean(fBgYs), mean(hYs * gYs), fBgYs0)
  expect_less_than(cf2(mean(hYs * gYs), fBgYs0), 1e-2)
})

## TODO: test sampling distributions