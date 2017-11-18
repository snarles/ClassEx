

#' Computes posterior moments for Bayesian multivariate regression
#' 
#' Posterior mean and covariance given a prior on Bvec
#'  of the form \code{Sigma_b \%x\% eye(pX)}
#'  and noise prior of form \code{Sigma_e \%x\% Sigma_t}
#' @param X covariate matrix
#' @param Y response matrix
#' @param Sigma_e Covariance within given time (row) of Y
#' @param Sigma_b Prior covariance of a column of B
#' @param Sigma_t Autocorrelation of a single response (column) of Y over time
#' @param computeCov Return the posterior covariance of B?
#' @param naive Do the naive computation (for testing)
#' @param matrix Return B as a matrix? (or vector)
#' @return a matrix/vector (if \code{computeCov = FALSE})
#' or a list (if \code{computeCov = TRUE})
#' @export
#' @examples
#' n <- 100; pX <- 20; pY <- 30
#' X <- randn(n, pX)
#' B <- randn(pX, pY)
#' Sigma_b <- eye(pY)
#' Sigma_e <- cor(randn(3 * pY, pY))
#' Sigma_t <- cor(randn(3 * n, n))
#' E <- sqrtm(Sigma_t) %*% randn(n, pY) %*% sqrtm(Sigma_e)
#' Y <- X %*% B + E
#' res <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, TRUE, TRUE)
#' f2(solve(t(X) %*% X, t(X) %*% Y), B)
#' f2(res$Mu, B)
post_moments <- function(X, Y, Sigma_e, Sigma_b, Sigma_t = eye(dim(X)[1]), 
                         computeCov = TRUE, naive = FALSE, matrix = TRUE, ...) {
  n <- dim(X)[1]; pX <- dim(X)[2]; pY <- dim(Y)[2]
  Omega_e <- solve(Sigma_e)
  Omega_t <- solve(Sigma_t)
  Omega_b <- diag(1/diag(Sigma_b))
  yVec <- as.numeric(Y)
  xtx <- t(X) %*% Omega_t %*% X
  if (naive) {
    ans_mu <- solve(Omega_e %x% xtx + Omega_b %x% eye(pX),
                    t(eye(pY) %x% X) %*% (Omega_e %x% Omega_t) %*% yVec)
    if (matrix) ans_mu <- matrix(ans_mu, pX, pY)
    if (computeCov) {
      ans_cov <- solve(Omega_e %x% xtx + Omega_b %x% eye(pX))
      return(list(Mu = ans_mu, Cov = ans_cov))
    }
    return(ans_mu)
  }
  homega_b <- diag(1/sqrt(diag(Sigma_b)))
  hsigma_b <- diag(sqrt(diag(Sigma_b)))
  res <- eigen(hsigma_b %*% Omega_e %*% hsigma_b)
  D_e <- diag(res$values)
  V_e <- homega_b %*% res$vectors
  iV_e <- solve(V_e)
  resX <- eigen(xtx)
  V_x <- resX$vectors
  D_x <- diag(resX$values)
  d <- 1/(diag(D_e) %x% diag(D_x) + 1)
  temp <- kron_v(iV_e %*% Omega_e, t(V_x) %*% t(X) %*% Omega_t, yVec)
  ans_mu <- kron_v(t(iV_e), V_x, d * temp)
  if (matrix) ans_mu <- matrix(ans_mu, pX, pY)
  if (computeCov) {
    ans_cov <- tkron_d_kron(iV_e, t(V_x), d)
    return(list(Mu = ans_mu, Cov = ans_cov))
  }
  ans_mu
}


#' Posterior predictive means
#' 
#' Returns a list with entry for every row of Xte, which is a list with mu, cov
#' @param X covariate matrix
#' @param Y response matrix
#' @param X_te covariates for responses to be predicted
#' @param Sigma_e Covariance within given time (row) of Y
#' @param Sigma_b Prior covariance of a column of B
#' @param Sigma_t Autocorrelation of a single response (column) of Y over time
#' @param naive Do the naive computation (for testing)
#' @param mc.cores How many cores to use in parallel (0 for single-core)
#' @return a list of list, one sublist per row of X_te,
#' containing posterior mean and covariance
#' @export
#' @examples
#' n <- 100; pX <- 20; pY <- 30; L <- 10
#' X <- randn(n, pX)
#' X_te <- randn(L, pX)
#' B <- randn(pX, pY)
#' Sigma_b <- eye(pY)
#' Sigma_e <- cor(randn(3 * pY, pY))
#' Sigma_t <- cor(randn(3 * n, n))
#' E <- sqrtm(Sigma_t) %*% randn(n, pY) %*% sqrtm(Sigma_e)
#' Y <- X %*% B + E
#' E2 <- randn(L, pY) %*% sqrtm(Sigma_e)
#' Y_te <- X_te %*% B + E2
#' res <- post_predictive(X, Y, X_te, Sigma_e, Sigma_b, Sigma_t)
#' Bols <- solve(t(X) %*% X, t(X) %*% Y)
#' Yh <- do.call(rbind, listcomb(res)$Mu)
#' f2(X_te %*% Bols, Y_te)
#' f2(Yh, Y_te)
post_predictive <- function(X, Y, X_te, Sigma_e, Sigma_b, Sigma_t = eye(dim(X)[1]), 
                            naive = FALSE, mc.cores = 0, ...) {
  n <- dim(X)[1]; pX <- dim(X)[2]; pY <- dim(Y)[2]
  L <- dim(X_te)[1]; ans <- as.list(numeric(L))
  Omega_e <- solve(Sigma_e)
  Omega_t <- solve(Sigma_t)
  Omega_b <- diag(1/diag(Sigma_b))
  yVec <- as.numeric(Y)
  xtx <- t(X) %*% Omega_t %*% X
  
  if (naive) {
    res <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, naive = TRUE, matrix = TRUE)
    B <- res$Mu
    SigmaB <- res$Cov
    for (i in 1:L) {
      x_star <- X_te[i, ]
      Ix <- eye(pY) %x% t(x_star)
      Mu <- as.numeric(t(x_star) %*% B)
      Cov <- Ix %*% SigmaB %*% t(Ix) + Sigma_e
      ans[[i]] <- list(Mu = Mu, Cov = Cov)
    }
  } else {
    homega_b <- diag(1/sqrt(diag(Sigma_b)))
    hsigma_b <- diag(sqrt(diag(Sigma_b)))
    res <- eigen(hsigma_b %*% Omega_e %*% hsigma_b)
    D_e <- diag(res$values)
    V_e <- homega_b %*% res$vectors
    iV_e <- solve(V_e)
    resX <- eigen(xtx)
    V_x <- resX$vectors
    D_x <- diag(resX$values)
    d <- 1/(diag(D_e) %x% diag(D_x) + 1)
    temp <- d * kron_v(iV_e %*% Omega_e, t(V_x) %*% t(X) %*% Omega_t, yVec)
    pre_moments <- function(i) {
      x_star <- X_te[i, ]
      Mu <- kron_v(t(iV_e), t(x_star) %*% V_x, temp)
      Cov <- tkron_d_kron(iV_e, t(V_x) %*% x_star, d) + Sigma_e
      list(Mu = Mu, Cov = Cov)
    }
    ans <- mclapply0(1:L, pre_moments, mc.cores = mc.cores)
  }
  ans
}



#' Ridge regression subroutine
#' 
#' computes \code{solve(xtx + lambda, t(X) \%*\% Y)}
#' as \code{1/lambda * t(X) \%*\% Y}
#' \code{V_X, d_X} are from svd of \code{X}
#' @param X covariate/design matrix
#' @param V_X from svd of X
#' @param d_X from svd of X
#' @param lambda ridge regularization
#' @param Y response matrix
#' @param naive Do the naive computation (for testing)
#' @return a matrix B
#' @export
fast_ridge <- function(X, V_X, d_X, lambda, Y, naive = FALSE) {
  if (naive) {
    return(solve(t(X) %*% X + lambda * eye(dim(X)[2]), t(X) %*% Y))
  }
  if (!is.null(X)) {
    xty <- t(X) %*% Y    
  } else {
    xty <- Y
  }
  temp1 <- 1/lambda * xty
  d2 <- (1/(d_X^2 + lambda)) - 1/lambda
  temp2 <- V_X %*% (d2 * (t(V_X) %*% xty))
  temp1 + temp2
}

#' Sampling distribution of ridge regression
#'
#' @param X covariate matrix
#' @param Y response matrix
#' @param Sigma_e Covariance within given time (row) of Y
#' @param lambdas Ridge regularization per column of B
#' @param Sigma_t Autocorrelation of a single response (column) of Y over time
#' @param computeCov Return the posterior covariance of B?
#' @param matrix Return B as a matrix? (or vector)
#' @return a matrix/vector (if \code{computeCov = FALSE})
#' or a list (if \code{computeCov = TRUE})
#' @export
samp_moments <- function(X, Y, Sigma_e, lambdas, Sigma_t = eye(dim(X)[1]), 
                         computeCov = TRUE, matrix = TRUE, ...) {
  n <- dim(X)[1]; pX <- dim(X)[2]; pY <- dim(Y)[2]
  res_X <- svd(X); V_X <- res_X$v; d_X <- res_X$d
  xtx <- t(X) %*% X
  xtsx <- t(X) %*% Sigma_t %*% X
  #temp <- lapply(1:pY, function(i) solve(xtx + lambdas[i] * eye(pX), t(X) %*% Y[, i]))
  temp <- lapply(1:pY, function(i) fast_ridge(X, V_X, d_X, lambdas[i], Y[, i]))
  if (matrix) {
    ans_mu <- do.call(cbind, temp)
  } else {
    ans_mu <- do.call(c, temp)
  }
  if (computeCov) {
    ans_cov <- zeros(pX * pY, pX * pY)
    #d1s <- lapply(1:pY, function(i) solve(xtx + lambdas[i] * eye(pX), xtsx))
    #d2s <- lapply(1:pY, function(i) solve(xtx + lambdas[i] * eye(pX)))
    for (i in 1:pY) {
      for (j in i:pY) {
        #temp <- Sigma_e[i, j] * d1s[[i]] %*% d2s[[j]]
        temp0 <- fast_ridge(X, V_X, d_X, lambdas[j], Sigma_t)
        temp <- Sigma_e[i, j] * fast_ridge(X, V_X, d_X, lambdas[i], t(temp0))
        ans_cov[(i-1)*pX+(1:pX),(j-1)*pX+(1:pX)] <- temp
        ans_cov[(j-1)*pX+(1:pX),(i-1)*pX+(1:pX)] <- t(temp)
      }
    }
    return(list(Mu = ans_mu, Cov = ans_cov))
  }
  ans_mu
}

#' Frequentist predictive distributions for ridge regression
#' 
#' Returns a list with entry for every row of Xte, which is a list with mu, cov
#' @param X covariate matrix
#' @param Y response matrix
#' @param X_te covariates for responses to be predicted
#' @param Sigma_e Covariance within given time (row) of Y
#' @param lambdas Ridge regularization per column of B
#' @param Sigma_t Autocorrelation of a single response (column) of Y over time
#' @param mc.cores How many cores to use in parallel (0 for single-core)
#' @return a list of list, one sublist per row of X_te,
#' containing posterior mean and covariance
#' @export
samp_predictive <- function(X, Y, X_te, Sigma_e, lambdas, Sigma_t = eye(dim(X)[1]), 
                            mc.cores = 0, ...) {
  n <- dim(X)[1]; pX <- dim(X)[2]; pY <- dim(Y)[2]
  L <- dim(X_te)[1]; ans <- as.list(numeric(L))
  res <- samp_moments(X, Y, Sigma_e, lambdas, Sigma_t, computeCov = TRUE, matrix = TRUE)
  BMu <- res$Mu
  BCov <- res$Cov
  pre_moments <- function(i) {
    x_star <- X_te[i, ]
    Mu <- as.numeric(t(x_star) %*% BMu)
    temp <- zeros(pY, pY * pX)
    for (ii in 1:pY) {
      temp[ii, ] <- t(x_star) %*% BCov[(ii - 1)*pX + (1:pX), ]
    }
    Cov <- temp %*% (eye(pY) %x% t(t(x_star))) + Sigma_e
    list(Mu = Mu, Cov = Cov)
  }
  ans <- mclapply0(1:L, pre_moments, mc.cores = mc.cores)
  ans
}
