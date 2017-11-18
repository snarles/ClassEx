#' Determines eigenprism weights
#' 
#' @param l Eigenvalues of X'X
#' @param tol Tolerance
#' @param sigma2 Compute weights for estimating sigma2 (TRUE) or signal strength (FALSE)?
#' @export
#' @references Janson 2015
eigenprism_w <- function(l, tol = 1e-4, sigma2 = FALSE) {
  n <- length(l)
  f <- function(w) max(sum(w^2), sum((w * l)^2))
  # null space directions
  res <- svd(cbind(1, l))
  P <- -(res$u %*% t(res$u))
  diag(P) <- diag(P) + 1
  A <- P[, -c(1, 2)]
  # Newton operators
  Ls <- list(diag(rep(1, n)), diag(l))
  Ns <- lapply(Ls, function(L) A %*% solve(t(A) %*% L %*% A) %*% t(A) %*% L)
  # initialization
  w <- 1/(l[1] - l[2]) * c(1, -1, rep(0, n - 2))
  if (sigma2) {
    w <- 1/(l[1] - l[2]) * c(-l[2], l[1], rep(0, n - 2))
  }
  f_old <- f(w)  
  iter_flag <- TRUE
  while(iter_flag) {
    # evaluate the gradient
    f_flag <- sum(w^2) < sum((w * l)^2)
    lterm <- l^(2 * f_flag)
    # constrained Newton direction
    g <- Ns[[1 + f_flag]] %*% w
    # t_minimize and t_eq
    t_min <- sum(w * g * lterm)/(sum(g^2 * lterm) + 1e-10)
    aa <- sum(g ^ 2) - sum((g * l)^2)
    bb <- -2 * (sum(g * w) - sum(l^2 * g * w))
    cc <- sum(w^2) - sum((w * l)^2)
    suppressWarnings({t_eq <- (-bb + c(-1, 1)*sqrt(bb^2 - 4 * aa * cc))/(2 * aa)})
    ts <- c(t_min, t_eq[!is.na(t_eq)])
    vals <- sapply(ts, function(t) f(w - t * g))
    t_best <- ts[order(vals)[1]]    
    w <- w - t_best * g
    f_new <- min(vals)
    if (f_old - f_new < tol) iter_flag <- FALSE
    f_old <- f_new
  }
  list(w = w[, 1], val = f_old)
}

#' Confidence intervals for ||beta||^2 or sigma2 (see Janson 2015)
#' 
#' @param X matrix of covariates
#' @param y vector response
#' @param alpha 1-coverage of confidence interval
#' @param sigma2 Compute weights for estimating sigma2 (TRUE) or signal strength (FALSE)?
#' @export
#' @references Janson 2015
eigenprism <- function(X, y, alpha = 0.05, sigma2 = FALSE) {
  n <- length(y)
  res <- svd(X)
  dim(res$u)
  z <- t(res$u) %*% y
  l <- res$d^2/p
  res <- eigenprism_w(l, sigma2 = sigma2)
  w <- res$w; val <- res$val
  T2 <- sum(z^2 * w)
  qq <- qnorm(1 - alpha/2)
  width <- qq * sqrt(2 * val) * sum(y^2)/n
  lower <- max(T2 - width, 0)
  upper <- T2 + width
  list(T2 = T2, lower = lower, upper = upper, w = w, l = l)
}

#' Confidence intervals for ||beta||^2 or sigma2 (see Janson 2015)
#' 
#' @param X matrix of covariates
#' @param Y matrix response
#' @param alpha 1-coverage of confidence interval
#' @param sigma2 Compute weights for estimating sigma2 (TRUE) or signal strength (FALSE)?
#' @export
#' @references Janson 2015
eigenprisms <- function(X, Y, alpha = 0.05, sigma2 = FALSE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  res <- svd(X)
  dim(res$u)
  z <- t(res$u) %*% Y
  l <- res$d^2/p
  res <- eigenprism_w(l, sigma2 = sigma2)
  w <- res$w; val <- res$val
  T2 <- apply(z^2 * w, 2, sum)
  qq <- qnorm(1 - alpha/2)
  width <- qq * sqrt(2 * val) * apply(Y^2, 2, sum)/n
  lower <- pmax(T2 - width, 0)
  upper <- T2 + width
  list(T2 = T2, lower = lower, upper = upper)
}