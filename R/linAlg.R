#' Inverse square root of a matrix
#' 
#' @param m matrix
#' @export
isqrtm <- function(m) {
  res <- eigen(m)
  d <- res$values
  if (min(d) < -1e-5) warning("Negative eigenvalues in isqrtm")
  d[d < 0] <- 0
  d[d > 0] <- 1/sqrt(d[d > 0])
  v <- res$vectors
  return (v %*% diag(d) %*% t(v))
}

#' Square root of a matrix
#' 
#' @param m matrix
#' @export
sqrtm <- function(m) {
  res <- eigen(m)
  d <- res$values
  if (min(d) < -1e-5) warning("Negative eigenvalues in sqrtm")
  d[d < 0] <- 0
  d[d > 0] <- sqrt(d[d > 0])
  v <- res$vectors
  return (v %*% diag(d) %*% t(v))
}

#' Frobenius norm squared
#' 
#' @param x vector
#' @param y (optional)
#' @export
f2 <- function(x, y = 0) sum((x - y)^2)

#' Computes \code{t(A \%x\% B) \%*\% diag(d) \%*\% (A \%*\% B)}
#' 
#' @import pracma
#' @import MASS
#' @import parallel
#' @param A matrix of size a x b
#' @param B matrix of size c x d
#' @param d vector: represents diagonal matrix of size a * c
#' @return square matrix of size a * c
#' @examples
#' a1 <- 5; a2 <- 3; b1 <- 4; b2 <- 6
#' A <- pracma::randn(a1, a2)
#' B <- pracma::randn(b1, b2)
#' cc <- rnorm(a1 * b1)
#' ans <- tkron_d_kron(A, B, cc)
#' ansn <- t(A %x% B) %*% diag(cc) %*% (A %x% B)
#' f2(ans, ansn)
#' @export
tkron_d_kron <- function(A, B, d) {
  a1 <- dim(A)[1]; a2 <- dim(A)[2]; b1 <- dim(B)[1]; b2 <- dim(B)[2]
  dmat <- matrix(d, b1, a1)
  # columns of dmat are diag(D1), diag(D2), ...
  C <- zeros(a2 * b2)
  if (a1 < b1) {
    bdbmat <- zeros(b2^2, a1)
    for (i in 1:a1) bdbmat[, i] <- as.numeric(t(B) %*% (dmat[, i] * B))
    for (i in 1:a2) {
      for (j in i:a2) {
        Cij <- matrix(bdbmat %*% (A[, i] * A[, j]), b2, b2)
        C[(i-1) * b2 + (1:b2), (j-1) * b2 + (1:b2)] <- Cij
        C[(j-1) * b2 + (1:b2), (i-1) * b2 + (1:b2)] <- t(Cij)
      }
    }
  } else {
    for (i in 1:a2) {
      for (j in i:a2) {
        dtemp <- as.numeric(dmat %*% (A[, i] * A[, j]))
        Cij <- t(B) %*% (dtemp * B)
        C[(i-1) * b2 + (1:b2), (j-1) * b2 + (1:b2)] <- Cij
        C[(j-1) * b2 + (1:b2), (i-1) * b2 + (1:b2)] <- t(Cij)
      }
    }  
  }
  C
}

#' Computes \code{(A \%x\% B) \%*\% v}
#' 
#' @param A matrix of size a x b
#' @param B matrix of size c x d
#' @param cc vector of length b * d
#' @return vector of size a * c
#' @examples
#' a1 <- 5; a2 <- 3; b1 <- 3; b2 <- 4
#' A <- pracma::randn(a1, a2)
#' B <- pracma::randn(b1, b2)
#' cc <- rnorm(a2 * b2)
#' ans <- kron_v(A, B, cc)
#' ansn <- (A %x% B) %*% cc
#' f2(ans, ansn)
#' @export
kron_v <- function(A, B, cc) {
  a1 <- dim(A)[1]; a2 <- dim(A)[2]; b1 <- dim(B)[1]; b2 <- dim(B)[2]
  as.numeric(B %*% matrix(cc, b2, a2) %*% t(A))
}

#' Splits a vector into approximately equally-sized parts
#' 
#' @param v Vector to be split
#' @param n.splits Number of splits
#' @return A list of vectors
#' @export
#' @examples
#' split_vec(1:100, 6)
split_vec <- function(v, n.splits) {
  if (n.splits == 0) n.splits <- 1
  ans <- list()
  n <- length(v)
  for (i in 1:n.splits) {
    li <- floor((i-1)/n.splits * n) + 1
    ui <- floor(i/n.splits * n)
    ans[[i]] <- v[li:ui]
  }
  ans
}

#' Multiple core matrix multiplication
#' 
#' Given A and B, computes AB by splitting columns of B
#' @param A matrix
#' @param B matrix with many columns
#' @return The product
#' @export
#' @examples
#' A <- randn(100, 100)
#' B <- randn(100, 1000)
#' C <- paramultiply(A, B, mc.cores = 3)
paramultiply <- function(A, B, mc.cores = 3) {
  rowpartition <- split_vec(1:dim(B)[2], mc.cores)
  res <- mclapply0(rowpartition, function(i) {
    A %*% B[, i]
  }, mc.cores = mc.cores)
  do.call(cbind, res)
}




