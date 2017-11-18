
#' Similar to attach, but overwrites global
#'
#' @param ll list with named values
#' @export
#' @examples
#' zattach(list(a = 5))
#' a
zattach <- function(ll) {
  for (i in 1:length(ll)) {
    assign(names(ll)[i], ll[[i]], envir=globalenv())
  }
}

#' Similar to do.call, but allows extra arguments
#' 
#' @param f function to be evaluated
#' @param ll list containing arguments
#' @export
#' @examples
#' do.call2(cbind, list(a = rep(0, 3), b = rep(1, 3), c = rep(2, 3)), a = c(1, 2, 3), d = c(4, 5, 4))
do.call2 <- function(f, ll, ...) {
  dots <- eval(substitute(alist(...)))
  ll <- modifyList(ll, dots)
  do.call(f, ll)
}

#' Converts list(list(fields)) into list(fields(list))
#' 
#' @param dots list of lists
#' @export
#' @examples
#' listcomb(list(list(a = 5, b = 3), list(a = 2, b = 4)))
listcomb <- function(dots) {
  nms <- names(dots[[1]])
  ans <- as.list(nms)
  names(ans) <- nms
  for (nm in nms) {
    ans[[nm]] <- lapply(1:length(dots), function(i) dots[[i]][[nm]])
  }
  ans
}

#' Calls mclapply and combines the end result using listcomb
#' 
#' @param x List of arguments
#' @param f Single-argument function
#' @param mc.cores Number of cores: 0 means use lapply
#' @import parallel
#' @export
#' @examples
#' lclapply(1:10, function(i) list(a = i, b = i^2), mc.cores = 1)
lclapply <- function(x, f, mc.cores = 0) {
  if (mc.cores == 0) {
    return(listcomb(lapply(x, f)))
  } else {
    return(listcomb(mclapply(x, f, mc.cores = mc.cores)))    
  }
}

#' Either uses lapply or mclapply
#' 
#' @param x List of arguments
#' @param f Single-argument function
#' @param mc.cores Number of cores: 0 means use lapply
#' @import parallel
#' @export
#' @examples
#' lclapply(1:10, function(i) list(a = i, b = i^2), mc.cores = 1)
mclapply0 <- function(x, f, mc.cores = 0) {
  if (mc.cores == 0) {
    return(lapply(x, f))
  } else {
    return(mclapply(x, f, mc.cores = mc.cores))    
  }
}
