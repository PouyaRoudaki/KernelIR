# Compute an n x m Gram matrix from:
#  - a kernlab kernel object (rbfdot, polydot, ...)
#  - OR a function(A, B = NULL) returning a Gram matrix
.compute_gram <- function(kernel, A, B = NULL) {
  A <- as.matrix(A)
  if (!is.null(B)) B <- as.matrix(B)

  if (is.null(kernel)) stop("kernel is NULL.")

  if (inherits(kernel, "kernel")) {
    if (is.null(B)) {
      return(kernlab::kernelMatrix(kernel, A))
    } else {
      return(kernlab::kernelMatrix(kernel, A, B))
    }
  }

  if (is.function(kernel)) {
    K <- kernel(A, B)
    if (!is.matrix(K)) K <- as.matrix(K)
    return(K)
  }

  stop("kernel must be a kernlab kernel object or a function(A, B = NULL) returning a Gram matrix.")
}

# Default RBF kernel via median heuristic on data Z
.default_rbf_kernel <- function(Z, bandwidth_scale = 1) {
  if (!is.numeric(bandwidth_scale) || length(bandwidth_scale) != 1L || bandwidth_scale <= 0) {
    stop("bandwidth_scale must be a single positive number.")
  }
  Z <- as.matrix(Z)
  gamma <- stats::median(stats::dist(Z))
  if (!is.finite(gamma) || gamma <= 0) gamma <- 1
  sigma0 <- 1 / (2 * gamma^2)
  kernlab::rbfdot(sigma = bandwidth_scale * sigma0)
}

# Center a square Gram matrix: H K H without forming H
.center_gram <- function(K) {
  rmean <- rowMeans(K)
  cmean <- colMeans(K)
  gmean <- mean(K)

  Kc <- K
  Kc <- sweep(Kc, 1, rmean, `-`)
  Kc <- sweep(Kc, 2, cmean, `-`)
  Kc <- Kc + gmean
  Kc
}
