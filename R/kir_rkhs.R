#' Kernel-integrated R2: RKHS estimator (KIR-RKHS)
#'
#' Computes a dependence/association score between \code{X} and \code{Y} using
#' an RBF kernel and an RKHS-based regression operator with Tikhonov regularization.
#'
#' @param X Numeric matrix/data frame of predictors (n x p).
#' @param Y Numeric matrix/data frame of responses (n x q).
#' @param bandwidth_scale Positive scalar multiplying the default RBF bandwidth.
#' @param eps Nonnegative ridge parameter (scaled by \code{n}) added to the centered
#'   Gram matrix diagonal for numerical stability.
#'
#' @returns A numeric scalar in principle near \eqn{[0, 1]}: larger values indicate
#'   stronger association.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' X <- matrix(rnorm(n * 2), n, 2)
#' Y <- X[, 1, drop = FALSE] + 0.2 * matrix(rnorm(n), n, 1)
#' kir_rkhs(X, Y)
kir_rkhs <- function(X, Y, bandwidth_scale = 1, eps = 1e-4) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  n <- nrow(X)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows (samples).")
  if (n < 2) stop("Need at least 2 samples.")
  if (!is.numeric(bandwidth_scale) || length(bandwidth_scale) != 1L || bandwidth_scale <= 0) {
    stop("bandwidth_scale must be a single positive number.")
  }
  if (!is.numeric(eps) || length(eps) != 1L || eps < 0) {
    stop("eps must be a single nonnegative number.")
  }

  # Median heuristic bandwidth from Y; guard against degenerate Y
  gamma <- stats::median(stats::dist(Y))
  if (!is.finite(gamma) || gamma <= 0) gamma <- 1

  sigma0 <- 1 / (2 * gamma^2)
  rbf <- kernlab::rbfdot(sigma = bandwidth_scale * sigma0)

  kX <- kernlab::kernelMatrix(rbf, X)
  kY <- kernlab::kernelMatrix(rbf, Y)

  # Center kX: H kX H (avoid constructing H)
  rmean <- rowMeans(kX)
  cmean <- colMeans(kX)
  gmean <- mean(kX)

  HkXH <- kX
  HkXH <- sweep(HkXH, 1, rmean, `-`)
  HkXH <- sweep(HkXH, 2, cmean, `-`)
  HkXH <- HkXH + gmean

  # Explicit inverse: inv(HkXH + n*eps*I)
  A <- HkXH
  diag(A) <- diag(A) + n * eps
  Ainv <- solve(A)

  B <- HkXH %*% Ainv
  M <- kY %*% B

  # Reuse terms efficiently
  kY1   <- rowSums(kY)
  kY2   <- kY * kY
  kY2_1 <- rowSums(kY2)

  M1    <- rowSums(M)
  diagMMt <- rowSums(M * M)

  # rowSums(kY^2 %*% B)
  kY2B_1 <- rowSums(kY2 %*% B)

  denom <- kY2_1 - (kY1 * kY1) / n
  nom <- kY2_1 + kY2B_1 -
    ((kY1 * kY1) / n + (2 / n) * kY1 * M1 + diagMMt)

  # Avoid division by ~0 in degenerate cases
  ratio <- nom / pmax(denom, .Machine$double.eps)
  1 - mean(ratio)
}

