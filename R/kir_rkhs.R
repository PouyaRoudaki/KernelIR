#' Kernel-integrated R2: RKHS estimator (KIR-RKHS)
#'
#' Computes a dependence/association score between \code{X} and \code{Y} using
#' an RBF kernel and an RKHS-based regression operator with Tikhonov regularization.
#'
#' @param X Numeric matrix/data frame of predictors (n x p).
#' @param Y Numeric matrix/data frame of responses (n x q).
#' @param kx Kernel for X (kernlab kernel object or function(A,B=NULL)->Gram). If NULL, defaults to RBF.
#' @param ky Kernel for Y (kernlab kernel object or function(A,B=NULL)->Gram). If NULL, defaults to RBF.
#' @param bandwidth_scale Positive scalar multiplying the default RBF bandwidth (used only when kernel is NULL).
#' @param eps Nonnegative ridge parameter (scaled by n) for stability.
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
kir_rkhs <- function(X,
                     Y,
                     kx = NULL,
                     ky = NULL,
                     bandwidth_scale = 1,
                     eps = 1e-4) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  n <- nrow(X)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows (samples).")
  if (n < 2) stop("Need at least 2 samples.")
  if (!is.numeric(eps) || length(eps) != 1L || eps < 0) {
    stop("eps must be a single nonnegative number.")
  }

  # Defaults (old behavior): RBF for each domain if not supplied
  if (is.null(kx)) kx <- .default_rbf_kernel(X, bandwidth_scale = bandwidth_scale)
  if (is.null(ky)) ky <- .default_rbf_kernel(Y, bandwidth_scale = bandwidth_scale)

  kX <- .compute_gram(kx, X, X)  # n x n
  kY <- .compute_gram(ky, Y, Y)  # n x n

  # Center kX: H kX H
  HkXH <- .center_gram(kX)

  # Ainv = (HkXH + n*eps*I)^{-1}
  A <- HkXH
  diag(A) <- diag(A) + n * eps
  Ainv <- solve(A)

  B <- HkXH %*% Ainv
  M <- kY %*% B

  # Reuse terms efficiently
  kY1   <- rowSums(kY)
  kY2   <- kY * kY
  kY2_1 <- rowSums(kY2)

  M1      <- rowSums(M)
  diagMMt <- rowSums(M * M)

  kY2B_1 <- rowSums(kY2 %*% B)

  denom <- kY2_1 - (kY1 * kY1) / n
  nom <- kY2_1 + kY2B_1 -
    ((kY1 * kY1) / n + (2 / n) * kY1 * M1 + diagMMt)

  ratio <- nom / pmax(denom, .Machine$double.eps)
  1 - mean(ratio)
}

