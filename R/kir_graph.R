#' kNN-based estimator with tie/duplicate-aware neighbors and exclusion of \eqn{i} from \eqn{N_j}
#'
#' Computes
#' \deqn{\widehat D = 1 - \frac{1}{n}\sum_{i=1}^n \frac{\widehat E_i}{\max(\widehat V_i, \varepsilon)}}
#' where \eqn{\widehat E_i} and \eqn{\widehat V_i} correspond to the quantities in Proposition 7.
#' The neighbor sets \eqn{N_j} are obtained via \code{get_knn()}, which is assumed to be
#' tie/duplicate-aware.
#'
#' For each \eqn{i}, the computation uses \eqn{v = K_Y[, i]} where \eqn{K_Y} is an RBF kernel
#' matrix on \eqn{Y}. For each \eqn{j \neq i}, we form \eqn{N_j} from \code{get_knn(X, Knn = Knn + 1)}
#' and then remove \eqn{i} from that neighbor list (to enforce \eqn{i \notin N_j}). If fewer than
#' \code{Knn} neighbors remain, the set is backfilled by sampling without replacement from the
#' remaining indices excluding \eqn{\{i, j\}} and the already-chosen neighbors.
#'
#' If \code{debias = TRUE}, the squared mean term \eqn{(\frac{1}{|N_j|}\sum_{\ell\in N_j} v_\ell)^2}
#' is debiased as \eqn{\bar v^2 - s^2/|N_j|}, where \eqn{s^2} is the sample variance over \eqn{N_j}.
#'
#' @param X Numeric matrix/data.frame of predictors with \code{n} rows (samples).
#' @param Y Numeric matrix/data.frame of responses with \code{n} rows (samples).
#' @param Knn Positive integer number of neighbors (per row). Internally clipped to \code{min(Knn, n-1)}.
#' @param eps Small positive constant \eqn{\varepsilon} used to stabilize the ratio when
#'   \eqn{\widehat V_i} is close to zero.
#' @param ky Kernel on Y (kernlab kernel object or function(A,B=NULL)->Gram). If NULL, defaults to RBF.
#' @param bandwidth_scale Positive scalar for default RBF bandwidth (used only when ky is NULL).
#'
#' @return A numeric scalar: the estimate \eqn{\widehat D}.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(rnorm(n * 2), n, 2)
#' Y <- X[, 1, drop = FALSE] + 0.3 * rnorm(n)
#' kir_graph(X, Y, Knn = 5)
kir_graph <- function(X,
                      Y,
                      Knn = 5L,
                      eps = 1e-12,
                      ky = NULL,
                      bandwidth_scale = 1) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  n <- nrow(X)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows (samples).")
  if (n < 2L) stop("Need at least 2 samples.")
  Knn <- as.integer(Knn)
  Knn <- min(Knn, n - 1L)
  if (Knn < 1L) stop("Knn must be >= 1.")
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0) {
    stop("eps must be a single positive number.")
  }

  # Default (old behavior): RBF on Y if not supplied
  if (is.null(ky)) ky <- .default_rbf_kernel(Y, bandwidth_scale = bandwidth_scale)

  Ky <- .compute_gram(ky, Y, Y)  # n x n

  Knn_big  <- min(Knn + 1L, n - 1L)
  Nmat_big <- get_knn(X, Knn = Knn_big)

  idx_all <- seq_len(n)

  ratios <- vapply(idx_all, function(i) {
    idx_minus_i <- idx_all[idx_all != i]
    v <- Ky[, i]

    sum_jk <- sum(vapply(idx_minus_i, function(j) {
      Nj <- Nmat_big[j, ]

      Nj <- Nj[Nj != i]
      if (length(Nj) >= Knn) {
        Nj <- Nj[seq_len(Knn)]
      } else {
        pool <- setdiff(idx_all, c(i, j, Nj))
        if (length(pool) > 0L) {
          Nj <- c(Nj, sample(pool, Knn - length(Nj), replace = FALSE))
        }
      }

      if (length(Nj) == 0L) 0 else sum((v[j] - v[Nj])^2)
    }, numeric(1)))

    Ehat_i <- sum_jk / (2 * (n - 1L) * Knn)

    term1  <- mean(v[idx_minus_i]^2)
    mean_v <- mean(v[idx_minus_i])
    Vhat_i <- term1 - mean_v^2

    Ehat_i / max(Vhat_i, eps)
  }, numeric(1))

  1 - mean(ratios)
}
