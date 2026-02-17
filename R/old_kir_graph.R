#' kNN-based estimator using tie/duplicate-aware kNN
#'
#' Computes the estimator
#' \deqn{\widehat D = 1 - \frac{1}{n}\sum_{i=1}^n \frac{\widehat E_i}{\widehat V_i}}
#' where \eqn{\widehat E_i} and \eqn{\widehat V_i} are defined in Proposition 7 and
#' the kNN sets \eqn{N_j} are obtained from \code{get_knn()}, which is assumed to be
#' tie/duplicate-aware.
#'
#' @param X Numeric matrix/data.frame of predictors with \code{n} rows (samples).
#' @param Y Numeric matrix/data.frame of responses with \code{n} rows (samples).
#' @param Knn Positive integer number of neighbors. Will be clipped to \code{min(Knn, n-1)}.
#' @param eps Small positive constant used to stabilize the division when \eqn{\widehat V_i}
#'   is close to zero.
#'
#' @returns A numeric scalar: the estimated value \eqn{\widehat D}.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(rnorm(n * 2), n, 2)
#' Y <- X[, 1, drop = FALSE] + 0.3 * rnorm(n)
#' kir_graph(X, Y, Knn = 5)
old_kir_graph <- function(X, Y, Knn = 5L, eps = 1e-12) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  n <- nrow(X)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows (samples).")
  if (n < 2L) stop("Need at least 2 samples.")
  if (!is.numeric(Knn) || length(Knn) != 1L || is.na(Knn) || Knn < 1) {
    stop("Knn must be a positive integer.")
  }
  Knn <- as.integer(Knn)
  Knn <- min(Knn, n - 1L)

  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0) {
    stop("eps must be a single positive number.")
  }

  # Default ky: Gaussian RBF with median heuristic on Y
  gamma <- stats::median(stats::dist(Y))
  if (!is.finite(gamma) || gamma <= 0) gamma <- 1
  sigma <- 1 / (2 * gamma^2)

  rbf <- kernlab::rbfdot(sigma = sigma)
  Ky  <- kernlab::kernelMatrix(rbf, Y, Y)  # Ky[a, b] = ky(Y_a, Y_b)

  # Neighbor indices: n x Knn (1-based), tie/duplicate-aware
  Nmat <- get_knn(X, Knn = Knn)
  if (!is.matrix(Nmat) || nrow(Nmat) != n || ncol(Nmat) != Knn) {
    stop("get_knn() must return a matrix of dimension n x Knn.")
  }

  idx_all <- seq_len(n)

  ratios <- vapply(idx_all, function(i) {
    idx <- idx_all[idx_all != i]

    # v[j] = ky(Y_j, Y_i)
    v <- Ky[, i]

    # term1 = (1/(n-1)) * sum_{j != i} ky(Y_j, Y_i)^2
    term1 <- mean(v[idx]^2)

    # term2 = (1/(n-1)) * sum_{j != i} [ avg_{s in N_j \ {i}} ky(Y_s, Y_i) ]^2
    avg_sq <- vapply(idx, function(j) {
      Nj <- Nmat[j, ]
      Nj <- Nj[Nj != i]   # N_j \ {i}
      if (length(Nj) == 0L) 0 else (mean(v[Nj]))^2
    }, numeric(1))

    term2 <- mean(avg_sq)
    Ehat_i <- term1 - term2

    # Vhat_i
    mean_v <- mean(v[idx])
    Vhat_i <- term1 - mean_v^2

    Ehat_i / max(Vhat_i, eps)
  }, numeric(1))

  1 - mean(ratios)
}
