#' kNN-based estimator (vectorised, fast) with tie/duplicate-aware neighbors
#'
#' A computationally efficient version of \code{kir_graph} that replaces the
#' inner loop over \eqn{j \neq i} with vectorised matrix operations. It
#' computes the same estimator
#' \deqn{\widehat D = 1 - \frac{1}{n}\sum_{i=1}^n
#'       \frac{\widehat E_i}{\max(\widehat V_i,\, \varepsilon)}}
#' as defined in Proposition 7, but is typically 10–100\eqn{\times} faster
#' than \code{kir_graph} for moderate-to-large \eqn{n} (e.g.\ \eqn{n = 1500}).
#'
#' @details
#' The neighbor matrix \eqn{N_{\cdot}} is computed once via \code{get_knn(X, Knn+1)}.
#' Inside the outer loop over \eqn{i}, point \eqn{i} is masked out of every
#' neighbor set \eqn{N_j} (\eqn{j \neq i}) by replacing its index with
#' \code{NA} and retaining the first \code{Knn} valid entries.  Because
#' \code{get_knn} returns \eqn{K+1} neighbors, at most one entry per row is
#' removed, so backfilling (random sampling from the remaining pool) is needed
#' only in the rare case that \eqn{i} appeared in \eqn{N_j}.
#'
#' Once the clean \eqn{(n-1) \times K} index matrix is available, the
#' double sum
#' \deqn{\sum_{j \neq i}\sum_{\ell \in N_j}(v_j - v_\ell)^2,
#'       \quad v = K_Y[\,,i]}
#' is evaluated in two vectorised steps: a gather \code{v[Nmat_clean]}
#' followed by \code{sum((v_j - v_neigh)^2)}, avoiding any inner R loop
#' over \eqn{j}.
#'
#' @param X Numeric matrix (or object coercible to one) of predictors with
#'   \eqn{n} rows (samples) and \eqn{d_X} columns (features).
#' @param Y Numeric matrix (or object coercible to one) of responses with
#'   \eqn{n} rows (samples) and \eqn{d_Y} columns.
#' @param Knn Positive integer. Number of nearest neighbors used to estimate
#'   the conditional variance structure. Internally clipped to
#'   \eqn{\min(K, n-1)}.
#' @param eps Small positive numeric scalar \eqn{\varepsilon > 0} used to
#'   stabilise the ratio \eqn{\widehat E_i / \widehat V_i} when
#'   \eqn{\widehat V_i \approx 0}. Defaults to \code{1e-12}.
#' @param ky Kernel on \eqn{Y}: either a \pkg{kernlab} kernel object or a
#'   function of the form \code{function(A, B = NULL)} returning the Gram
#'   matrix. If \code{NULL} (default), an RBF kernel with median-heuristic
#'   bandwidth scaled by \code{bandwidth_scale} is used.
#' @param bandwidth_scale Positive numeric scalar applied to the
#'   median-heuristic bandwidth when \code{ky = NULL}. Defaults to \code{1}.
#'
#' @return A single numeric scalar: the estimated dependence measure
#'   \eqn{\widehat D \in (-\infty, 1]}, where values close to 1 indicate
#'   strong dependence between \eqn{X} and \eqn{Y} and values close to 0
#'   indicate near independence.
#'
#' @seealso \code{\link{get_knn}} for the tie/duplicate-aware neighbor search.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' X <- matrix(rnorm(n * 2), n, 2)
#'
#' # Strong dependence: Y is a noisy function of X
#' Y_dep <- X[, 1, drop = FALSE] + 0.3 * rnorm(n)
#' kir_graph(X, Y_dep, Knn = 5)   # close to 1
#'
#' # Near independence
#' Y_ind <- matrix(rnorm(n), n, 1)
#' kir_graph(X, Y_ind, Knn = 5)   # close to 0
#'
#' # Custom kernel
#' if (requireNamespace("kernlab", quietly = TRUE)) {
#'   k <- kernlab::rbfdot(sigma = 0.5)
#'   kir_graph(X, Y_dep, Knn = 5, ky = k)
#' }
kir_graph <- function(X, Y, Knn = 5L, eps = 1e-12,
                           ky = NULL, bandwidth_scale = 1) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows.")
  if (n < 2L)       stop("Need at least 2 samples.")
  Knn <- min(as.integer(Knn), n - 1L)
  if (Knn < 1L)     stop("Knn must be >= 1.")
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0)
    stop("eps must be a single positive number.")

  # ---- kernel matrix ----
  if (is.null(ky)) ky <- .default_rbf_kernel(Y, bandwidth_scale = bandwidth_scale)
  Ky <- .compute_gram(ky, Y, Y)   # n x n

  # ---- neighbor matrix: n x (Knn+1), 1-based ----
  Knn_big  <- min(Knn + 1L, n - 1L)
  Nmat_big <- get_knn(X, Knn = Knn_big)   # n x Knn_big

  # ---- precompute Vhat for all i at once ----
  col_sum     <- colSums(Ky)
  col_ssq     <- colSums(Ky^2)
  diag_Ky     <- diag(Ky)
  mean_excl   <- (col_sum  - diag_Ky)  / (n - 1L)
  meansq_excl <- (col_ssq  - diag_Ky^2) / (n - 1L)
  Vhat_all    <- pmax(meansq_excl - mean_excl^2, eps)

  # ---- precompute affected rows ----
  # affected_rows[[i]] : rows j (!=i) whose Nmat_big row contains i
  # affected_col[[i]]  : which column of Nmat_big holds i in that row
  affected_rows <- vector("list", n)
  affected_col  <- vector("list", n)

  for (c in seq_len(Knn_big)) {
    col_vals <- Nmat_big[, c]
    for (j in seq_len(n)) {
      i_val <- col_vals[j]
      if (i_val >= 1L && i_val <= n && i_val != j) {
        affected_rows[[i_val]] <- c(affected_rows[[i_val]], j)
        affected_col[[i_val]]  <- c(affected_col[[i_val]],  c)
      }
    }
  }

  # ---- base neighbor matrix (first Knn cols only) ----
  Nmat_base <- Nmat_big[, seq_len(Knn), drop = FALSE]   # n x Knn

  Ehat_all <- numeric(n)

  for (i in seq_len(n)) {
    v  <- Ky[, i]
    js <- seq_len(n)[-i]

    # Copy base; patch only the ~Knn affected rows
    Nmat_i <- Nmat_base

    aff_j <- affected_rows[[i]]
    aff_c <- affected_col[[i]]

    if (length(aff_j) > 0L) {
      for (k in seq_along(aff_j)) {
        j0 <- aff_j[k]
        c0 <- aff_c[k]
        if (c0 <= Knn) {
          # i sits in the first Knn cols of row j0 -> replace with the lookahead col
          repl <- if (Knn_big > Knn) Nmat_big[j0, Knn + 1L] else NA_integer_
          if (!is.na(repl) && repl >= 1L && repl <= n && repl != i && repl != j0) {
            Nmat_i[j0, c0] <- repl
          } else {
            # rare backfill: sample one unused point
            pool <- setdiff(seq_len(n), c(i, j0, Nmat_i[j0, ]))
            Nmat_i[j0, c0] <- if (length(pool) > 0L) sample(pool, 1L) else repl
          }
        }
        # c0 > Knn: i was only in the lookahead slot -> base row already clean
      }
    }

    # Vectorised squared-difference sum over all j != i
    Nmat_js <- Nmat_i[js, , drop = FALSE]          # (n-1) x Knn
    v_j     <- v[js]                               # length n-1
    v_neigh <- matrix(v[Nmat_js], nrow = n - 1L, ncol = Knn)

    Ehat_all[i] <- sum((v_j - v_neigh)^2) / (2 * (n - 1L) * Knn)
  }

  1 - mean(Ehat_all / Vhat_all)
}
