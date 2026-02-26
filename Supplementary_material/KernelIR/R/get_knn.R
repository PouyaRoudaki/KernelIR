#' K-Nearest Neighbors with Random Tie-Breaking
#'
#' Finds k-nearest neighbors for each point in a dataset, with special handling
#' for duplicate points and random sampling when multiple points are equidistant
#' (tied) at the k-th nearest neighbor distance.
#'
#' @param X A numeric matrix where each row represents a data point and each
#'   column represents a feature (coordinate). If not a matrix, it will be coerced to one.
#' @param Knn An integer specifying the number of nearest neighbors to find
#'   for each point.
#'
#' @return An integer matrix of dimension \code{n x Knn} where \code{n} is the
#'   number of rows in \code{X}. Each row \code{i} contains the indices
#'   (in 1-based R indexing) of the \code{Knn} nearest neighbors of point \code{i}.
#'
#' @details
#' This function combines the fast k-nearest neighbor search from the RANN package
#' with custom C++ code to handle two special cases:
#' \itemize{
#'   \item \strong{Duplicate points}: When multiple points have identical coordinates,
#'     they are grouped and neighbors are randomly sampled from the group.
#'   \item \strong{Tied distances}: When multiple points are equidistant at the k-th
#'     nearest neighbor position, the function randomly samples uniformly
#'     from the full tied set.
#' }
#'
#' The random sampling ensures unbiased neighbor selection and is useful for
#' statistical applications where deterministic tie-breaking could introduce bias.
#'
#' @note
#' This function requires the compiled C++ function
#' \code{knn_tie_duplicate_aware_cpp} to be available.
#'
#' @examples
#' ### 1) No duplicates, no ties (deterministic)
#' X <- matrix(c(0, 1, 3, 6, 10), ncol = 1)
#' get_knn(X, Knn = 3)
#'
#' ### 2) Symmetric tie case
#' X <- matrix(c(-2, -1, 0, 1, 2), ncol = 1)
#' set.seed(1)
#' get_knn(X, Knn = 3)
#' set.seed(2)
#' get_knn(X, Knn = 3)
#'
#' ### 3) Small duplicate group
#' X <- matrix(c(0, 0, 1, 2, 3), ncol = 1)
#' set.seed(1)
#' get_knn(X, Knn = 3)
#'
#' ### 4) Large duplicate block (> Knn)
#' X <- matrix(c(0, 0, 0, 0, 5, 6, 7), ncol = 1)
#' set.seed(1)
#' nn1 <- get_knn(X, Knn = 3)
#' set.seed(2)
#' nn2 <- get_knn(X, Knn = 3)
#' nn1
#' nn2
#'
#' ### 5) Large boundary tie set
#' X <- matrix(c(0, 1, 1, 1, 1, 5), ncol = 1)
#' set.seed(1)
#' get_knn(X, Knn = 3)
#'
#' ### 6) Evenly spaced grid (frequent ties)
#' X <- matrix(c(0, 2, 4, 6, 8, 10), ncol = 1)
#' set.seed(1)
#' get_knn(X, Knn = 3)
#'
#' ### 7) Small n edge case
#' X <- matrix(c(0, 1, 2, 3), ncol = 1)
#' get_knn(X, Knn = 3)
#'
#' ### 8) Duplicates + symmetric ties combined
#' X <- matrix(c(-1, 0, 0, 0, 1), ncol = 1)
#' set.seed(1)
#' get_knn(X, Knn = 3)
#'
#' @seealso
#' \code{\link[RANN]{nn2}}
#'
#' @useDynLib KernelIR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RANN nn2
#' @export
get_knn <- function(X, Knn) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  if (n < 2L) stop("Need at least 2 points.")
  Knn <- as.integer(Knn)
  Knn <- min(Knn, n - 1L)
  if (Knn < 1L) stop("Knn must be >= 1.")

  # k = Knn + 2 gives: self + Knn neighbors + 1 lookahead for tie detection
  k <- min(n, Knn + 2L)

  nn_result <- RANN::nn2(X, query = X, k = k)

  knn_tie_duplicate_aware_cpp(
    X = X,
    nn_idx = nn_result$nn.idx,
    nn_dists = nn_result$nn.dists,
    Knn = Knn,
    tol = 1e-12
  )
}
