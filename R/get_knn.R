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
#'     nearest neighbor position, the function randomly samples which points to include
#'     as neighbors.
#' }
#'
#' The random sampling ensures unbiased neighbor selection and is useful for
#' statistical applications where deterministic tie-breaking could introduce bias.
#'
#' @note
#' This function requires the compiled C++ function \code{handle_duplicates_and_ties_cpp}
#' to be available. Use \code{Rcpp::sourceCpp()} to compile the associated C++ code
#' before calling this function.
#'
#' @examples
#' \dontrun{
#' # Create sample data with some duplicate points
#' set.seed(123)
#' X <- matrix(rnorm(1000), ncol = 5)
#' X [ 10:15, ] <- X[ 1, ]  # Add duplicates
#'
#' # Find 5 nearest neighbors with random tie-breaking
#' nn_indices <- get_knn(X, Knn = 5)
#'
#' # The results will differ between runs when ties exist
#' nn_indices2 <- get_knn(X, Knn = 5)
#' sum(nn_indices != nn_indices2)  # Shows differences due to random tie-breaking
#' }
#'
#' @seealso
#' \code{\link[RANN]{nn2}} for the underlying nearest neighbor search,
#' \code{\link[base]{sample}} for the random sampling behavior
#'
#' @useDynLib IKMD, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RANN nn2
#' @export
get_knn <- function(X, Knn) {
  if (!is.matrix(X)) X <- as.matrix(X)

  # Use RANN for initial neighbor search
  nn_result <- RANN::nn2(X, query = X, k = Knn + 2)
  # Note that we have we need at least Knn +1 because the first element is the
  # the node itself and in order to flag ties at the boundary case we need to
  # check whether distance of Knn+1 and Knn +2 are same or not. If they are same
  # then we need to find all the nodes that have same distance and take a sample
  # of them.

  # Use Rcpp for handling duplicates and ties
  nn_indices <- handle_duplicates_and_ties(
    X,
    nn_result$nn.idx,
    nn_result$nn.dists,
    Knn
  )

  return(nn_indices)
}
