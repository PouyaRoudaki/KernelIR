# =========================
# 0) Setup (your libraries)
# =========================
library(KPC)      # KPCRKHS, KMAc
library(KernelIR) # kir_rkhs, kir_graph

library(doParallel)
library(foreach)

# =========================
# 1) Params
# =========================
n <- 100
n_sim <- 500
B <- 1000                 # permutations per dataset

#n <- 10
#n_sim <- 10
#B <- 10                 # permutations per dataset

alpha <- 0.05

lambda_vals <- seq(0, 1, 0.1)
alternatives <- c("SO3")

set.seed(7)

# parallel setup: respect Slurm allocation
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

# =========================
# 2) Data generator
# =========================

# Some auxiliary functions
wrap_to_pi <- function(a) atan2(sin(a), cos(a))

R1 <- function(x) {
  matrix(c(
    1, 0, 0,
    0, cos(x), -sin(x),
    0, sin(x),  cos(x)
  ), nrow = 3, byrow = TRUE)
}

R3 <- function(z) {
  matrix(c(
    cos(z), -sin(z), 0,
    sin(z),  cos(z), 0,
    0,       0,      1
  ), nrow = 3, byrow = TRUE)
}

vec_mat <- function(A) as.numeric(A)  # column-major vectorization (R default)

generate_data_so3_mv <- function(lambda, seed, n = 100) {
  set.seed(seed)

  X <- matrix(rnorm(n * 3), ncol = 3)
  E <- matrix(rnorm(n * 2), ncol = 2)

  # paper model:
  # Y = R1(X1) R3(X2 X3)

  a <- wrap_to_pi(X[,1] + lambda * E[,1])
  b <- wrap_to_pi(X[,2] * X[,3] + lambda * E[,2])

  Y <- matrix(NA_real_, nrow = n, ncol = 9)

  for (i in seq_len(n)) {
    A <- R1(a[i]) %*% R3(b[i])
    Y[i, ] <- as.numeric(A)
  }

  list(X = X, Y = Y)
}


# =========================
# 3) Generic permutation p-value (right tail)
# =========================
perm_pval_righttail <- function(obs, null) {
  (1 + sum(null >= obs)) / (length(null) + 1)
}

# ========================
# 3.5) Specify the kernel for simulation
# ========================

library(kernlab)

wrap_to_pi <- function(a) atan2(sin(a), cos(a))

so3_theta_from_vec <- function(y1, y2) {
  # y1,y2 are length-9 vectors encoding 3x3 rotation matrices in column-major order
  A <- matrix(y1, nrow = 3, ncol = 3)
  B <- matrix(y2, nrow = 3, ncol = 3)

  R <- t(B) %*% A
  cth <- (sum(diag(R)) - 1) / 2
  cth <- max(-1, min(1, cth))
  acos(cth)
}

so3_char_kernel_scalar <- function(theta, eps = 1e-12) {
  s <- sin(theta)
  if (abs(s) < eps) return(pi^2 / 8)  # correct limit at 0 and pi
  (pi * theta * (pi - theta)) / (8 * s)
}

# ---- kernlab kernel object for SO(3)^m (m blocks, each 9 numbers) ----

so3_kernlab <- function(m = 1L, combine = c("product", "mean"), eps = 1e-12) {
  combine <- match.arg(combine)

  kfun <- function(x, y) {
    # x,y are numeric vectors of length 9*m
    stopifnot(length(x) == 9L * m, length(y) == 9L * m)

    vals <- numeric(m)
    for (b in seq_len(m)) {
      idx <- ((b - 1L) * 9L + 1L):(b * 9L)
      th  <- so3_theta_from_vec(x[idx], y[idx])
      vals[b] <- so3_char_kernel_scalar(th, eps = eps)
    }

    if (combine == "product") prod(vals) else mean(vals)
  }

  class(kfun) <- "kernel"
  kfun
}




# =========================
# 4) One place to define ALL test statistics
#    Each returns a single numeric.
# =========================
stats_list <- list(

  kpcrkhs = function(X, Y) {
    ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")

    # In KPCRKHS you used Z = X; so kxz should be kernel on Z (same as X here)
    KPC::KPCRKHS(Y = Y, Z = X, ky = ky_so3)
  },

  kmagraph_5nn  = function(X, Y) {
    ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")
    KPC::KMAc(Y = Y, X = X, Knn = 5,  k = ky_so3)
  },

  kir_rkhs = function(X, Y) {
    ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")
    KernelIR::kir_rkhs(X = X, Y = Y, ky = ky_so3, kx = NULL)
  },

  kir_graph_5nn = function(X, Y) {
    ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")
    KernelIR::kir_graph(X = X, Y = Y, Knn = 5,  ky = ky_so3)
  }
)

stat_names <- names(stats_list)

# =========================
# 5) Core engine: rejections for ONE dataset using ONE permutation bank
# =========================
one_dataset_rejections <- function(x, y, B, alpha, seed_perm = NULL) {
  if (!is.null(seed_perm)) set.seed(seed_perm)

  n_y <- if (is.matrix(y)) nrow(y) else length(y)

  # permutation indices: n x B
  perm_idx <- replicate(B, sample.int(n_y))

  # observed stats
  obs <- vapply(stats_list, function(f) f(x, y), numeric(1))

  # null stats: B x M
  null_mat <- matrix(NA_real_, nrow = B, ncol = length(stats_list))
  colnames(null_mat) <- stat_names

  for (b in seq_len(B)) {
    yb <- if (is.matrix(y)) y[perm_idx[, b], , drop = FALSE] else y[perm_idx[, b]]
    null_mat[b, ] <- vapply(stats_list, function(f) f(x, yb), numeric(1))
  }

  # right-tail permutation p-values
  pvals <- vapply(seq_along(obs),
                  function(j) perm_pval_righttail(obs[j], null_mat[, j]),
                  numeric(1))
  names(pvals) <- stat_names

  as.numeric(pvals < alpha)
}

# =========================
# 6) Simulation wrapper: returns power estimates (colMeans)
# =========================
compute_power <- function(type, lambda) {

  results <- foreach(
    i = 1:n_sim, .combine = rbind,
    .packages = c("KPC","KernelIR","kernlab"),
    .export = c(
      "generate_data_so3_mv",
      "wrap_to_pi","R1","R3","vec_mat",
      "perm_pval_righttail","one_dataset_rejections",
      "n","B","alpha",
      "so3_kernlab","so3_theta_from_vec","so3_char_kernel_scalar"
    )
  ) %dopar% {

    # define stats_list on the worker (so it definitely sees so3_kernlab)
    stats_list_local <- list(
      kpcrkhs = function(X, Y) {
        ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")
        KPC::KPCRKHS(Y = Y, Z = X, ky = ky_so3, eps = 1e-4)
      },

      kmagraph_5nn  = function(X, Y) {
        ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")
        KPC::KMAc(Y = Y, X = X, Knn = 5,  k = ky_so3)
      },

      kir_rkhs = function(X, Y) {
        ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")
        KernelIR::kir_rkhs(X = X, Y = Y, ky = ky_so3, kx = NULL, eps = 1e-4)
      },

      kir_graph_5nn = function(X, Y) {
        ky_so3 <- so3_kernlab(m = ncol(Y)/9, combine = "product")
        KernelIR::kir_graph(X = X, Y = Y, Knn = 5,  ky = ky_so3)
      }
    )
    stat_names_local <- names(stats_list_local)

    # data
    dat <- generate_data_so3_mv(lambda, seed = i,n)
    x <- dat$X
    y <- dat$Y

    # permutation test computed locally (so it uses local stats_list)
    if (!is.null(100000 + i)) set.seed(100000 + i)
    n_y <- nrow(y)
    perm_idx <- replicate(B, sample.int(n_y))

    obs <- vapply(stats_list_local, function(f) f(x, y), numeric(1))

    null_mat <- matrix(NA_real_, nrow = B, ncol = length(stats_list_local))
    colnames(null_mat) <- stat_names_local
    for (b in seq_len(B)) {
      yb <- y[perm_idx[, b], , drop = FALSE]
      null_mat[b, ] <- vapply(stats_list_local, function(f) f(x, yb), numeric(1))
    }

    pvals <- vapply(seq_along(obs),
                    function(j) perm_pval_righttail(obs[j], null_mat[, j]),
                    numeric(1))
    rej <- as.numeric(pvals < alpha)
    names(rej) <- stat_names_local
    rej
  }

  colMeans(results)
}

# =========================
# 7) Main loop
# =========================
results_all <- list()


for (type in alternatives) {
  powers_all <- vector("list", length(lambda_vals))
  names(powers_all) <- as.character(lambda_vals)

  for (lambda in lambda_vals) {
    cat(sprintf("Running %s | lambda = %.1f\n", type, lambda))
    powers_all[[as.character(lambda)]] <- compute_power(type, lambda)
  }

  powers_df <- do.call(rbind, powers_all)
  powers_df <- data.frame(
    alternative = type,
    lambda = lambda_vals,
    powers_df,
    row.names = NULL
  )

  saveRDS(powers_df, file = paste0("power_", type, ".rds"))
  results_all[[type]] <- powers_df
}

parallel::stopCluster(cl)





