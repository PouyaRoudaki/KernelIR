# =============================================================================
# Experiment 7.5: Convergence Rate Adapts to Intrinsic Dimension (Theorem 2)
#
# Theorem 2:
#   |D̂^{K-NN}(Y,X) - D(Y,X)| = O_P( 1/sqrt(n) + (K/n)^{β/d} + (log n)^{β/α}/n^2 )
#
# Three experiments (each panel in the figure):
#
#   Exp 1 — Fix d=5, vary D ∈ {5,10,20,50}
#     All RMSE curves should share the same slope ~−β/d.
#     Message: ambient dim D does NOT affect convergence rate.
#
#   Exp 2 — Fix D=20, vary d ∈ {1,2,5,10,20}
#     RMSE curves show clearly different slopes −min(β/d, 1/2) (β=1).
#     Message: intrinsic dim d DOES drive convergence rate.
#
#   Exp 3 — Fix d=5, D=20, vary K ∈ {3,5,10,20}
#     Parallel curves (same slope), shifted by K^{β/d}.
#     Message: K shifts the constant, not the rate.
#
# Data model — both X and Y are D-dimensional with intrinsic dimension d:
#   Z_x      ~ N(0, I_d)                         (latent X signal)
#   Z_y      = Z_x + 0.5 * N(0, I_d)            (noisy copy of Z_x)
#   Q_x, Q_y ~ Haar(O(D))                        (independent random orthogonal)
#   X        = Z_x %*% t(Q_x[,1:d])  ∈ R^{n×D}  (on d-dim subspace)
#   Y        = Z_y %*% t(Q_y[,1:d])  ∈ R^{n×D}  (on d-dim subspace)
#
# Since KNN distances on a d-dimensional subspace of R^D equal distances in R^d
# (Q is an isometry), the effective KNN dimension is d, not D.
#
# D(Y,X) depends on d (more active dimensions → richer dependence structure), so
# we compute a separate oracle for each distinct d value used.  The oracle is
# evaluated at D=d (minimal ambient = intrinsic) to maximise accuracy at n=10000;
# D(Y,X) is the same for D=d and any D>d because the pairwise distances are
# preserved by the orthogonal embedding.
# =============================================================================

library(KernelIR)
library(doParallel)
library(foreach)

# =============================================================================
# 1) Parameters
# =============================================================================
d_fixed_exp1 <- 5L
D_vary_vals  <- c(5L, 10L, 20L, 50L)   # Exp 1: ambient dim varies, d fixed

D_fixed_exp2 <- 20L
d_vary_vals  <- c(1L, 2L, 5L, 10L, 20L) # Exp 2: intrinsic dim varies, D fixed

d_fixed_exp3 <- 5L
D_fixed_exp3 <- 20L
K_vary_vals  <- c(3L, 5L, 10L, 20L)    # Exp 3: K varies, d and D fixed

K_fixed      <- 5L
n_vals       <- c(50, 100, 200, 500, 1000, 2000, 5000)
#n_vals       <- c(50, 100, 200, 500, 1000)  # smaller n_vals for quick testing
n_sim        <- 200
#n_sim        <- 10
n_oracle_n   <- 10000L
#n_oracle_n   <- 200L
n_oracle_sim <- 30L
#n_oracle_sim <- 1L

out_dir <- "results_convergence_rate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(7)

n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

# =============================================================================
# 2) Data generator — both X and Y are D-dimensional with intrinsic dim d
# =============================================================================
random_orthogonal <- function(D) {
  A      <- matrix(rnorm(D * D), D, D)
  qr_res <- qr(A)
  Q      <- qr.Q(qr_res)
  Q * rep(sign(diag(qr.R(qr_res))), each = D)
}

generate_data <- function(n, D, d, seed) {
  set.seed(seed)
  Z_x <- matrix(rnorm(n * d), n, d)
  Z_y <- Z_x + 0.5 * matrix(rnorm(n * d), n, d)   # noisy copy of Z_x in R^d
  Q_x <- random_orthogonal(D)
  Q_y <- random_orthogonal(D)
  X   <- Z_x %*% t(Q_x[, 1:d, drop = FALSE])       # n × D, rank d
  Y   <- Z_y %*% t(Q_y[, 1:d, drop = FALSE])       # n × D, rank d
  list(X = X, Y = Y)
}

# =============================================================================
# 3) Oracles — one per distinct d value
#
#    D(Y,X) depends on d but not on D, so we evaluate at D = d (cheapest) for
#    maximum accuracy at n = 10000.
# =============================================================================
all_d_oracle <- sort(unique(c(d_fixed_exp1, d_vary_vals, d_fixed_exp3)))
cat(sprintf("Computing oracles (D = d, n = %d, %d reps each)...\n",
            n_oracle_n, n_oracle_sim))

oracle_by_d <- sapply(all_d_oracle, function(d) {
  ests <- foreach(
    i         = seq_len(n_oracle_sim),
    .combine  = c,
    .packages = "KernelIR",
    .export   = c("generate_data", "random_orthogonal", "n_oracle_n", "K_fixed")
  ) %dopar% {
    dat <- generate_data(n_oracle_n, D = as.integer(d), d = as.integer(d),
                         seed = 90000L + i)
    KernelIR::kir_graph(dat$X, dat$Y, Knn = K_fixed)
  }
  cat(sprintf("  D*(d=%d) = %.4f (sd = %.4f)\n", d, mean(ests), sd(ests)))
  mean(ests)
})
names(oracle_by_d) <- as.character(all_d_oracle)

# =============================================================================
# 4) Core parallel worker
# =============================================================================
run_cell <- function(d, D, Knn, n) {
  foreach(
    i         = seq_len(n_sim),
    .combine  = c,
    .packages = "KernelIR",
    .export   = c("generate_data", "random_orthogonal", "n_sim")
  ) %dopar% {
    dat <- generate_data(n, D, d, seed = i)
    KernelIR::kir_graph(dat$X, dat$Y, Knn = as.integer(Knn))
  }
}

# Returns an n_vals × 3 matrix: columns rmse, ci_lo, ci_hi (95% CI).
# CI is derived via the delta method on the mean squared error:
#   SE(MSE) = sd(sq_err) / sqrt(n_sim)
#   CI for RMSE = sqrt( MSE ± 1.96 * SE(MSE) )
rmse_grid <- function(d_val, D_val, Knn_val, D_star) {
  t(vapply(n_vals, function(n) {
    cat(sprintf("    n = %4d ...", n))
    ests   <- run_cell(d = d_val, D = D_val, Knn = Knn_val, n = n)
    sq_err <- (ests - D_star)^2
    mse    <- mean(sq_err)
    se_mse <- sd(sq_err) / sqrt(n_sim)
    rmse   <- sqrt(mse)
    ci_lo  <- sqrt(max(mse - 1.96 * se_mse, 0))
    ci_hi  <- sqrt(mse + 1.96 * se_mse)
    cat(sprintf("  RMSE = %.4f  [%.4f, %.4f]\n", rmse, ci_lo, ci_hi))
    c(rmse = rmse, ci_lo = ci_lo, ci_hi = ci_hi)
  }, numeric(3)))
}

# =============================================================================
# 5) Experiment 1: fix d=5, vary D
#    All curves should share slope ≈ −β/d = −0.2 (D is irrelevant)
# =============================================================================
cat(sprintf("\n=== Exp 1: d=%d fixed, D varies ===\n", d_fixed_exp1))
D_star_exp1 <- oracle_by_d[as.character(d_fixed_exp1)]

results_D_vary <- lapply(setNames(D_vary_vals, as.character(D_vary_vals)), function(D) {
  cat(sprintf("  D = %d\n", D))
  rmse_grid(d_val = d_fixed_exp1, D_val = D, Knn_val = K_fixed,
            D_star = D_star_exp1)
})

# =============================================================================
# 6) Experiment 2: fix D=20, vary d
#    Curves should show different slopes −min(β/d, 1/2); smaller d → steeper
# =============================================================================
cat(sprintf("\n=== Exp 2: D=%d fixed, d varies ===\n", D_fixed_exp2))

results_d_vary <- lapply(setNames(d_vary_vals, as.character(d_vary_vals)), function(d) {
  cat(sprintf("  d = %d\n", d))
  D_star_d <- oracle_by_d[as.character(d)]
  rmse_grid(d_val = d, D_val = D_fixed_exp2, Knn_val = K_fixed,
            D_star = D_star_d)
})

# =============================================================================
# 7) Experiment 3: fix d=5, D=20, vary K
#    All curves should share slope ≈ −0.2; larger K shifts curve up
# =============================================================================
cat(sprintf("\n=== Exp 3: d=%d, D=%d fixed, K varies ===\n",
            d_fixed_exp3, D_fixed_exp3))
D_star_exp3 <- oracle_by_d[as.character(d_fixed_exp3)]

results_K_vary <- lapply(setNames(K_vary_vals, as.character(K_vary_vals)), function(K) {
  cat(sprintf("  K = %d\n", K))
  rmse_grid(d_val = d_fixed_exp3, D_val = D_fixed_exp3, Knn_val = K,
            D_star = D_star_exp3)
})

# =============================================================================
# 8) Save
# =============================================================================
saveRDS(
  list(
    oracle_by_d    = oracle_by_d,
    n_vals         = n_vals,
    K_fixed        = K_fixed,
    d_fixed_exp1   = d_fixed_exp1,
    D_vary_vals    = D_vary_vals,
    results_D_vary = results_D_vary,
    D_fixed_exp2   = D_fixed_exp2,
    d_vary_vals    = d_vary_vals,
    results_d_vary = results_d_vary,
    d_fixed_exp3   = d_fixed_exp3,
    D_fixed_exp3   = D_fixed_exp3,
    K_vary_vals    = K_vary_vals,
    results_K_vary = results_K_vary
  ),
  file = file.path(out_dir,"convergence_rate.rds")
)
cat("Saved convergence_rate.rds\n")

parallel::stopCluster(cl)
