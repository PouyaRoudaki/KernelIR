# =========================
# 0) Setup (your libraries)
# =========================
library(XICOR)
library(FOCI)
library(FORD)
library(KPC)
library(KernelIR)

library(doParallel)
library(foreach)

# =========================
# 1) Params
# =========================
n <- 100
n_sim <- 500
B <- 1000                 # number of permutations per dataset

#n <- 10
#n_sim <- 10
#B <- 10                 # number of permutations per dataset

alpha <- 0.05

lambda_vals <- seq(0, 1, 0.1)
alternatives <- c("heteroscedastic")

set.seed(7)

# parallel setup: respect Slurm allocation
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

# =========================
# 2) Data generator (yours)
# =========================
generate_data <- function(type, lambda, seed) {
  set.seed(seed)

  x <- runif(n, -1, 1)
  eps <- rnorm(n)

  y <- switch(type,
              linear = 0.5 * x + 3 * lambda * eps,

              step = {
                f <- numeric(n)
                f[x < -0.5] <- -3
                f[x >= -0.5 & x < 0] <- 2
                f[x >= 0 & x < 0.5] <- -4
                f[x >= 0.5] <- -3
                f + 10 * lambda * eps
              },

              W = abs(x + 0.5) * (x < 0) + abs(x - 0.5) * (x >= 0) + 0.75 * lambda * eps,

              Sin = cos(8 * pi * x) + 3 * lambda * eps,

              circular = {
                Z <- sample(c(-1, 1), n, replace = TRUE)
                Z * sqrt(pmax(0, 1 - x^2)) + 0.9 * lambda * eps
              },

              spiral = {
                theta <- runif(n, 0, 4*pi)
                r <- theta / (4*pi)
                x2 <- r * cos(theta) + lambda * 0.1 * eps
                y2 <- r * sin(theta) + lambda * 0.1 * eps
                x <<- x2
                y2
              },

              x_shape = {
                Z <- sample(c(-1, 1), n, replace = TRUE)
                Z * x + 0.1 * lambda * eps
              },

              two_moons = {
                theta1 <- runif(n/2, 0, pi)
                theta2 <- runif(n/2, 0, pi)
                x2 <- c(cos(theta1), 1 - cos(theta2))
                y2 <- c(sin(theta1), -sin(theta2) + 0.5)
                x2 <- x2 + lambda * rnorm(n)
                y2 + lambda * rnorm(n)
              },

              two_sines = {
                Z <- sample(c(-1, 1), n, replace = TRUE)
                Z * sin(4 * pi * x) + 0.1 * lambda * eps
              },

              lissajous_curve = {
                A <- 1; Bc <- 1
                a <- 3; b <- 2
                delta <- pi / 2
                t <- runif(n, 0, 2*pi)
                x2 <- A * sin(a * t + delta) + lambda * 0.1 * eps
                y2 <- Bc * sin(b * t) + lambda * 0.1 * eps
                x <<- x2
                y2
              },

              heteroscedastic = {
                sigma_x <- ifelse(abs(x) <= 0.5, 1, 0)
                3 * (sigma_x * (1 - lambda) + lambda) * eps
              },

              heteroscedastic_sin = {
                cos(20*pi*(1 + 10 * lambda * eps) * x^2)
              },

              hetero_band = {
                base_sigma <- 0.3
                high_sigma <- 2.4
                sigma <- ifelse(abs(x) < 0.2,
                                (1 - lambda) * high_sigma + lambda * base_sigma,
                                base_sigma)
                rnorm(n, 0, sigma)
              },

              hetero_smooth = {
                base_sigma <- 0.25
                bump <- 2.25 * exp(-(x - 0.2)^2 / (2 * 0.18^2))
                sigma <- base_sigma + (1 - lambda) * bump
                rnorm(n, 0, sigma)
              },

              stop("Invalid type")
  )

  list(X = x, Y = y)
}

# =========================
# 3) Generic permutation p-value
# =========================
perm_pval_righttail <- function(obs, null) {
  # "greater" tail: large statistic => dependence
  # include obs in denominator for exactness
  (1 + sum(null >= obs)) / (length(null) + 1)
}

# =========================
# 4) One place to define ALL test statistics
#    Each returns a *single numeric*.
# =========================
stats_list <- list(
  xi        = function(x, y) xicor(y = y, x = x),

  codec     = function(x, y) FOCI::codec(Y = y, Z = x),

  irdc      = function(x, y) FORD::irdc(Y = y, X = x),

  kpcrkhs   = function(x, y) KPC::KPCRKHS(Y = y, Z = x, eps = 1e-4),
  kmagraph_5nn  = function(x, y) KPC::KMAc(Y = y, X = x, Knn = 5),

  kir_rkhs      = function(x, y) KernelIR::kir_rkhs(Y = y, X = x, eps = 1e-4),
  kir_graph_5nn = function(x, y) KernelIR::kir_graph(Y = y, X = x, Knn = 5)
)
stat_names <- names(stats_list)

# =========================
# 5) Core engine: compute rejections for ONE dataset using ONE set of perms
# =========================
one_dataset_rejections <- function(x, y, B, alpha, seed_perm = NULL) {
  if (!is.null(seed_perm)) set.seed(seed_perm)

  # permutation indices: n x B
  perm_idx <- replicate(B, sample.int(length(y)))

  # observed stats
  obs <- vapply(stats_list, function(f) f(x, y), numeric(1))

  # null stats: B x M
  null_mat <- matrix(NA_real_, nrow = B, ncol = length(stats_list))
  colnames(null_mat) <- stat_names

  for (b in seq_len(B)) {
    yb <- y[perm_idx[, b]]
    null_mat[b, ] <- vapply(stats_list, function(f) f(x, yb), numeric(1))
  }

  # right-tail permutation p-values (consistent across all tests)
  pvals <- vapply(seq_along(obs), function(j) perm_pval_righttail(obs[j], null_mat[, j]), numeric(1))
  names(pvals) <- stat_names

  # return rejection indicators
  as.numeric(pvals < alpha)
}

# =========================
# 6) Simulation wrapper: returns power estimates (colMeans)
# =========================
compute_power <- function(type, lambda) {
  results <- foreach(i = 1:n_sim, .combine = rbind,
                     .packages = c("XICOR","FOCI","FORD","KPC","KernelIR"),
                     .export = c("generate_data","stats_list","stat_names",
                                 "one_dataset_rejections","perm_pval_righttail",
                                 "n","B","alpha")) %dopar% {

                                   dat <- generate_data(type, lambda, seed = i)
                                   x <- dat$X
                                   y <- dat$Y

                                   # IMPORTANT: one permutation bank per dataset
                                   # If you want reproducible across cores, tie perm seed to i (and maybe type/lambda)
                                   rej <- one_dataset_rejections(x, y, B = B, alpha = alpha, seed_perm = 100000 + i)

                                   # return named vector
                                   out <- rej
                                   names(out) <- stat_names
                                   out
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

stopCluster(cl)



