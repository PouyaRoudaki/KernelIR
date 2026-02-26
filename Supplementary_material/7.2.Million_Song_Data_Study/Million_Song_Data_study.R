# =========================
# 0) Setup
# =========================

library(FORD)
library(FOCI)
library(KPC)
library(KernelIR)

library(doParallel)
library(foreach)

# =========================
# 1) Params
# =========================
n_sim   <- 100
B       <- 100

alpha   <- 0.01
eps_reg <- 1e-4
Knn     <- 5

sample_sizes <- c(50, 100, 150, 200, 250, 500, 750, 1000, 1250, 1500)

set.seed(7)

# parallel setup
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

# =========================
# 2) Load & preprocess full dataset
#    YearPredictionMSD: first column = year (Y), remaining 90 = features (X)
# =========================
cat("Loading YearPredictionMSD...\n")
msd_full <- read.csv("YearPredictionMSD.txt", header = FALSE)
#msd_full <- YearPredictionMSD
Y_full   <- as.numeric(msd_full[, 1])
X_full   <- as.matrix(msd_full[, -1])
N_full   <- nrow(X_full)
cat(sprintf("Full dataset: %d songs, %d features\n", N_full, ncol(X_full)))

# =========================
# 3) Test statistics
# =========================
stats_list <- local({
  .eps <- eps_reg
  .knn <- Knn
  list(
    codec         = function(x, y) FOCI::codec(Y = y, Z = x),
    irdc          = function(x, y) FORD::irdc(Y = y, X = x),

    kpcrkhs       = function(x, y) KPC::KPCRKHS(Y = y, Z = x, eps = .eps),
    kmagraph_5nn  = function(x, y) KPC::KMAc(Y = y, X = x, Knn = .knn),

    kir_rkhs      = function(x, y) KernelIR::kir_rkhs(Y = y, X = x, eps = .eps),
    kir_graph_5nn = function(x, y) KernelIR::kir_graph(Y = y, X = x, Knn = .knn)
  )
})
stat_names <- names(stats_list)

# =========================
# 4) Permutation p-value (right tail)
# =========================
perm_pval_righttail <- function(obs, null) {
  (1 + sum(null >= obs)) / (length(null) + 1)
}

# =========================
# 5) One-dataset rejection indicators
# =========================
one_dataset_rejections <- function(x, y, B, alpha, seed_perm = NULL) {
  if (!is.null(seed_perm)) set.seed(seed_perm)

  perm_idx <- replicate(B, sample.int(length(y)))
  obs      <- vapply(stats_list, function(f) f(x, y), numeric(1))

  null_mat <- matrix(NA_real_, nrow = B, ncol = length(stats_list))
  colnames(null_mat) <- stat_names

  for (b in seq_len(B)) {
    yb <- y[perm_idx[, b]]
    null_mat[b, ] <- vapply(stats_list, function(f) f(x, yb), numeric(1))
  }

  pvals <- vapply(seq_along(obs), function(j)
    perm_pval_righttail(obs[j], null_mat[, j]), numeric(1))
  names(pvals) <- stat_names

  as.numeric(pvals < alpha)
}

# =========================
# 6) Power estimation for a fixed sample size n
# =========================
compute_power_msd <- function(n, X_full, Y_full) {
  N_full <- nrow(X_full)

  results <- foreach(i = 1:n_sim, .combine = rbind,
                     .packages = c("XICOR", "FORD", "FOCI", "KPC", "KernelIR"),
                     .export   = c("stats_list", "stat_names",
                                   "one_dataset_rejections", "perm_pval_righttail",
                                   "B", "alpha", "eps_reg", "Knn")) %dopar% {

                                     idx <- sample.int(N_full, size = n, replace = FALSE)
                                     x   <- X_full[idx, , drop = FALSE]
                                     y   <- Y_full[idx]

                                     rej <- one_dataset_rejections(x, y, B = B, alpha = alpha, seed_perm = 200000 + i)
                                     names(rej) <- stat_names
                                     rej
                                   }

  colMeans(results)
}

# =========================
# 7) Main loop over sample sizes
# =========================
powers_list <- vector("list", length(sample_sizes))
names(powers_list) <- as.character(sample_sizes)

for (n in sample_sizes) {
  cat(sprintf("Running n = %d\n", n))
  powers_list[[as.character(n)]] <- compute_power_msd(n, X_full, Y_full)
}

powers_df <- do.call(rbind, powers_list)
powers_df <- data.frame(
  n = sample_sizes,
  powers_df,
  row.names = NULL
)

dir.create("resultsMSD", showWarnings = FALSE)
saveRDS(powers_df, file = paste("resultsMSD/power_MSD_n",n,".rds"))
#cat("Done. Results saved to results/power_MSD.rds\n")
print(powers_df)

stopCluster(cl)
