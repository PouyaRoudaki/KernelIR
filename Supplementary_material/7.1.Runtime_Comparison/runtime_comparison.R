# ── Libraries ─────────────────────────────────────────────────────────────────
library(KernelIR)
library(XICOR)
library(FORD)
library(FOCI)
library(KPC)
library(foreach)
library(doParallel)

options(digits.secs = 6)   # ensure microsecond printing

# ── Settings ───────────────────────────────────────────────────────────────────
n_grid  <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
reps    <- 100L
seed0   <- 1L
out_dir <- "results_asymp_time"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Microsecond timer using Sys.time() ─────────────────────────────────────────
time_safe <- function(expr) {
  t0 <- as.numeric(Sys.time())
  tryCatch(expr, error = function(e) NULL)
  (as.numeric(Sys.time()) - t0)    # microseconds
}

# ── One replicate: returns named vector of timings (microseconds) ─────────────
one_rep <- function(n, i, seed0) {
  set.seed(seed0 + i + 10000L * n)
  X <- rnorm(n)
  Y <- rnorm(n)

  c(
    xicor_us       = time_safe(xicor(X, Y)),
    codec_us       = time_safe(codec(Y = Y, Z = X)),
    irdc_us        = time_safe(irdc(Y = Y, X = X)),
    kpc_graph_1_us = time_safe(KPCgraph(Y = Y, Z = X, X = NULL, Knn = 1L)),
    kpc_graph_5_us = time_safe(KPCgraph(Y = Y, Z = X, X = NULL, Knn = 5L)),
    kpc_rkhs_us    = time_safe(KPCRKHS(Y = Y, Z = X, X = NULL)),
    kir_graph_1_us = time_safe(kir_graph(Y = Y, X = X, Knn = 1L)),
    kir_graph_5_us = time_safe(kir_graph(Y = Y, X = X, Knn = 5L)),
    kir_rkhs_us    = time_safe(kir_rkhs(Y = Y, X = X))
  )
}

# ── Parallel setup ─────────────────────────────────────────────────────────────
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
message("Using ", n_cores, " core(s).")

cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

parallel::clusterEvalQ(cl, {
  library(XICOR)
  library(FORD)
  library(FOCI)
  library(KPC)
  library(KernelIR)
  NULL
})

parallel::clusterExport(cl, varlist = c("time_safe", "one_rep", "seed0"))

# ── Main loop: compute + save per n only ───────────────────────────────────────
for (n in n_grid) {

  message(format(Sys.time(), "%H:%M:%OS6"), "  n = ", n)

  mat <- tryCatch(
    foreach(
      i         = seq_len(reps),
      .combine  = "rbind",
      .inorder  = FALSE,
      .packages = c("XICOR", "FORD", "FOCI", "KPC", "KernelIR")
    ) %dopar% {
      one_rep(n, i, seed0)
    },
    error = function(e) {
      warning("n = ", n, " failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(mat)) next

  res_n <- cbind(
    n   = n,
    rep = seq_len(nrow(mat)),
    mat
  )

  saveRDS(
    res_n,
    file = file.path(out_dir, sprintf("results_n%05d.rds", n))
  )

  rm(res_n, mat)
  gc(FALSE)

  message("  -> saved results_n", sprintf("%05d", n), ".rds")
}

# ── Save metadata ──────────────────────────────────────────────────────────────
meta <- list(
  n_grid    = n_grid,
  reps      = reps,
  seed0     = seed0,
  n_cores   = n_cores,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%OS6")
)

saveRDS(meta, file = file.path(out_dir, "meta.rds"))

message("All done. Results saved to ", out_dir)

parallel::stopCluster(cl)
