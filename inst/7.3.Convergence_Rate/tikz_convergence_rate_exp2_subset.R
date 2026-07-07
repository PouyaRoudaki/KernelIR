# =============================================================================
# TikZ version of convergence_rate_exp2_subset.pdf
# Same plot as plot_convergence_rate_exp2_subset.R, rendered through the
# tikzDevice so the figure is emitted as native TikZ (LaTeX) code.
# Output: convergence_rate_exp2_subset.tex  (+ raster .png sidecars if any)
# =============================================================================

library(ggplot2)
library(scales)
library(tikzDevice)

res <- readRDS("KernelIR/Supplementary_material/7.3.Convergence_Rate/convergence_rate.rds")

n_vals      <- res$n_vals
D_fixed     <- res$D_fixed_exp2
K_fixed     <- res$K_fixed

d_subset <- c(1L, 5L, 10L)

# =============================================================================
# Tidy data frame (subset only)
# =============================================================================
df_d <- do.call(rbind, lapply(as.character(d_subset), function(kv) {
  mat <- res$results_d_vary[[kv]]
  data.frame(key   = kv,
             n     = n_vals,
             rmse  = mat[, "rmse"],
             ci_lo = mat[, "ci_lo"],
             ci_hi = mat[, "ci_hi"],
             stringsAsFactors = FALSE)
}))

df_d$key <- factor(df_d$key, levels = as.character(d_subset))

d_colors <- setNames(c("#d62728", "#2ca02c", "#1f77b4"), as.character(d_subset))
d_shapes <- setNames(c(16, 15, 8),                       as.character(d_subset))

# TeX-friendly labels (tikzDevice passes these straight through to LaTeX)
d_label_exprs <- paste0("$d = ", d_subset, "$")
plot_title    <- sprintf("$d_0 = %d,\\ \\ K = %d$", D_fixed, K_fixed)

# =============================================================================
# Plot
# =============================================================================
p <- ggplot(df_d, aes(x = n, y = rmse, colour = key, shape = key)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = key),
              alpha = 0.12, color = NA, show.legend = FALSE) +
  scale_fill_manual(values = d_colors) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8, stroke = 0.6) +
  scale_colour_manual(values = d_colors, labels = d_label_exprs,
                      name = "$d$") +
  scale_shape_manual( values = d_shapes, labels = d_label_exprs,
                      name = "$d$") +
  scale_x_log10(breaks = n_vals, labels = as.character(n_vals),
                expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  labs(x     = "Sample size $n$",
       y     = "RMSE",
       title = plot_title) +
  theme_bw(base_size = 9) +
  theme(
    plot.title        = element_text(size = 8, hjust = 0.5),
    axis.title        = element_text(size = 7),
    axis.text         = element_text(size = 5.5),
    axis.text.x       = element_text(angle = 35, hjust = 1),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 6.5, face = "bold"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_blank(),
    legend.margin     = margin(2, 4, 2, 4),
    legend.position   = "right",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(4, 6, 2, 4)
  )

# =============================================================================
# Save as TikZ (same location / size as the .pdf)
# =============================================================================
tikz("convergence_rate_exp2_subset.tex", width = 3.0, height = 2.6,
     standAlone = TRUE)
print(p)
dev.off()
cat("Saved convergence_rate_exp2_subset.tex\n")
