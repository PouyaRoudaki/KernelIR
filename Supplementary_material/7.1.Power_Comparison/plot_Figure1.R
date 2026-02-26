library(tidyr)
library(ggplot2)
library(tikzDevice)
library(cowplot)

#dir.create("tikz_combined", recursive = TRUE, showWarnings = FALSE)

# === Shared aesthetics ===
all_methods <- c(
  "$\\xi_n$", "$T_n$", "$\\nu_n$",
  "$\\hat{\\eta}^{RKHS}$", "$\\hat{\\eta}^{K=5}$",
  "$\\hat{D}^{RKHS}$", "$\\hat{D}^{K=5}$"
)

shape_vals <- c(15, 18, 4, 16, 1, 17, 2)
names(shape_vals) <- all_methods

color_vals <- c(
  "#66a61e", "#1f77b4", "#f68d62",
  "#1b9e77", "#7570b3", "#085f02", "#e7298a"
)
names(color_vals) <- all_methods

base_theme <- theme_bw(base_size = 9) +
  theme(
    axis.title        = element_text(size = 6),
    axis.text         = element_text(size = 5),
    legend.text       = element_text(size = 15),
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_blank(),
    plot.title        = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),

    # minimize outer margins: (t, r, b, l)
    plot.margin       = margin(0.5, 0.1, 0.1, 0.1)
  )

# === Helper ===
make_plot <- function(df_long, x_label, show_legend = FALSE) {
  p <- ggplot(df_long, aes(x = lambda, y = power, group = method)) +
    geom_line(aes(color = method, linetype = method), linewidth = 0.2) +
    geom_point(aes(color = method, shape = method), size = 1, stroke = 0.7) +
    scale_shape_manual(values = shape_vals, breaks = all_methods, drop = FALSE) +
    scale_color_manual(values = color_vals, breaks = all_methods, drop = FALSE) +
    scale_linetype_manual(
      values = c("solid","solid","solid","solid","solid","solid","solid"),
      breaks = all_methods, drop = FALSE
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = x_label, y = "Power", color = NULL, shape = NULL, linetype = NULL) +
    base_theme

  if (!show_legend) p <- p + theme(legend.position = "none")
  else              p <- p + theme(legend.position = "bottom",
                                   legend.direction = "horizontal")
  p
}

# === Load heteroskedastic ===
keep_hetero <- c("lambda","xi","codec","irdc","kpcrkhs","kmagraph_5nn","kir_rkhs","kir_graph_5nn")
name_map_h  <- c(xi="$\\xi_n$", codec="$T_n$", irdc="$\\nu_n$",
                 kpcrkhs="$\\hat{\\eta}^{RKHS}$", kmagraph_5nn="$\\hat{\\eta}^{K=5}$",
                 kir_rkhs="$\\hat{D}^{RKHS}$", kir_graph_5nn="$\\hat{D}^{K=5}$")
name_map_h  <- c(
  xi = "$\\xi_n$",
  codec = "$T_n$",
  irdc = "$\\nu_n$",
  kpcrkhs = "$\\hat{\\eta}^{RKHS}$",
  kmagraph_5nn = "$\\hat{\\eta}^{K=5}$",
  kir_rkhs = "$\\hat{D}^{RKHS}$",
  kir_graph_5nn = "$\\hat{D}^{K=5}$"
)


df_h <- readRDS("power_heteroscedastic.rds")
df_h <- df_h[, keep_hetero[keep_hetero %in% names(df_h)], drop = FALSE]
mc <- setdiff(names(df_h), "lambda")
names(df_h)[match(mc, names(df_h))] <- unname(name_map_h[mc])
df_h_long <- pivot_longer(df_h, -lambda, names_to = "method", values_to = "power")
df_h_long$method <- factor(df_h_long$method, levels = all_methods)

# === Load SO3 ===
keep_so3   <- c("lambda","kpcrkhs","kmagraph_5nn","kir_rkhs","kir_graph_5nn")
name_map_s <- c(kpcrkhs="$\\hat{\\eta}^{RKHS}$", kmagraph_5nn="$\\hat{\\eta}^{K=5}$",
                kir_rkhs="$\\hat{D}^{RKHS}$", kir_graph_5nn="$\\hat{D}^{K=5}$")
name_map_s <- c(
  kpcrkhs = "$\\hat{\\eta}^{RKHS}$",
  kmagraph_5nn = "$\\hat{\\eta}^{K=5}$",
  kir_rkhs = "$\\hat{D}^{RKHS}$",
  kir_graph_5nn = "$\\hat{D}^{K=5}$"
)

df_s <- readRDS("power_SO3.rds")
df_s <- df_s[, keep_so3[keep_so3 %in% names(df_s)], drop = FALSE]
mc <- setdiff(names(df_s), "lambda")
names(df_s)[match(mc, names(df_s))] <- unname(name_map_s[mc])
df_s_long <- pivot_longer(df_s, -lambda, names_to = "method", values_to = "power")
df_s_long$method <- factor(df_s_long$method, levels = all_methods)

# === Build plots ===
p_h <- make_plot(df_h_long, "Homoscedasticity Level $(\\lambda)$", show_legend = FALSE)
p_s <- make_plot(df_s_long, "Noise Scale $(\\lambda)$",            show_legend = FALSE)

# === Legend source: use df_h_long which has ALL 7 methods ===
p_legend_src <- make_plot(df_h_long, "", show_legend = TRUE) +
  guides(
    color    = guide_legend(nrow = 1, byrow = TRUE,
                            override.aes = list(size = 1.5, linewidth = 0.4)),
    shape    = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme(
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.box           = "horizontal",   # important
    legend.text          = element_text(size = 5),
    legend.background    = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.box.background= element_rect(fill = "white", color = "black", linewidth = 1),
    legend.key.width     = unit(0.4, "cm"),
    legend.key.height    = unit(0.35, "cm"),
    legend.spacing.x     = unit(-0.15, "cm"),
    legend.margin        = margin(4, 2, 4, 2)
  )

# Robust legend extraction
grob      <- ggplotGrob(p_legend_src)
leg_index <- which(sapply(grob$grobs, function(x) grepl("guide-box", x$name)))
if (length(leg_index) == 0) stop("Legend not found in grob!")
legend_row <- grob$grobs[[leg_index[1]]]

# Assemble: legend on top, plots below
combined <- plot_grid(
  legend_row,
  plot_grid(p_h, p_s, ncol = 2, align = "h"),
  ncol = 1,
  rel_heights = c(0.10, 1)
)

# === Output single tikz file ===
tikz("power_comparison.tex",
     width = 3, height = 1.5, standAlone = FALSE)
print(combined)
dev.off()


