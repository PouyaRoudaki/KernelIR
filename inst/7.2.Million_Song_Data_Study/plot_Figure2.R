library(tidyr)
library(ggplot2)
library(cowplot)
library(tikzDevice)

#dir.create("../tikz_combined", recursive = TRUE, showWarnings = FALSE)

# === Method registry ===
all_methods <- c(
  "$T_n$", "$\\nu_n$",
  "$\\hat{\\eta}^{RKHS}$", "$\\hat{\\eta}^{K=5}$",
  "$\\hat{D}^{RKHS}$", "$\\hat{D}^{K=5}$"
)
#all_methods <- c("ti", "nu", "eraaaa", "egaaaa", "draaaa", "dgaaaa")

shape_vals <- c(18, 4, 16, 1, 17, 2)
names(shape_vals) <- all_methods

color_vals <- c(
  "#1f77b4", "#f68d62",
  "#1b9e77", "#7570b3", "#085f02", "#e7298a"
)
names(color_vals) <- all_methods

base_theme <- theme_bw(base_size = 9) +
  theme(
    axis.title        = element_text(size = 6),
    axis.text         = element_text(size = 5),
    legend.text       = element_text(size = 5),
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_blank(),
    plot.title        = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(0.5, 0.1, 0.1, 0.1)
  )

# === Helper ===
make_plot <- function(df_long, x_label, y_label, y_scale = "linear", show_legend = FALSE) {
  p <- ggplot(df_long, aes(x = n, y = value, group = method)) +
    geom_line(aes(color = method, linetype = method), linewidth = 0.2) +
    geom_point(aes(color = method, shape = method), size = 1, stroke = 0.7) +
    scale_shape_manual(values = shape_vals, breaks = all_methods, drop = FALSE) +
    scale_color_manual(values = color_vals, breaks = all_methods, drop = FALSE) +
    scale_linetype_manual(
      values = rep("solid", length(all_methods)),
      breaks = all_methods, drop = FALSE
    ) +
    labs(x = x_label, y = y_label, color = NULL, shape = NULL, linetype = NULL) +
    base_theme

  if (y_scale == "log") {
    p <- p + scale_y_log10() + scale_x_log10()
  } else {
    p <- p + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_x_log10()
  }

  if (!show_legend) p <- p + theme(legend.position = "none")
  else              p <- p + theme(legend.position = "bottom",
                                   legend.direction = "horizontal")
  p
}

# === Load power data (MSD_data) ===
#name_map_power <- c(
#  codec       = "ti",
#  irdc        = "nu",
#  kpcrkhs     = "eraaaa",
#  kmagraph_5nn= "egaaaa",
#  kir_rkhs    = "draaaa",
#  kir_graph_5nn = "dgaaaa"
#)

name_map_power <- c(
  codec_us        = "$T_n$",
  irdc_us         = "$\\nu_n$",
  kpc_rkhs_us     = "$\\hat{\\eta}^{RKHS}$",
  kpc_graph_5_us  = "$\\hat{\\eta}^{K=5}$",
  kir_rkhs_us     = "$\\hat{D}^{RKHS}$",
  kir_graph_5_us  = "$\\hat{D}^{K=5}$"
)

df_power <- MSD_data
keep_power <- c("n", names(name_map_power))
df_power <- df_power[, keep_power[keep_power %in% names(df_power)], drop = FALSE]
mc <- setdiff(names(df_power), "n")
names(df_power)[match(mc, names(df_power))] <- unname(name_map_power[mc])

df_power_long <- pivot_longer(df_power, -n, names_to = "method", values_to = "value")
df_power_long$method <- factor(df_power_long$method, levels = all_methods)

# === Load time data (time_df) ===

#name_map_time <- c(
#  codec_us        = "ti",
#  irdc_us         = "nu",
#  kpc_rkhs_us     = "eraaaa",
#  kpc_graph_5_us  = "egaaaa",
#  kir_rkhs_us     = "draaaa",
#  kir_graph_5_us  = "dgaaaa"
#)

name_map_time <- c(
  codec_us        = "$T_n$",
  irdc_us         = "$\\nu_n$",
  kpc_rkhs_us     = "$\\hat{\\eta}^{RKHS}$",
  kpc_graph_5_us  = "$\\hat{\\eta}^{K=5}$",
  kir_rkhs_us     = "$\\hat{D}^{RKHS}$",
  kir_graph_5_us  = "$\\hat{D}^{K=5}$"
)



df_time <- time_df
keep_time <- c("n", names(name_map_time))
df_time <- df_time[, keep_time[keep_time %in% names(df_time)], drop = FALSE]
mc <- setdiff(names(df_time), "n")
names(df_time)[match(mc, names(df_time))] <- unname(name_map_time[mc])

df_time_long <- pivot_longer(df_time, -n, names_to = "method", values_to = "value")
df_time_long$method <- factor(df_time_long$method, levels = all_methods)

# === Build plots ===
p_power <- make_plot(df_power_long, x_label = "Sample Size $n$",
                     y_label = "Power", y_scale = "linear", show_legend = FALSE)

p_time  <- make_plot(df_time_long,  x_label = "Sample Size $n$",
                     y_label = "Time ($\\mu$s)", y_scale = "log", show_legend = FALSE)

# === Legend from power plot (has all 6 methods) ===
p_legend_src <- make_plot(df_power_long, "", "", show_legend = TRUE) +
  guides(
    color    = guide_legend(nrow = 1, byrow = TRUE,
                            override.aes = list(size = 1.5, linewidth = 0.4)),
    shape    = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme(
    legend.position       = "bottom",
    legend.direction      = "horizontal",
    legend.box            = "horizontal",
    legend.text           = element_text(size = 5),
    legend.background     = element_rect(fill = "white", color = "white", linewidth = 1),
    legend.box.background = element_rect(fill = "white", color = "white", linewidth = 1),
    legend.key.width      = unit(0.4, "cm"),
    legend.key.height     = unit(0.35, "cm"),
    legend.spacing.x      = unit(-0.15, "cm"),
    legend.margin         = margin(4, 2, 4, 2)
  )

grob      <- ggplotGrob(p_legend_src)
leg_index <- which(sapply(grob$grobs, function(x) grepl("guide-box", x$name)))
if (length(leg_index) == 0) stop("Legend not found in grob!")
legend_row <- grob$grobs[[leg_index[1]]]

# === Assemble ===
combined <- plot_grid(
  legend_row,
  plot_grid(p_power, p_time, ncol = 2, align = "h"),
  ncol = 1,
  rel_heights = c(0.10, 1)
)

# === Output ===
tikz("power_MSD_time_final.tex",
     width = 3, height = 1.5, standAlone = FALSE)
print(combined)
dev.off()
