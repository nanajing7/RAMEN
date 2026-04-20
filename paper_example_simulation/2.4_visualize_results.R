library(ggplot2)
library(patchwork)
library(scales)

# -----------------------------
# Load data
# -----------------------------
save_dir  <- "xxxx"
load("part2a_varyGamma.RData")
load("part2b_varyT.RData")
load("part2c_varyMask.RData")
mask_summary <- mask_summary[mask_summary$mask_frac <= 0.8,]

# -----------------------------
# Common theme
# -----------------------------
my_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    panel.grid.minor = element_line(color = "grey85"),
    panel.grid.major = element_line(color = "grey80")
  )

# -----------------------------
# Panel (a): gamma
# -----------------------------
p_gamma <- ggplot(rmse_summary, aes(x = gamma, y = mean_rmse)) +
  geom_line(linewidth = 1.2, color = "#1f5aa6") +
  geom_point(size = 3.5, color = "#1f5aa6") +
  theme_minimal(base_size = 16) +
  labs(
    x = expression(gamma),
    y = "Out-of-sample RMSE",
    title = "(a) Temporal Smoothing"
  ) +
  my_theme

# -----------------------------
# Panel (b): T
# -----------------------------
p_T <- ggplot(T_summary, aes(x = T_periods, y = mean_rmse)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, fill = "#1f5aa6") +
  geom_line(linewidth = 1.2, color = "#1f5aa6") +
  geom_point(size = 3.5, color = "#1f5aa6") +
  scale_x_continuous(breaks = T_grid) +
  theme_minimal(base_size = 16) +
  labs(
    x = "T",
    y = NULL,
    title = "(b) Time Horizon"
  ) +
  my_theme

# -----------------------------
# Panel (c): mask fraction
# -----------------------------
p_mask <- ggplot(mask_summary, aes(x = mask_frac, y = mean_rmse)) +
  geom_line(linewidth = 1.2, color = "#1f5aa6") +
  geom_point(size = 3.5, color = "#1f5aa6") +
  scale_x_continuous(
    breaks = mask_frac_grid,
    labels = percent(mask_frac_grid)
  ) +
  theme_minimal(base_size = 16) +
  labs(
    x = "% held-out dyads",
    y = NULL,
    title = "(c) Holdout Fraction"
  ) +
  my_theme

# -----------------------------
# Combine panels
# -----------------------------
p_combined <- (p_gamma | p_T | p_mask) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 24, hjust = 0.5)
    )
  )

p_combined
