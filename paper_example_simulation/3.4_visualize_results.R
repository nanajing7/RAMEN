library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(patchwork)

rm(list = ls())

# ------------------------------------------------------------
# Style settings (colorblind-friendly)
# ------------------------------------------------------------
method_colors <- c(
  "ALS + bootstrap" = "#0072B2",
  "MCMC" = "#D55E00"
)

method_shapes <- c(
  "ALS + bootstrap" = 16,
  "MCMC" = 15
)

method_linetypes <- c(
  "ALS + bootstrap" = "solid",
  "MCMC" = "dashed"
)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
save_dir <- "xxx"
load("part3a.RData")
load("part3b.RData")
load("part3c.RData")

# restrict K
timing_C <- timing_C[timing_C$K <= 3, ]

# ------------------------------------------------------------
# Common y-axis
# ------------------------------------------------------------
y_max <- max(
  timing_A$t_ramen_boot, timing_A$t_amen,
  timing_B$t_ramen_boot, timing_B$t_amen,
  timing_C$t_ramen_boot, timing_C$t_amen,
  na.rm = TRUE
)

y_min <- 0

# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------
plot_df_A <- timing_A %>%
  select(N, t_ramen_boot, t_amen) %>%
  pivot_longer(
    cols = c(t_ramen_boot, t_amen),
    names_to = "method",
    values_to = "elapsed"
  ) %>%
  mutate(
    method = recode(
      method,
      t_ramen_boot = "ALS + bootstrap",
      t_amen = "MCMC"
    ),
    method = factor(method, levels = c("MCMC", "ALS + bootstrap"))
  )

p1 <- ggplot(plot_df_A,
             aes(x = N, y = elapsed,
                 color = method,
                 linetype = method,
                 shape = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_linetypes) +
  labs(
    x = "Number of row nodes (N)",
    y = "Computation time (seconds)",
    color = NULL, shape = NULL, linetype = NULL,
    title = "(a) Network Size"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 18)
  )

# ------------------------------------------------------------
# Panel B
# ------------------------------------------------------------
plot_df_B <- timing_B %>%
  select(TT, t_ramen_boot, t_amen) %>%
  pivot_longer(
    cols = c(t_ramen_boot, t_amen),
    names_to = "method",
    values_to = "elapsed"
  ) %>%
  mutate(
    method = recode(
      method,
      t_ramen_boot = "ALS + bootstrap",
      t_amen = "MCMC"
    ),
    method = factor(method, levels = c("MCMC", "ALS + bootstrap"))
  )

p2 <- ggplot(plot_df_B,
             aes(x = TT, y = elapsed,
                 color = method,
                 linetype = method,
                 shape = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_linetypes) +
  labs(
    x = "Number of time periods (T)",
    y = NULL,
    color = NULL, shape = NULL, linetype = NULL,
    title = "(b) Time Horizon"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 18)
  )

# ------------------------------------------------------------
# Panel C
# ------------------------------------------------------------
plot_df_C <- timing_C %>%
  select(K, t_ramen_boot, t_amen) %>%
  pivot_longer(
    cols = c(t_ramen_boot, t_amen),
    names_to = "method",
    values_to = "elapsed"
  ) %>%
  mutate(
    method = recode(
      method,
      t_ramen_boot = "ALS + bootstrap",
      t_amen = "MCMC"
    ),
    method = factor(method, levels = c("MCMC", "ALS + bootstrap"))
  )

p3 <- ggplot(plot_df_C,
             aes(x = K, y = elapsed,
                 color = method,
                 linetype = method,
                 shape = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_linetypes) +
  labs(
    x = "Number of latent factors (K)",
    y = NULL,
    color = NULL, shape = NULL, linetype = NULL,
    title = "(c) Latent Dimension"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 18)
  )

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------
combined_plot <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    ##title = "Computational Efficiency: MCMC vs. ALS",
    theme = theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(size = 18, hjust = 0.5)
    )
  )

combined_plot
