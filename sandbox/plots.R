library(tidyverse)
library(here)
library(ggpubr)

source(here("get_truth.R"))

sim_truth <- get_truth(n_obs = 1e7)
EY_A1_Z1 <- sim_truth$EY_A1_Z1
EY_A1_Z0 <- sim_truth$EY_A1_Z0
EY_A0_Z1 <- sim_truth$EY_A0_Z1
EY_A0_Z0 <- sim_truth$EY_A0_Z0
# compute true NIE via empirical substitution estimator
psi_NDE_true <- mean(EY_A1_Z0 - EY_A0_Z0)
psi_NIE_true <- mean(EY_A1_Z1 - EY_A1_Z0)

sim_results <- readRDS(here('data/tmle3mediate_2020-05-01_11:14:45.rds'))

# combine simulation results into one big df with n_obs column
sim_results <- names(sim_results) %>% map_dfr(
  function(this_n) {
    df <- sim_results[[this_n]]
    n_num <- as.numeric(stringr::str_replace(this_n, "n_", ""))
    cbind(n_obs = n_num, df)
  }
)

rename_sim_types <- function(sym_types) {
  sapply(
    sym_types, switch,
    corr_NDE = "Efficient",
    mis_e_NDE = "Eff.(E--mis.)",
    mis_z_NDE = "Eff.(M--mis.)",
    corr_NIE = "Efficient",
    mis_e_NIE = "Eff.(E--mis.)",
    mis_z_NIE = "Eff.(M--mis.)"
  )
}

sim_statistics <- sim_results %>%
  group_by(n_obs, type, sim_type) %>%
  summarize(
    est_bias = if (first(type) == "NDE") {
      mean(tmle_est) - psi_NDE_true
    } else {
      mean(tmle_est) - psi_NIE_true
    },
    est_sd = sd(tmle_est),
    est_MSE = est_bias^2 + est_sd^2,
  ) %>%
  mutate(
    sim_type = rename_sim_types(sim_type),
    ## stats for plotting
    abs_bias = abs(est_bias),
    sqrt_n_abs_bias = sqrt(n_obs) * abs(est_bias),
    sqrt_n_sd = sqrt(n_obs) * est_sd,
    n_MSE = n_obs*est_MSE
  ) %>%
  ungroup

sim_statistics_long <- sim_statistics %>%
  gather(statistic, value, -c(n_obs, type, sim_type))

n_obs <- (cumsum(rep(sqrt(100), 8))^2)[-1]

make_sim_statistics_plot <- function(sim_statistics_long, est_type) {
  filtered_sim_stats <- sim_statistics_long %>%
    filter(
      type == est_type,
      statistic %in% c("abs_bias", "sqrt_n_abs_bias", "sqrt_n_sd", "n_MSE")
    )

  ggplot(filtered_sim_stats, aes(n_obs, value)) +
    geom_point() + geom_line(linetype = "dashed") +
    facet_grid(statistic ~ sim_type, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=45, vjust=1, hjust=1)
    ) +
    scale_x_continuous(breaks=n_obs)
}

NDE_stats_plot <- make_sim_statistics_plot(sim_statistics_long, 'NDE') +
  ggtitle('NDE')

NIE_stats_plot <- make_sim_statistics_plot(sim_statistics_long, 'NIE') +
  ggtitle('NIE')

ggsave(
  'plots/simulation_statistics.png',
  ggarrange(NDE_stats_plot, NIE_stats_plot, ncol = 2),
  width = 13,
  device = png()
)
