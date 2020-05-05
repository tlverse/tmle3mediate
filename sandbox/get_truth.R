library(here)
library(tidyverse)
source(here("dgp_funs.R"))
source(here("01_setup_data.R"))

get_truth <- function(n_obs = 1e7) { # number of observations
  # compute large data set for true values
  dgp <- make_dgp()
  data_truth <- sim_data(n_obs)

  w_names <- str_subset(colnames(data_truth), "W")
  z_names <- str_subset(colnames(data_truth), "Z")
  W <- data_truth[, w_names]
  Z <- data_truth[, z_names]

  Z_1_probs <- dgp$z_mech(W, 1)
  Z_0_probs <- dgp$z_mech(W, 0)

  # Z1 counterfactuals
  Z1_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[1]])
  Z1_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[1]])

  # Z2 counterfactuals
  Z2_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[2]])
  Z2_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[2]])

  # Z3 counterfactuals
  Z3_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[3]])
  Z3_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[3]])

  Z1 <- cbind(Z1_1, Z2_1, Z3_1)
  Z0 <- cbind(Z1_0, Z2_0, Z3_0)

  # compute Y counterfactuals
  # Y(1, 1)
  EY_A1_Z1 <- dgp$m_mech(W, 1, Z1, eps_sd = 0)
  # Y(0, 0)
  EY_A0_Z0 <- dgp$m_mech(W, 0, Z0, eps_sd = 0)
  # Y(1, 0)
  EY_A1_Z0 <- dgp$m_mech(W, 1, Z0, eps_sd = 0)
  # Y(0, 1)
  EY_A0_Z1 <- dgp$m_mech(W, 0, Z1, eps_sd = 0)

  # compute TRUE M under counterfactual regimes
  m_Ais1 <- dgp$m_mech(W, 1, Z, eps_sd = 0)
  m_Ais0 <- dgp$m_mech(W, 0, Z, eps_sd = 0)

  # output: true values of nuisance parameters
  return(list(
    EY_A1_Z1 = EY_A1_Z1,
    EY_A1_Z0 = EY_A1_Z0,
    EY_A0_Z1 = EY_A0_Z1,
    EY_A0_Z0 = EY_A0_Z0
  ))
}
