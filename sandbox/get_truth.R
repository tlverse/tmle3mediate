library(here)
library(tidyverse)
source(here("dgp_funs.R"))
source(here("01_setup_data.R"))

get_truth <- function(n_obs = 1e7) {
  ## set up DGP
  dgp <- make_dgp()
  data_truth <- sim_data(n_obs)
  W <- data_truth[, str_subset(names(data_truth), "W")]
  A <- data_truth$A
  Z <- data_truth[, str_subset(names(data_truth), "Z")]

  ## generate mediator counterfactuals
  z_mech_nat <- dgp$z_mech(a = data_truth$A, w = W)
  z_mech_a1 <- dgp$z_mech(a = 1, w = W)
  z_mech_a0 <- dgp$z_mech(a = 0, w = W)
  Z_1_1 <- rbinom(n_obs, 1, prob = z_mech_a1[[1]])
  Z_1_0 <- rbinom(n_obs, 1, prob = z_mech_a0[[1]])
  Z_2_1 <- rbinom(n_obs, 1, prob = z_mech_a1[[2]])
  Z_2_0 <- rbinom(n_obs, 1, prob = z_mech_a1[[2]])
  Z_3_1 <- rbinom(n_obs, 1, prob = z_mech_a1[[3]])
  Z_3_0 <- rbinom(n_obs, 1, prob = z_mech_a1[[3]])

  ## compute counterfactuals Y(1) and Y(0)
  EY_A1_Z1 <- dgp$m_mech(z = list(Z_1_1, Z_2_1, Z_3_1), a = 1, w = W)
  EY_A0_Z0 <- dgp$m_mech(z = list(Z_1_0, Z_2_0, Z_3_0), a = 0, w = W)
  EY_A1_Z0 <- dgp$m_mech(z = list(Z_1_0, Z_2_0, Z_3_0), a = 1, w = W)
  EY_A0_Z1 <- dgp$m_mech(z = list(Z_1_1, Z_2_1, Z_3_1), a = 0, w = W)

  ## compute true mediation mechanism under counterfactual regimes
  m_Ais1 <- dgp$m_mech(z = Z, a = 1, w = W)
  m_Ais0 <- dgp$m_mech(z = Z, a = 0, w = W)

  ## output
}
