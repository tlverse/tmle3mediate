context("TML estimator for natural (in)direct effects")

library(data.table)
library(stringr)
library(dplyr)
library(hal9001)
library(sl3)
library(tmle3)

# source in utils.R
testthat::source_test_helpers()

set.seed(7128816)

################################################################################
# setup learners for the nuisance parameters
################################################################################

# instantiate some learners
mean_lrnr <- Lrnr_mean$new()
fglm_contin_lrnr <- Lrnr_glm_fast$new()
fglm_binary_lrnr <- Lrnr_glm_fast$new(family = binomial())
hal_contin_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5
)
hal_binary_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5,
  family = "binomial"
)
cv_hal_contin_lrnr <- Lrnr_cv$new(hal_contin_lrnr, full_fit = TRUE)
cv_hal_binary_lrnr <- Lrnr_cv$new(hal_binary_lrnr, full_fit = TRUE)

################################################################################
# setup data and simulate to test with estimators

# get data and column names for sl3 tasks (for convenience)
data <- make_simulated_data()
z_names <- colnames(data)[str_detect(colnames(data), "Z")]
w_names <- colnames(data)[str_detect(colnames(data), "W")]

# create node list and learner list
node_list <- list(
  W = c("W_1", "W_2", "W_3"),
  A = "A",
  Z = c("Z_1", "Z_2", "Z_3"),
  Y = "Y"
)
learner_list <- list(
  Y = cv_hal_contin_lrnr,
  A = cv_hal_binary_lrnr
)

################################################################################

## instantiate tmle3 spec for NIE
tmle_spec_NIE <- tmle_NIE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 100 # TODO: use default when convergence bug fixed
)

set.seed(71281)
## define data (from tmle3_Spec base class)
tmle_task_NIE <- tmle_spec_NIE$make_tmle_task(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init_NIE <- tmle_spec_NIE$make_initial_likelihood(tmle_task_NIE, learner_list)

## define update method (submodel and loss function)
updater_NIE <- tmle_spec_NIE$make_updater()
likelihood_targeted_NIE <- Targeted_Likelihood$new(likelihood_init_NIE, updater_NIE)

## define param
tmle_params_NIE <- tmle_spec_NIE$make_params(tmle_task_NIE, likelihood_targeted_NIE)
updater_NIE$tmle_params <- tmle_params_NIE

## fit tmle update
tmle_fit_NIE <- fit_tmle3(tmle_task_NIE, likelihood_targeted_NIE, tmle_params_NIE, updater_NIE)

## one-line call with faster with tmle3 wrapper
set.seed(71281)
tmle_fit_NIE <- tmle3(tmle_spec_NIE, data, node_list, learner_list)

################################################################################

## instantiate tmle3 spec for NDE
tmle_spec_NDE <- tmle_NDE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 100 # TODO: use default when convergence bug fixed
)

## define data (from tmle3_Spec base class)
tmle_task_NDE <- tmle_spec_NDE$make_tmle_task(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init_NDE <- tmle_spec_NDE$make_initial_likelihood(tmle_task_NDE, learner_list)

## define update method (submodel and loss function)
updater_NDE <- tmle_spec_NDE$make_updater()
likelihood_targeted_NDE <- Targeted_Likelihood$new(likelihood_init_NDE, updater_NDE)

## define param
tmle_params_NDE <- tmle_spec_NDE$make_params(tmle_task_NDE, likelihood_targeted_NDE)
updater_NDE$tmle_params <- tmle_params_NDE

## fit tmle update
tmle_fit_NDE <- fit_tmle3(tmle_task_NDE, likelihood_targeted_NDE, tmle_params_NDE, updater_NDE)

## one-line call with faster with tmle3 wrapper
set.seed(71281)
tmle_fit_NDE <- tmle3(tmle_spec_NDE, data, node_list, learner_list)

################################################################################
## instantiate tmle3 spec for ATE
tmle_spec_ATE <- tmle3::tmle_ATE(1, 0)

set.seed(71281)
tmle_fit_ATE <- tmle3(tmle_spec_ATE, data, node_list, learner_list)

ATE_from_NDE_NIE <- tmle_fit_NDE$summary$tmle_est +
  tmle_fit_NIE$summary$tmle_est

test_that(
  "ATE estimate from sum of NDE and NIE TMLE estimates matches tmle3 ATE",
  expect_equal(ATE_from_NDE_NIE,
               tmle_fit_ATE$summary$tmle_est,
               tolerance = tmle_fit_ATE$summary$se)
)

################################################################################
get_sim_truth_NIE_NDE <- function(n_obs = 1e7) { # number of observations
  # compute large data set for true values
  data <- make_simulated_data(
    n_obs = n_obs
  )
  w_names <- str_subset(colnames(data), "W")
  z_names <- str_subset(colnames(data), "Z")
  W <- data[, ..w_names]
  Z <- data[, ..z_names]

  dgp <- make_dgp()
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

# simulate data and extract components for computing true parameter value
sim_truth <- get_sim_truth_NIE_NDE()

EY_A1_Z1 <- sim_truth$EY_A1_Z1
EY_A1_Z0 <- sim_truth$EY_A1_Z0
EY_A0_Z1 <- sim_truth$EY_A0_Z1
EY_A0_Z0 <- sim_truth$EY_A0_Z0

# compute true NIE via empirical substitution estimator
psi_NDE_true <- mean(EY_A1_Z0 - EY_A0_Z0)
psi_NIE_true <- mean(EY_A1_Z1 - EY_A1_Z0)

test_that("TMLE estimate of NIE for the simulation matches true value",
  expect_equal(tmle_fit_NIE$summary$tmle_est,
               psi_NIE_true,
               tolerance = tmle_fit_NIE$summary$se))

test_that("TMLE estimate of NDE for the simulation matches true value",
          expect_equal(tmle_fit_NDE$summary$tmle_est,
                       psi_NDE_true,
                       tolerance = tmle_fit_NDE$summary$se))
