context("TML estimator for natural (in)direct effects")

library(data.table)
library(stringr)
library(dplyr)
library(hal9001)
library(sl3)
library(tmle3)
library(here)

# source in utils.R
testthat::source_test_helpers(here("tests/testthat"))

set.seed(7128816)

################################################################################
# setup learners for the nuisance parameters
################################################################################

# instantiate some learners
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
data <- make_simulated_data(1000)

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
  max_iter = 100
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
  max_iter = 100
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
## compare estimates to the truth

# simulate data and extract components for computing true parameter value
sim_truth <- get_sim_truth_NIE_NDE()

psi_NDE_true <- sim_truth$NDE
psi_NIE_true <- sim_truth$NIE

test_that("TMLE estimate of NIE for the simulation matches true value",
  expect_equal(tmle_fit_NIE$summary$tmle_est,
               psi_NIE_true,
               tolerance = tmle_fit_NIE$summary$se))

test_that("TMLE estimate of NDE for the simulation matches true value",
          expect_equal(tmle_fit_NDE$summary$tmle_est,
                       psi_NDE_true,
                       tolerance = tmle_fit_NDE$summary$se))
