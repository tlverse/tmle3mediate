library(data.table)
library(stringr)
library(dplyr)
library(hal9001)
library(sl3)
library(tmle3)

# source in helpers
testthat::source_test_helpers(here::here("tests/testthat"))
set.seed(712861)

################################################################################
# setup learners for the nuisance parameters

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

###############################################################################
# setup data and simulate to test with estimators

# create continuous Y and binary Y data and set column names for sl3 tasks
data_cont_y <- make_simulated_data(1000, binary_outcome = FALSE)
data_binary_y <- make_simulated_data(1000, binary_outcome = TRUE)

# create node list and learner list
node_list <- list(
  W = c("W_1", "W_2", "W_3"),
  A = "A",
  Z = c("Z_1", "Z_2", "Z_3"),
  Y = "Y"
)

###############################################################################

# define data (from tmle3_Spec base class)
run_NIE <- function(data, binary_outcome = FALSE, max_iter = 50) {
  if (binary_outcome) {
    learner_list <- list(
      Y = cv_hal_binary_lrnr,
      A = cv_hal_binary_lrnr
    )
  } else {
    learner_list <- list(
      Y = cv_hal_contin_lrnr,
      A = cv_hal_binary_lrnr
    )
  }

  # instantiate tmle3 spec for NIE
  tmle_spec_NIE <- tmle_NIE(
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr,
    max_iter = max_iter
  )

  tmle_task_NIE <- tmle_spec_NIE$make_tmle_task(data, node_list)

  # define likelihood (from tmle3_Spec base class)
  likelihood_init_NIE <- tmle_spec_NIE$make_initial_likelihood(
    tmle_task_NIE, learner_list
  )

  # define update method (submodel and loss function)
  updater_NIE <- tmle_spec_NIE$make_updater()
  likelihood_targeted_NIE <- Targeted_Likelihood$new(
    likelihood_init_NIE, updater_NIE
  )

  # define param
  tmle_params_NIE <- tmle_spec_NIE$make_params(
    tmle_task_NIE, likelihood_targeted_NIE
  )
  updater_NIE$tmle_params <- tmle_params_NIE

  # fit tmle update
  tmle_fit_NIE <- fit_tmle3(
    tmle_task_NIE, likelihood_targeted_NIE,
    tmle_params_NIE, updater_NIE
  )
  return(tmle_fit_NIE)
}

# run the above function for continuous and binary outcomes
tmle_fit_NIE_cont_y <- run_NIE(data = data_cont_y,
                               binary_outcome = FALSE,
                               max_iter = 10)
tmle_fit_NIE_binary_y <- run_NIE(data = data_binary_y,
                                 binary_outcome = TRUE,
                                 max_iter = 10)

# instantiate tmle3 spec for NIE outside of wrapper function
tmle_spec_NIE <- tmle_NIE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 10
)

# create learner lists and pass in to tmle3 wrapper with Spec
learner_list_cont <- list(
  Y = cv_hal_contin_lrnr,
  A = cv_hal_binary_lrnr
)
tmle_fit_NIE_cont_y_onestep <- tmle3(tmle_spec_NIE, data_cont_y,
                                     node_list, learner_list_cont)

learner_list_binary <- list(
  Y = cv_hal_binary_lrnr,
  A = cv_hal_binary_lrnr
)
tmle_fit_NIE_binary_y_onestep <- tmle3(tmle_spec_NIE, data_binary_y,
                                       node_list, learner_list_binary)

###############################################################################

run_NDE <- function(data, binary_outcome = FALSE, max_iter = 50) {
  if (binary_outcome) {
    learner_list <- list(
      Y = cv_hal_binary_lrnr,
      A = cv_hal_binary_lrnr
    )
  } else {
    learner_list <- list(
      Y = cv_hal_contin_lrnr,
      A = cv_hal_binary_lrnr
    )
  }

  # instantiate tmle3 spec for NDE
  tmle_spec_NDE <- tmle_NDE(
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr,
    max_iter = max_iter
  )

  # define data (from tmle3_Spec base class)
  tmle_task_NDE <- tmle_spec_NDE$make_tmle_task(data, node_list)

  # define likelihood (from tmle3_Spec base class)
  likelihood_init_NDE <- tmle_spec_NDE$make_initial_likelihood(
    tmle_task_NDE, learner_list
  )

  # define update method (submodel and loss function)
  updater_NDE <- tmle_spec_NDE$make_updater()
  likelihood_targeted_NDE <- Targeted_Likelihood$new(
    likelihood_init_NDE, updater_NDE
  )

  # define param
  tmle_params_NDE <- tmle_spec_NDE$make_params(
    tmle_task_NDE, likelihood_targeted_NDE
  )
  updater_NDE$tmle_params <- tmle_params_NDE

  # fit tmle update
  tmle_fit_NDE <- fit_tmle3(
    tmle_task_NDE, likelihood_targeted_NDE, tmle_params_NDE, updater_NDE
  )
  return(tmle_fit_NDE)
}

# run the above function for continuous and binary outcomes
tmle_fit_NDE_cont_y <- run_NDE(data = data_cont_y,
                               binary_outcome = FALSE,
                               max_iter = 10)
tmle_fit_NDE_binary_y <- run_NDE(data = data_binary_y,
                                 binary_outcome = TRUE,
                                 max_iter = 10)


# one-line call with faster with tmle3 wrapper
tmle_spec_NDE <- tmle_NDE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 10
)

tmle_fit_NDE_cont_y_onestep <- tmle3(
  tmle_spec_NDE, data_cont_y, node_list, learner_list_cont
)

tmle_fit_NDE_binary_y_onestep <- tmle3(
  tmle_spec_NDE, data_binary_y, node_list, learner_list_binary
)

###############################################################################
# instantiate tmle3 spec for ATE
tmle_spec_ATE <- tmle3::tmle_ATE(1, 0)

###############################################################################
# fit in the continuous case and test
tmle_fit_ATE_cont_y <- tmle3(
  tmle_spec_ATE, data_cont_y, node_list, learner_list_cont
)

ATE_from_NDE_NIE_cont_y <- tmle_fit_NDE_cont_y$summary$tmle_est +
  tmle_fit_NIE_cont_y$summary$tmle_est

test_that("TMLE for ATE equals sum of TMLEs for NDE + NIE for continuous Y", {
  expect_equal(ATE_from_NDE_NIE_cont_y, tmle_fit_ATE_cont_y$summary$tmle_est,
               tolerance = tmle_fit_ATE_cont_y$summary$se)
})

###############################################################################
# fit in the binary case and test

tmle_fit_ATE_binary_y <- tmle3(
  tmle_spec_ATE, data_binary_y, node_list, learner_list_binary
)

ATE_from_NDE_NIE_binary_y <- tmle_fit_NDE_binary_y$summary$tmle_est +
  tmle_fit_NIE_binary_y$summary$tmle_est

test_that("TMLE for ATE matches sum of TMLEs for NDE + NIE for binary Y", {
  expect_equal(ATE_from_NDE_NIE_binary_y,
               tmle_fit_ATE_binary_y$summary$tmle_est,
               tolerance = 1.96 * tmle_fit_ATE_binary_y$summary$se)
})

# NOTE: tests using truth based on simulation disabled for speed
if (FALSE) {
  #############################################################################
  # compare estimates to the truth

  # simulate data and extract components for computing true parameter value for
  # binary and continuous outcome
  sim_truth_cont_y <- get_sim_truth_NIE_NDE(
    n_obs = 1e7, binary_outcome = FALSE, EIC = FALSE
  )
  sim_truth_binary_y <- get_sim_truth_NIE_NDE(
    n_obs = 1e7, binary_outcome = TRUE, EIC = FALSE
  )

  #############################################################################
  # continuous testing
  psi_NDE_true_cont_y <- sim_truth_cont_y$NDE
  psi_NIE_true_cont_y <- sim_truth_cont_y$NIE

  test_that("TMLE estimate of NIE for simulation setting matches true value", {
    expect_equal(tmle_fit_NIE_cont_y$summary$tmle_est, psi_NIE_true_cont_y,
                 tolerance = tmle_fit_NIE_cont_y$summary$se)
  })

  test_that("TMLE estimate of NDE for simulation setting matches true value", {
    expect_equal(tmle_fit_NDE_cont_y$summary$tmle_est, psi_NDE_true_cont_y,
                 tolerance = tmle_fit_NDE_cont_y$summary$se)
  })

  #############################################################################
  # binary testing

  psi_NDE_true_binary_y <- sim_truth_binary_y$NDE
  psi_NIE_true_binary_y <- sim_truth_binary_y$NIE

  test_that("TMLE estimate of NIE for the simulation matches true value", {
    expect_equal(tmle_fit_NIE_binary_y$summary$tmle_est, psi_NIE_true_binary_y,
                 tolerance = tmle_fit_NIE_binary_y$summary$se)
  })

  test_that("TMLE estimate of NDE for the simulation matches true value", {
    expect_equal(tmle_fit_NDE_binary_y$summary$tmle_est, psi_NDE_true_binary_y,
                 tolerance = tmle_fit_NDE_binary_y$summary$se)
  })
}
