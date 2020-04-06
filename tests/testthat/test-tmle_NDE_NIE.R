context("TML estimator for natural (in)direct effects")

library(data.table)
library(stringr)
library(dplyr)
library(hal9001)
library(sl3)
library(tmle3)
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
################################################################################
make_simulated_data <- function(n_obs = 10000, # no. observations
                                n_w = 3) { # no. baseline covariates
  # baseline covariate -- simple, binary
  W_1 <- rbinom(n_obs, 1, prob = 0.50)
  W_2 <- rbinom(n_obs, 1, prob = 0.65)
  W_3 <- rbinom(n_obs, 1, prob = 0.35)
  W <- cbind(W_1, W_2, W_3)

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = (rowSums(W) / 4 + 0.1)))

  # mediators to affect the outcome
  ## 1st mediator (binary)
  z1_prob <- 1 - plogis((A + W[, 1]) / (A + W[, 1]^3 + 0.5))
  Z_1 <- rbinom(n_obs, 1, prob = z1_prob)
  ## 2nd mediator (binary)
  z2_prob <- plogis((A - 1) + W[, 2] / (W[, 3] + 3))
  Z_2 <- rbinom(n_obs, 1, prob = z2_prob)
  ## 3rd mediator (binary)
  z3_prob <- plogis((A - 1) + 2 * W[, 1]^3 - 1 / (2 * W[, 1] + 0.5))
  Z_3 <- rbinom(n_obs, 1, prob = z3_prob)
  ## build matrix of mediators
  Z <- cbind(Z_1, Z_2, Z_3)

  # create outcome as a function of A, Z, W + white noise
  Y <- Z_1 + Z_2 - Z_3 + exp(A + Z_3 / (1 + rowSums(W)^2)) +
    rnorm(n_obs, mean = 0, sd = 0.5)

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c(
    "Y", paste("Z", 1:3, sep = "_"), "A",
    paste("W", seq_len(dim(W)[2]), sep = "_")
  ))
  return(data)
}

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

## define data (from tmle3_Spec base class)
tmle_task_NIE <- tmle_spec_NIE$make_tmle_task_NIE(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init_NIE <- tmle_spec_NIE$make_initial_likelihood(tmle_task_NIE, learner_list)

## define update method (submodel and loss function)
updater_NIE <- tmle_spec_NIE$make_updater_NIE()
likelihood_targeted_NIE <- Targeted_Likelihood$new(likelihood_init_NIE, updater_NIE)

## define param
tmle_params_NIE <- tmle_spec_NIE$make_params(tmle_task_NIE, likelihood_targeted_NIE)
updater_NIE$tmle_params_NIE <- tmle_params_NIE

## fit tmle update
tmle_fit_NIE <- fit_tmle3(tmle_task_NIE, likelihood_targeted_NIE, tmle_params_NIE, updater_NIE)
tmle_fit_NIE

## one-line call with faster with tmle3 wrapper
set.seed(71281)
tmle_fit_NIE <- tmle3(tmle_spec_NIE, data, node_list, learner_list)
tmle_fit_NIE

################################################################################

## instantiate tmle3 spec for NDE
tmle_spec_NDE <- tmle_NDE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 100 # TODO: use default when convergence bug fixed
)

## define data (from tmle3_Spec base class)
tmle_task_NDE <- tmle_spec_NDE$make_tmle_task_NDE(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init_NDE <- tmle_spec_NDE$make_initial_likelihood(tmle_task_NDE, learner_list)

## define update method (submodel and loss function)
updater_NDE <- tmle_spec_NDE$make_updater_NDE()
likelihood_targeted_NDE <- Targeted_Likelihood$new(likelihood_init_NDE, updater_NDE)

## define param
tmle_params_NDE <- tmle_spec_NDE$make_params(tmle_task_NDE, likelihood_targeted_NDE)
updater_NDE$tmle_params_NDE <- tmle_params_NDE

## fit tmle update
tmle_fit_NDE <- fit_tmle3(tmle_task_NDE, likelihood_targeted_NDE, tmle_params_NDE, updater_NDE)

## one-line call with faster with tmle3 wrapper
set.seed(71281)
tmle_fit_NDE <- tmle3(tmle_spec_NDE, data, node_list, learner_list)
tmle_fit_NDE

################################################################################
## instantiate tmle3 spec for ATE
tmle_spec_ATE <- tmle3::tmle_ATE(1, 0)

set.seed(71281)
tmle_fit_ATE <- tmle3(tmle_spec_ATE, data, node_list, learner_list)
tmle_fit_ATE

################################################################################
get_sim_truth_NIE_NDE <- function(n_obs = 1e7, # number of observations
                                  n_w = 3) { # number of baseline covariates
  # compute large data set for true values
  data <- make_simulated_data(
    n_obs = n_obs,
    n_w = n_w
  )
  w_names <- str_subset(colnames(data), "W")
  z_names <- str_subset(colnames(data), "Z")
  W <- data[, ..w_names]
  Z <- data[, ..z_names]

  # Z_1 counterfactuals
  z1_prob_1 <- 1 - plogis((1 + W$W_1) / (1 + W$W_1^3 + 0.5))
  z1_prob_0 <- 1 - plogis((0 + W$W_1) / (0 + W$W_1^3 + 0.5))
  Z_1_1 <- rbinom(n_obs, 1, prob = z1_prob_1)
  Z_1_0 <- rbinom(n_obs, 1, prob = z1_prob_0)

  # Z_2 counterfactuals
  z2_prob_1 <- plogis((1 - 1)^3 + W$W_2 / (W$W_3 + 3))
  z2_prob_0 <- plogis((0 - 1)^3 + W$W_2 / (W$W_3 + 3))
  Z_2_1 <- rbinom(n_obs, 1, prob = z2_prob_1)
  Z_2_0 <- rbinom(n_obs, 1, prob = z2_prob_0)

  ## Z_3 counterfactuals
  z3_prob_1 <- plogis((1 - 1)^2 + 2 * W$W_1^3 - 1 / (2 * W$W_1 + 0.5))
  z3_prob_0 <- plogis((0 - 1)^2 + 2 * W$W_1^3 - 1 / (2 * W$W_1 + 0.5))
  Z_3_1 <- rbinom(n_obs, 1, prob = z3_prob_1)
  Z_3_0 <- rbinom(n_obs, 1, prob = z3_prob_0)

  # compute counterfactuals Y(1) and Y(0)
  EY_A1_Z1 <- Z_1_1 + Z_2_1 - Z_3_1 +
    exp(1 + Z_3_1 / (1 + rowSums(W)^2))

  EY_A0_Z0 <- Z_1_0 + Z_2_0 - Z_3_0 +
    exp(0 + Z_3_0 / (1 + rowSums(W)^2))

  EY_A1_Z0 <- Z_1_0 + Z_2_0 - Z_3_0 +
    exp(1 + Z_3_0 / (1 + rowSums(W)^2))

  EY_A0_Z1 <- Z_1_1 + Z_2_1 - Z_3_1 +
    exp(0 + Z_3_1 / (1 + rowSums(W)^2))

  # compute TRUE M under counterfactual regimes
  m_Ais1 <- Z$Z_1 + Z$Z_2 - Z$Z_3 + exp(1 + Z$Z_3 / (1 + rowSums(W)^2))
  m_Ais0 <- Z$Z_1 + Z$Z_2 - Z$Z_3 + exp(0 + Z$Z_3 / (1 + rowSums(W)^2))

  # compute E(Y | A = a, W, Z) for A = 0,1 and all levels of (W,Z)
  EYa_ZW <- data %>%
    mutate(
      m_Ais1 = m_Ais1,
      m_Ais0 = m_Ais0
    ) %>%
    group_by(W_1, W_2, W_3, Z_1, Z_2, Z_3) %>%
    summarize(A1 = mean(m_Ais1), A0 = mean(m_Ais0))

  # output: true values of nuisance parameters
  return(list(
    EY_A1_Z1 = EY_A1_Z1,
    EY_A1_Z0 = EY_A1_Z0,
    EY_A0_Z1 = EY_A0_Z1,
    EY_A0_Z0 = EY_A0_Z0,
    EY1_ZW = EYa_ZW$A1,
    EY0_ZW = EYa_ZW$A0
  ))
}

# simulate data and extract components for computing true parameter value
sim_truth <- get_sim_truth_NIE_NDE()

EY_A1_Z1 <- sim_truth$EY_A1_Z1
EY_A1_Z0 <- sim_truth$EY_A1_Z0
EY_A0_Z1 <- sim_truth$EY_A0_Z1
EY_A0_Z0 <- sim_truth$EY_A0_Z0

EY1_ZW <- sim_truth$EY1_ZW
EY0_ZW <- sim_truth$EY0_ZW

# compute true NIE via empirical substitution estimator
psi_NDE_true <- mean(EY_A1_Z0 - EY_A0_Z0)
psi_NDE_true
psi_NIE_true <- mean(EY_A1_Z1 - EY_A1_Z0)
psi_NIE_true
ATE <- psi_NDE_true + psi_NIE_true
ATE
ATE_2 <- mean(EY1_ZW - EY0_ZW)
ATE_2

test_that("TMLE estimate of NIE for the simulation matches expected", 
          expect_equal(tmle_fit_NIE$summary$tmle_est, psi_NIE_true, tolerance = tmle_fit_NIE$summary$se))

test_that("TMLE estimate of NDE for the simulation matches expected", 
          expect_equal(tmle_fit_NDE$summary$tmle_est, psi_NDE_true, tolerance = tmle_fit_NDE$summary$se))

test_that("TMLE estimate of NDE for the simulation matches expected", 
          expect_equal(tmle_fit_NDE$summary$tmle_est + tmle_fit_NIE$summary$tmle_est, 
                       ATE, 
                       tolerance = 0.01))


