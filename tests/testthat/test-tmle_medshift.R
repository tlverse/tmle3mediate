context("TML estimator for incremental propensity score interventions")

library(data.table)
library(stringr)
library(hal9001)
library(sl3)
library(tmle3)
set.seed(7128816)

# delta used for stochastic intervention
delta_ipsi <- 0.5

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

# get data and column names for sl3 tasks (for convenience)
data <- make_simulated_data(n_obs = 10000, binary_outcome = FALSE)
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

## instantiate tmle3 spec for stochastic mediation
tmle_spec <- tmle3mediate::tmle_medshift(
  delta = delta_ipsi,
  e_learners = cv_hal_binary_lrnr,
  phi_learners = cv_hal_contin_lrnr
)

## define data (from tmle3_Spec base class)
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

## define update method (submodel and loss function)
updater <- tmle_spec$make_updater()
likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

## define param
tmle_params <- tmle_spec$make_params(tmle_task, likelihood_targeted)
updater$tmle_params <- tmle_params

## fit tmle update
tmle_fit <- fit_tmle3(tmle_task, likelihood_targeted, tmle_params, updater)

## one-line call with faster with tmle3 wrapper
set.seed(71281)
tmle_fit <- tmle3(tmle_spec, data, node_list, learner_list)
tmle_fit



################################################################################
get_sim_truth_medshift <- function(n_obs = 1e7, # number of observations
                                   delta = 0.5) { # value of shift parameter

  # compute large data set for true values
  data <- make_simulated_data(
    n_obs = n_obs,
    binary_outcome = FALSE
  )

  w_names <- str_subset(colnames(data), "W")
  z_names <- str_subset(colnames(data), "Z")
  W <- data[, ..w_names]
  Z <- data[, ..z_names]
  Y <- data$Y

  dgp <- make_dgp()

  # compute TRUE G under counterfactual regimes
  g_Ais1 <- dgp$g_mech(W)
  g_Ais0 <- 1 - g_Ais1

  # compute TRUE SHIFTED G under counterfactual regimes
  g_shifted_Ais1 <- (delta * g_Ais1) / (delta * g_Ais1 + g_Ais0)
  g_shifted_Ais0 <- 1 - g_shifted_Ais1

  # compute TRUE M under counterfactual regimes
  m_Ais1 <- dgp$m_mech_cont(W, 1, Z, eps_sd = 0)
  m_Ais0 <- dgp$m_mech_cont(W, 0, Z, eps_sd = 0)

  # output: true values of nuisance parameters
  return(list(
    g_obs_true = data.table(
      A1 = g_Ais1,
      A0 = g_Ais0
    ),
    g_shifted_true = data.table(
      A1 = g_shifted_Ais1,
      A0 = g_shifted_Ais0
    ),
    m_true = data.table(
      A1 = m_Ais1,
      A0 = m_Ais0
    ),
    EY_true = mean(Y)
  ))
}

# simulate data and extract components for computing true parameter value
sim_truth <- get_sim_truth_medshift()
m_A1 <- sim_truth$m_true$A1
m_A0 <- sim_truth$m_true$A0
g_shifted_A1 <- sim_truth$g_shifted_true$A1
g_shifted_A0 <- sim_truth$g_shifted_true$A0
EY <- sim_truth$EY_true

# compute true parameter value based on the substitution estimator
true_param <- mean(m_A1 * g_shifted_A1) + mean(m_A0 * g_shifted_A0)

test_that(
  "TMLE estimate of medshift param for the simulation matches true value",
  expect_equal(tmle_fit$summary$tmle_est,
    true_param,
    tolerance = tmle_fit$summary$se
  )
)
