context("TML estimator for incremental propensity score interventions")

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
make_simulated_data <- function(n_obs = 1000, # no. observations
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
  z1_prob <- 1 - plogis((A^2 + W[, 1]) / (A + W[, 1]^3 + 0.5))
  Z_1 <- rbinom(n_obs, 1, prob = z1_prob)
  ## 2nd mediator (binary)
  z2_prob <- plogis((A - 1)^3 + W[, 2] / (W[, 3] + 3))
  Z_2 <- rbinom(n_obs, 1, prob = z2_prob)
  ## 3rd mediator (binary)
  z3_prob <- plogis((A - 1)^2 + 2 * W[, 1]^3 - 1 / (2 * W[, 1] + 0.5))
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

## instantiate tmle3 spec for stochastic mediation
tmle_spec <- tmle_NIE(
  e_learners = cv_hal_binary_lrnr,
  psi_Z_learners = cv_hal_contin_lrnr,
  max_iter = 100 # TODO: use default when convergence bug fixed
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
get_sim_truth_NIE <- function(n_obs = 1e7, # number of observations
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
  Y1 <- Z_1_1 + Z_2_1 - Z_3_1 +
    exp(1 + Z_3_1 / (1 + rowSums(W)^2))
  Y0 <- Z_1_0 + Z_2_0 - Z_3_0 +
    exp(0 + Z_3_0 / (1 + rowSums(W)^2))

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

  # compute p(z | A = 0, w)
  WZ_vals <- expand.grid(
    W_1 = c(0,1), W_2 = c(0,1), W_3 = c(0,1),
    Z_1 = c(0,1), Z_2 = c(0,1), Z_3 = c(0,1)
  )
  pZ_A0 <- apply(WZ_vals, MARGIN = 1, function(wz) {
    w1 <- wz["W_1"]
    w2 <- wz["W_2"]
    w3 <- wz["W_3"]
    z1 <- wz["Z_1"]
    z2 <- wz["Z_2"]
    z3 <- wz["Z_3"]
    W_subset <- data %>% filter(W_1 == w1, W_2 == w2, W_3 == w3, A == 0)
    n_W <- nrow(W_subset)
    n_Z <- nrow(W_subset %>% filter(Z_1 == z1, Z_2 == z2, Z_3 == z3))
    pZ_A0 = n_Z / n_W
  })

  pZ_A1 <- apply(WZ_vals, MARGIN = 1, function(wz) {
    w1 <- wz["W_1"]
    w2 <- wz["W_2"]
    w3 <- wz["W_3"]
    z1 <- wz["Z_1"]
    z2 <- wz["Z_2"]
    z3 <- wz["Z_3"]
    W_subset <- data %>% filter(W_1 == w1, W_2 == w2, W_3 == w3, A == 1)
    n_W <- nrow(W_subset)
    n_Z <- nrow(W_subset %>% filter(Z_1 == z1, Z_2 == z2, Z_3 == z3))
    pZ_A0 = n_Z / n_W
  })

  # compute p(W)
  pW <- data %>% group_by(W_1, W_2, W_3) %>%
    summarize(pW = n() / n_obs)

  WZ_vals <- WZ_vals %>% left_join(pW)

  # output: true values of nuisance parameters
  return(list(
    Y1 = Y1,
    Y0 = Y0,
    EY1_ZW = EYa_ZW$A1,
    EY0_ZW = EYa_ZW$A0,
    pZ_A0 = pZ_A0,
    pZ_A1 = pZ_A1,
    pW = WZ_vals$pW
  ))
}

# simulate data and extract components for computing true parameter value
sim_truth_NIE <- get_sim_truth_NIE()

Y1 <- sim_truth_NIE$Y1
Y0 <- sim_truth_NIE$Y0
EY1_ZW <- sim_truth_NIE$EY1_ZW
EY0_ZW <- sim_truth_NIE$EY0_ZW
pZ_A0 <- sim_truth_NIE$pZ_A0
pZ_A1 <- sim_truth_NIE$pZ_A1

pW <- sim_truth_NIE$pW

# compute true NIE via empirical substitution estimator
ATE <- mean(Y1 - Y0)
psi_NDE_true <- sum((EY1_ZW - EY0_ZW)*pZ_A0*pW)
psi_NIE_true <- sum(sum((EY1_ZW)*(pZ_A1 - pZ_A0))*pW)
