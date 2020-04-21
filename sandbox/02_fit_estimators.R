###############################################################################
# fit estimators
###############################################################################
fit_estimators <- function(data, cv_folds = 5) {
  ## NOTE: need to force `family = "gaussian"` to use `lambda.min.ratio`
  # hal_lrnr <- Lrnr_hal9001$new(max_degree = NULL,
  #                              n_folds = 5,
  #                              fit_type = "glmnet",
  #                              use_min = TRUE,
  #                              type.measure = "deviance",
  #                              standardize = FALSE,
  #                              family = "gaussian",
  #                              lambda.min.ratio = 1 / nrow(data),
  #                              nlambda = 500,
  #                              yolo = FALSE)

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

  # meta-learner to ensure predicted probabilities do not go outside [0,1]
  logistic_metalearner <- make_learner(Lrnr_solnp,
                                       metalearner_logistic_binomial,
                                       loss_loglik_binomial)

  # set nuisance regression learners based on ID's successful simulations
  sl <- Lrnr_sl$new(learners = list(hal_lrnr),
                    metalearner = Lrnr_nnls$new())
                    #metalearner = logistic_metalearner)
  node_list <- list(
    W = c("W1", "W2", "W3"),
    A = "A",
    Z = c("Z1", "Z2", "Z3"),
    Y = "Y"
  )
  learner_list <- list(
    Y = cv_hal_contin_lrnr,
    A = cv_hal_binary_lrnr
  )


  # compute TMLE under different misspecification settings
  ## 1) all nuisance functions correctly specified

  tmle_spec_NIE <- tmle_NIE(
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr,
    max_iter = 100 # TODO: use default when convergence bug fixed
  )

  NIE_est_corr <- tmle3(tmle_spec_NIE, data, node_list, learner_list)
  NIE_est_corr <- NIE_est_corr$summary

  tmle_spec_NDE <- tmle_NDE(
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr,
    max_iter = 100 # TODO: use default when convergence bug fixed
  )

  NDE_est_corr <- tmle3(tmle_spec_NDE, data, node_list, learner_list)
  NDE_est_corr <- NDE_est_corr$summary

  ##############################################################

  ## 2) e misspecified

  tmle_spec_NIE <- tmle_NIE(
    e_learners = mean_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr,
    max_iter = 100 # TODO: use default when convergence bug fixed
  )

  NIE_est_mis_e <- tmle3(tmle_spec_NIE, data, node_list, learner_list)
  NIE_est_mis_e <- NIE_est_mis_e$summary

  tmle_spec_NDE <- tmle_NDE(
    e_learners = mean_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr,
    max_iter = 100 # TODO: use default when convergence bug fixed
  )

  NDE_est_mis_e <- tmle3(tmle_spec_NDE, data, node_list, learner_list)
  NDE_est_mis_e <- NDE_est_mis_e$summary


  ###################################################################

  ## 3) z misspecified

  tmle_spec_NIE <- tmle_NIE(
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = mean_lrnr,
    max_iter = 100 # TODO: use default when convergence bug fixed
  )

  NIE_est_mis_z <- tmle3(tmle_spec_NIE, data, node_list, learner_list)
  NIE_est_mis_z <- NIE_est_mis_z$summary

  tmle_spec_NDE <- tmle_NDE(
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = mean_lrnr,
    max_iter = 100 # TODO: use default when convergence bug fixed
  )

  NDE_est_mis_z <- tmle3(tmle_spec_NDE, data, node_list, learner_list)
  NDE_est_mis_z <- NDE_est_mis_z$summary


  # bundle estimators in list
  estimates <- list(corr_NIE = NIE_est_corr, corr_NDE = NDE_est_corr,
                    mis_e_NIE = NIE_est_mis_e, mis_e_NDE = NDE_est_mis_e,
                    mis_z_NIE = NIE_est_mis_z, mis_z_NDE = NDE_est_mis_z)

  sim_out <- bind_rows(estimates, .id = "sim_type")
  return(sim_out)
}
