est_nde_nie <- function(data, m_learners, g_learners, e_learners,
                        psi_Z_learners) {
  node_list <- list(
    W = c("W1", "W2", "W3"),
    A = "A",
    Z = c("Z1", "Z2", "Z3"),
    Y = "Y"
  )
  learner_list <- list(
    Y = m_learners,
    A = g_learners
  )
  tmle_spec_NIE <- tmle_NIE(
    e_learners = e_learners,
    psi_Z_learners = psi_Z_learners,
    max_iter = 50
  )
  NIE_est <- tmle3(tmle_spec_NIE, data, node_list, learner_list)

  tmle_spec_NDE <- tmle_NDE(
    e_learners = e_learners,
    psi_Z_learners = psi_Z_learners,
    max_iter = 50
  )
  NDE_est <- tmle3(tmle_spec_NDE, data, node_list, learner_list)

  return(list(nie = NIE_est, nde = NDE_est))
}

###############################################################################
# fit estimators
###############################################################################
fit_estimators <- function(data, cv_folds = 5) {
  ## setup learners

  # misspecified learner
  mean_lrnr <- Lrnr_cv$new(Lrnr_mean$new(), full_fit = TRUE)

  # hal9001
  hal_contin_lrnr <- Lrnr_hal9001$new(n_folds = cv_folds,
                                      fit_type = "glmnet",
                                      yolo = FALSE)
  hal_binary_lrnr <- Lrnr_hal9001$new(n_folds = cv_folds,
                                      fit_type = "glmnet",
                                      family = "binomial",
                                      yolo = FALSE)
  cv_hal_contin_lrnr <- Lrnr_cv$new(hal_contin_lrnr, full_fit = TRUE)
  cv_hal_binary_lrnr <- Lrnr_cv$new(hal_binary_lrnr, full_fit = TRUE)

  # TODO: check misspecification in Wenjing's paper again
  # TODO: check three robustness conditions in sim
  # TODO: sometimes have to force hal learners to go deep enough into sequence of lambdas
  #       using binomial family, glmnet ignores the lambda.min.ratio argument

  # compute TMLE under different misspecification settings
  ## 1) all nuisance functions correctly specified

  correct <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = cv_hal_binary_lrnr,
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr
  )

  ##############################################################

  ## 2) e misspecified

  mis_e <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = cv_hal_binary_lrnr,
    e_learners = mean_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr
  )

  ###################################################################

  ## 3) m misspecified

  mis_m <- est_nde_nie(
    data,
    m_learners = mean_lrnr,
    g_learners = cv_hal_binary_lrnr,
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = mean_lrnr
  )

  ###################################################################

  ## 4) g misspecified

  mis_g <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = mean_lrnr,
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr
  )

  # helper to add steps to the tmle summary
  summary_with_steps <- function(tmle_out) {
    summary_list <- tmle_out$summary
    summary_list$steps <- tmle_out$steps
    summary_list
  }

  # bundle estimators in list
  estimates <- list(corr_NIE = summary_with_steps(correct$nie),
                    corr_NDE = summary_with_steps(correct$nde),
                    mis_e_NIE = summary_with_steps(mis_e$nie),
                    mis_e_NDE = summary_with_steps(mis_e$nde),
                    mis_m_NIE = summary_with_steps(mis_m$nie),
                    mis_m_NDE = summary_with_steps(mis_m$nde),
                    mis_g_NIE = summary_with_steps(mis_g$nie),
                    mis_g_NDE = summary_with_steps(mis_g$nde))

  sim_out <- bind_rows(estimates, .id = "sim_type")
  return(sim_out)
}
