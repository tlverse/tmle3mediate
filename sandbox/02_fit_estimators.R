est_nde_nie <- function(data, m_learners, g_learners, e_learners,
                        psi_Z_learners) {
  node_list <- list(
    W = c("W_1", "W_2", "W_3"),
    A = "A",
    Z = c("Z_1", "Z_2", "Z_3"),
    Y = "Y"
  )
  learner_list <- list(
    Y = m_learners,
    A = g_learners
  )
  tmle_spec_NIE <- tmle_NIE(
    e_learners = e_learners,
    psi_Z_learners = psi_Z_learners,
    max_iter = 100
  )
  NIE_est <- tmle3(tmle_spec_NIE, data, node_list, learner_list)

  tmle_spec_NDE <- tmle_NDE(
    e_learners = e_learners,
    psi_Z_learners = psi_Z_learners,
    max_iter = 100
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
  fglm_lrnr <- Lrnr_cv$new(Lrnr_glm_fast$new(), full_fit = TRUE)

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

  # compute TMLE under different misspecification settings
  ## 1) all nuisance functions correctly specified

  correct <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = cv_hal_binary_lrnr,
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr
  )

  ###################################################################

  ## 2) g and e misspecified
  ## Lemma 1. (i) of 2012 Zheng and van der Laan satisfied

  mis_i <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = mean_lrnr,
    e_learners = mean_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr
  )

  ##############################################################

  ## 3) e and psi_Z misspecified
  ## Lemma 1. (ii) of 2012 Zheng and van der Laan satisfied

  mis_ii <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = cv_hal_contin_lrnr,
    e_learners = mean_lrnr,
    psi_Z_learners = mean_lrnr
  )

  ###################################################################

  ## 4) m and psi_Z misspecified
  ## Lemma 1. (iii) of 2012 Zheng and van der Laan satisfied

  mis_iii <- est_nde_nie(
    data,
    m_learners = fglm_lrnr,
    g_learners = cv_hal_binary_lrnr,
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = mean_lrnr
  )

  ###################################################################

  ## 5) e misspecified
  ## Lemma 1. (i) and (ii) of 2012 Zheng and van der Laan satisfied

  mis_e <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = cv_hal_binary_lrnr,
    e_learners = mean_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr
  )

  ###################################################################

  ## 6) g misspecified
  ## Lemma 1. (i) of 2012 Zheng and van der Laan satisfied

  mis_g <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = mean_lrnr,
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = cv_hal_contin_lrnr
  )

  ###################################################################

  ## 7) psi_Z misspecified
  ## Lemma 1. (ii) and (iii) of 2012 Zheng and van der Laan satisfied

  mis_psi_Z <- est_nde_nie(
    data,
    m_learners = cv_hal_contin_lrnr,
    g_learners = cv_hal_binary_lrnr,
    e_learners = cv_hal_binary_lrnr,
    psi_Z_learners = mean_lrnr
  )

  ###################################################################

  ## 8) m misspecified
  ## Lemma 1. (iii) of 2012 Zheng and van der Laan satisfied

  mis_m <- est_nde_nie(
    data,
    m_learners = fglm_lrnr,
    g_learners = cv_hal_binary_lrnr,
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
                    mis_i_NIE = summary_with_steps(mis_i$nie),
                    mis_i_NDE = summary_with_steps(mis_i$nde),
                    mis_ii_NIE = summary_with_steps(mis_ii$nie),
                    mis_ii_NDE = summary_with_steps(mis_ii$nde),
                    mis_iii_NIE = summary_with_steps(mis_iii$nie),
                    mis_iii_NDE = summary_with_steps(mis_iii$nde),
                    mis_e_NIE = summary_with_steps(mis_e$nie),
                    mis_e_NDE = summary_with_steps(mis_e$nde),
                    mis_g_NIE = summary_with_steps(mis_g$nie),
                    mis_g_NDE = summary_with_steps(mis_g$nde),
                    mis_psi_Z_NIE = summary_with_steps(mis_psi_Z$nie),
                    mis_psi_Z_NDE = summary_with_steps(mis_psi_Z$nde),
                    mis_m_NIE = summary_with_steps(mis_m$nie),
                    mis_m_NDE = summary_with_steps(mis_m$nde))

  sim_out <- bind_rows(estimates, .id = "sim_type")
  return(sim_out)
}
