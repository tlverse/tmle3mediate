###############################################################################
# fit estimators
###############################################################################
fit_estimators <- function(data, cv_folds = 5) {
  ## NOTE: need to force `family = "gaussian"` to use `lambda.min.ratio`
  hal_lrnr <- Lrnr_hal9001$new(max_degree = NULL,
                               n_folds = 5,
                               fit_type = "glmnet",
                               use_min = TRUE,
                               type.measure = "deviance",
                               standardize = FALSE,
                               family = "gaussian",
                               lambda.min.ratio = 1 / nrow(data),
                               nlambda = 500,
                               yolo = FALSE)
  mean_lrnr <- Lrnr_mean$new()

  # meta-learner to ensure predicted probabilities do not go outside [0,1]
  logistic_metalearner <- make_learner(Lrnr_solnp,
                                       metalearner_logistic_binomial,
                                       loss_loglik_binomial)

  # set nuisance regression learners based on ID's successful simulations
  sl <- Lrnr_sl$new(learners = list(myglmnet_lrnr, hal_lrnr),
                    metalearner = Lrnr_nnls$new())
                    #metalearner = logistic_metalearner)

  # compute TMLE under different misspecification settings
  ## 1) all nuisance functions correctly specified
  est_corr <- estimators(data = data,
                         g_stack = sl,
                         e_stack = sl,
                         m_stack = sl,
                         cv_folds = cv_folds)

  ## 2) g misspecified
  est_misg <- estimators(data = data,
                         g_stack = mean_lrnr,
                         e_stack = sl,
                         m_stack = sl,
                         cv_folds = cv_folds)

  ## 3) e misspecified
  est_mise <- estimators(data = data,
                         g_stack = sl,
                         e_stack = mean_lrnr,
                         m_stack = sl,
                         cv_folds = cv_folds)

  ## 4) m misspecified
  est_mism <- estimators(data = data,
                         g_stack = sl,
                         e_stack = sl,
                         m_stack = mean_lrnr,
                         cv_folds = cv_folds)

  # bundle estimators in list
  estimates <- list(corr = est_corr, misg = est_misg, mise = est_mise,
                    mism = est_mism)
  sim_out <- bind_rows(estimates, .id = "sim_type")
  return(sim_out)
}
