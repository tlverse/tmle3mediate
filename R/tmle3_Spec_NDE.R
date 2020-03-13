#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_NDE <- R6::R6Class(
  classname = "tmle3_Spec_NDE",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(e_learners, phi_learners,
                          max_iter = 1e4, step_size = 1e-6,
                          ...) {
      options <- list(
        e_learners = e_learners,
        phi_learners = phi_learners,
        max_iter = max_iter,
        step_size = step_size,
        ...
      )
      do.call(super$initialize, options)
    },
    make_tmle_task = function(data, node_list, ...) {
      # get variable types by guessing
      variable_types <- self$options$variable_types

      # build custom NPSEM including mediators with helper function
      npsem <- stochastic_mediation_npsem(node_list)

      # set up TMLE task based on NPSEM and return
      tmle_task <- tmle3_Task$new(data, npsem, variable_types)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # build likelihood using helper function and return
      likelihood <- stochastic_mediation_likelihood(tmle_task, learner_list)
      # note that the likelihood build here is the same as in
      # tmle3::point_tx_likelihood()
      return(likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      # add derived likelihood factors to targeted likelihood object
      # can resuse medshift e_task
      lf_e <- tmle3::define_lf(
        tmle3::LF_derived, "E", self$options$e_learners,
        targeted_likelihood, make_e_task
      )
      lf_phi <- tmle3::define_lf(
        tmle3::LF_derived, "phi", self$options$phi_learners,
        targeted_likelihood, make_NDE_phi_task
      )
      targeted_likelihood$add_factors(lf_e)
      targeted_likelihood$add_factors(lf_phi)

      # compute a tmle3 "by hand"
      tmle_params <- tmle3::define_param(Param_medshift, targeted_likelihood,
        shift_param = self$options$delta_shift
      )
      tmle_params <- list(tmle_params)
      return(tmle_params)
    },
    make_updater = function() {
      # default to ULFM approach
      updater <- tmle3_Update$new(
        one_dimensional = TRUE,
        constrain_step = TRUE,
        maxit = self$options$max_iter,
        delta_epsilon = self$options$step_size,
        cvtmle = TRUE
      )
    }
  ),
  active = list(),
  private = list()
)

#' Natural Direct Effect
#'
#' O=(W,A,Z,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Z=Mediator
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
tmle_NDE <- function(treatment_level, control_level) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_ATE$new(treatment_level, control_level)
}

################################################################################

#' Make task for derived likelihood factor phi(W) for NDE
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'  NPSEM for use in constructing elements of TML estimator.
#' @param likelihood A trained \code{Likelihood} object from \code{tmle3},
#'  constructed via the helper function \code{stochastic_mediation_likelihood}.
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom uuid UUIDgenerate
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
#
make_NDE_phi_task <- function(tmle_task, likelihood) {
  # create treatment and control tasks for intervention conditions
  treatment_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 1)
    )
  control_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 0)
    )

  # create counterfactual outcomes and construct pseudo-outcome
  m1 <- likelihood$get_likelihood(treatment_task, "Y")
  m0 <- likelihood$get_likelihood(control_task, "Y")
  m_diff <- m1 - m0

  # create regression task for pseudo-outcome and baseline covariates
  phi_data <- data.table::as.data.table(list(
    m_diff = m_diff,
    tmle_task$get_tmle_node("W")
  ))
  phi_task <- sl3::sl3_Task$new(
    data = phi_data,
    outcome = "m_diff",
    covariates = tmle_task$npsem[["W"]]$variables,
    outcome_type = "continuous"
  )

  # subset data: control observations only
  control_row_index <- tmle_task$get_tmle_node("A") == 0
  
  return(phi_task$subset_task(control_row_index))
}
