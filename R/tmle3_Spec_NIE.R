#' Defines a TML Estimator for Outcome under Joint Static Intervention on
#' Treatment and Mediator
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
#
tmle3_Spec_NIE <- R6::R6Class(
  classname = "tmle3_Spec_NIE",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(e_learners, psi_Z_learners,
                          max_iter = 1e4, step_size = 1e-6,
                          ...) {
      options <- list(
        e_learners = e_learners,
        psi_Z_learners = psi_Z_learners,
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
      # note that the likelihood built here is the same as in
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
      
      lf_psi_Z0 <- tmle3::define_lf(
        tmle3::LF_derived, "psi_Z0", self$options$psi_Z_learners,
        targeted_likelihood, make_NIE_psi_Z0_task
      )
      
      lf_psi_Z1 <- tmle3::define_lf(
        tmle3::LF_derived, "psi_Z1", self$options$psi_Z_learners,
        targeted_likelihood, make_NIE_psi_Z1_task
      )
      
      targeted_likelihood$add_factors(lf_e)
      targeted_likelihood$add_factors(lf_psi_Z0)
      targeted_likelihood$add_factors(lf_psi_Z1)
      
      # create param
      tmle_params <- Param_NIE$new(targeted_likelihood)
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


################################################################################

#' Outcome under Joint Static Intervention on Treatment and Mediator
#'
#' O = (W, A, Z, Y)
#' W = Covariates (possibly multivariate)
#' A = Treatment (binary or categorical)
#' Z = Mediators (binary or categorical; possibly multivariate)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param e_learners A \code{Stack} object, or other learner class (inheriting
#'   from \code{Lrnr_base}), containing a single or set of instantiated learners
#'   from the \code{sl3} package, to be used in fitting a cleverly parameterized
#'   propensity score that includes the mediators, i.e., e = P(A | Z, W).
#' @param psi_Z_learners A \code{Stack} object, or other learner class
#'   (inheriting from \code{Lrnr_base}), containing a single or set of
#'   instantiated learners from the \code{sl3} package, to be used in fitting a
#'   reduced regression useful for computing the psi_Z nuisance parameter, i.e.,
#'   psi_Z(W) = E[m(A = 1, Z, W) - m(A = 0, Z, W) | A = 0, W].
#' @param max_iter A \code{numeric} setting the total number of iterations to be
#'   used in the targeted procedure based on universal least favorable
#'   submodels.
#' @param step_size A \code{numeric} giving the step size (\code{delta_epsilon}
#'   in \code{tmle3}) to be used in the targeted procedure based on universal
#'   least favorable submodels.
#' @param ... Additional arguments (currently unused).
#'
#' @export
#
tmle_NIE <- function(e_learners, psi_Z_learners,
                     max_iter = 1e4, step_size = 1e-6,
                     ...) {
  # this is a factory function
  tmle3_Spec_NIE$new(
    e_learners, psi_Z_learners,
    max_iter, step_size,
    ...
  )
}

################################################################################

#' Make task for derived likelihood factor psi_Z(W) for NIE
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
make_NIE_psi_Z0_task <- function(tmle_task, likelihood) {
  # create treatment and control tasks for intervention conditions
  # TODO: remove node name hard-coding
  treatment_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 1)
    )
  
  # create counterfactual outcomes and construct pseudo-outcome
  m1 <- likelihood$get_likelihood(treatment_task, "Y")
  
  # create regression task for pseudo-outcome and baseline covariates
  psi_Z_data <- data.table::as.data.table(list(
    m1 = m1,
    tmle_task$get_tmle_node("W")
  ))
  
  psi_Z_task <- sl3::sl3_Task$new(
    data = psi_Z_data,
    outcome = "m1",
    covariates = tmle_task$npsem[["W"]]$variables,
    outcome_type = "continuous"
  )
  
  # subset data: control observations only
  control_row_index <- tmle_task$get_tmle_node("A") == 0

  return(psi_Z_task$subset_task(control_row_index))
}

make_NIE_psi_Z1_task <- function(tmle_task, likelihood) {
  # create treatment and control tasks for intervention conditions
  # TODO: remove node name hard-coding
  treatment_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 1)
    )
  
  # create counterfactual outcomes and construct pseudo-outcome
  m1 <- likelihood$get_likelihood(treatment_task, "Y")
  
  # create regression task for pseudo-outcome and baseline covariates
  psi_Z_data <- data.table::as.data.table(list(
    m1 = m1,
    tmle_task$get_tmle_node("W")
  ))
  
  psi_Z_task <- sl3::sl3_Task$new(
    data = psi_Z_data,
    outcome = "m1",
    covariates = tmle_task$npsem[["W"]]$variables,
    outcome_type = "continuous"
  )
  
  # subset data: control observations only
  treatment_row_index <- tmle_task$get_tmle_node("A") == 1
  
  return(psi_Z_task$subset_task(treatment_row_index))
}
