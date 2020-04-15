#' Defines a TML Estimator for the Natural Indirect Effect
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
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
        targeted_likelihood,
        psi_Z_task_factory(subset_value = 0, param_type = "NIE")
      )

      lf_psi_Z1 <- tmle3::define_lf(
        tmle3::LF_derived, "psi_Z1", self$options$psi_Z_learners,
        targeted_likelihood,
        psi_Z_task_factory(subset_value = 1, param_type = "NIE")
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

###############################################################################

#' Defines a TML Estimator for the Natural Indirect Effect
#'
#' O = (W, A, Z, Y)
#' W = Covariates (possibly multivariate)
#' A = Treatment (binary or categorical)
#' Z = Mediators (binary or categorical; possibly multivariate)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param e_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'   inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'   instantiated learners from \pkg{sl3}, to be used in fitting a cleverly
#'   parameterized propensity score that conditions on the mediators, i.e.,
#'   e = P(A | Z, W).
#' @param psi_Z_learners A \code{\link[sl3]{Stack}} (or other learner class
#'   that inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or
#'   set of instantiated learners from \pkg{sl3}, to be used in a regression of
#'   a pseudo-outcome on the baseline covariates, i.e., psi_Z(W) =
#'   E[m(A = 1, Z, W) - m(A = 0, Z, W) | A = 0, W].
#' @param max_iter A \code{numeric} setting the maximum iterations allowed in
#'   the targeting step based on universal least favorable submodels.
#' @param step_size A \code{numeric} giving the step size (\code{delta_epsilon}
#'   in \code{tmle3}) to be used in the targeting step based on universal least
#'   favorable submodels.
#' @param ... Additional arguments (currently unused).
#'
#' @export
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
