#' TML Estimator for the Natural Direct Effect
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
tmle3_Spec_NDE <- R6::R6Class(
  classname = "tmle3_Spec_NDE",
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
      lf_psi_Z <- tmle3::define_lf(
        tmle3::LF_derived, "psi_Z", self$options$psi_Z_learners,
        targeted_likelihood, make_NDE_psi_Z_task
      )
      targeted_likelihood$add_factors(list(lf_e, lf_psi_Z))

      # create param
      tmle_params <- Param_NDE$new(targeted_likelihood)
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

#' TML Estimator for the Natural Direct Effect
#'
#' O = (W, A, Z, Y)
#' W = Covariates (possibly multivariate)
#' A = Treatment (binary or categorical)
#' Z = Mediators (binary or categorical; possibly multivariate)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param e_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, to be used in fitting a cleverly
#'  parameterized propensity score that conditions on the mediators, i.e.,
#'  \eqn{e = P(A \mid Z, W)}.
#' @param psi_Z_learners A \code{\link[sl3]{Stack}} (or other learner class
#'  that inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or
#'  set of instantiated learners from \pkg{sl3}, to be used in a regression of
#'  a pseudo-outcome on the baseline covariates, i.e.,
#'  \eqn{psi_Z(W) = E[m(A = 1, Z, W) - m(A = 0, Z, W) \mid A = 0, W]}.
#' @param max_iter A \code{numeric} setting the maximum iterations allowed in
#'  the targeting step based on universal least favorable submodels.
#' @param step_size A \code{numeric} giving the step size (\code{delta_epsilon}
#'  in \code{\link[tmle3]{tmle3}}) to be used in the targeting step based on
#'  universal least favorable submodels.
#' @param ... Additional arguments (currently unused).
#'
#' @export
tmle_NDE <- function(e_learners, psi_Z_learners,
                     max_iter = 1e4, step_size = 1e-6,
                     ...) {
  # this is a factory function
  tmle3_Spec_NDE$new(
    e_learners, psi_Z_learners,
    max_iter, step_size,
    ...
  )
}

###############################################################################

#' Make Task for Derived Likelihood Factor psi_Z(W) for Natural Direct Effect
#'
#' @param tmle_task A \code{[tmle3]{tmle3_Task}} specifying the data and NPSEM
#'  for use in constructing components required for TML estimation.
#' @param likelihood A trained \code{[tmle3]{Likelihood}}, constructed via the
#'  \code{\link{stochastic_mediation_likelihood}} helper.
#'
#' @importFrom data.table as.data.table data.table setnames
#' @importFrom uuid UUIDgenerate
#' @importFrom sl3 sl3_Task
#'
#' @name make_NDE_psi_Z_task
#'
#' @keywords internal
make_NDE_psi_Z_task <- function(tmle_task, likelihood) {
  mediator_task_fun_factory(param_type = "NDE", subset_value = 0)(
    tmle_task, likelihood
  )
}
