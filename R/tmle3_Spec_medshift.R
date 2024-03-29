#' TML Estimator for the Population Intervention (In)direct Effects
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
tmle3_Spec_medshift <- R6::R6Class(
  classname = "tmle3_Spec_medshift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(shift_type = "ipsi", delta,
                          e_learners, phi_learners,
                          max_iter = 1e4, step_size = 1e-6,
                          ...) {
      options <- list(
        shift_type = shift_type,
        delta_shift = delta,
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
      return(likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      # add derived likelihood factors to targeted likelihood object
      lf_e <- tmle3::define_lf(
        tmle3::LF_derived, "E", self$options$e_learners,
        targeted_likelihood, make_e_task
      )
      lf_phi <- tmle3::define_lf(
        tmle3::LF_derived, "phi", self$options$phi_learners,
        targeted_likelihood, make_phi_task
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

###############################################################################

#' TML Estimator for the Population Intervention (In)direct Effects
#'
#' O = (W, A, Z, Y)
#' W = Covariates (possibly multivariate)
#' A = Treatment (binary or categorical)
#' Z = Mediators (binary or categorical; possibly multivariate)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param shift_type A \code{character} defining the type of shift to be
#'  applied to the exposure -- an incremental propensity score intervention.
#' @param delta A \code{numeric}, specifying the magnitude of the shift.
#' @param e_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, to be used in fitting a cleverly
#'  parameterized propensity score that conditions on the mediators, i.e.,
#'  \eqn{e = P(A \mid Z, W)}.
#' @param phi_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, to be used in a regression of a
#'  pseudo-outcome on the baseline covariates, i.e.,
#'  \eqn{phi(W) = E[m(A = 1, Z, W) - m(A = 0, Z, W) | W)]}.
#' @param max_iter A \code{numeric} setting the maximum iterations allowed in
#'  the targeting step based on universal least favorable submodels.
#' @param step_size A \code{numeric} giving the step size (\code{delta_epsilon}
#'  in \code{\link[tmle3]{tmle3}}) to be used in the targeting step based on
#'  universal least favorable submodels.
#' @param ... Additional arguments (currently unused).
#'
#' @export
tmle_medshift <- function(shift_type = "ipsi",
                          delta, e_learners, phi_learners,
                          max_iter = 1e4, step_size = 1e-6,
                          ...) {
  # this is a factory function
  tmle3_Spec_medshift$new(
    shift_type, delta,
    e_learners, phi_learners,
    max_iter, step_size,
    ...
  )
}

###############################################################################

#' Stochastic Mediation NPSEM
#'
#' @param node_list A \code{list} object specifying the different nodes in the
#'  nonparametric structural equation model (NPSEM).
#' @param variable_types Used to define how variables are handled. Optional.
#'
#' @importFrom tmle3 define_node
#'
#' @keywords internal
stochastic_mediation_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task
  npsem <- list(
    tmle3::define_node("W", node_list$W, variable_type = variable_types$W),
    tmle3::define_node("A", node_list$A, c("W"),
      variable_type = variable_types$A
    ),
    tmle3::define_node("Z", node_list$Z, c("A", "W"),
      variable_type = variable_types$Z
    ),
    tmle3::define_node("Y", node_list$Y, c("Z", "A", "W"),
      variable_type = variable_types$Y, scale = TRUE
    )
  )
  return(npsem)
}

###############################################################################

#' Mediation Likelihood Factors for Population Intervention (In)Direct Effects
#'
#' @param tmle_task A \code{[tmle3]{tmle3_Task}} specifying the data and NPSEM
#'  for use in constructing components required for TML estimation.
#' @param likelihood A trained \code{[tmle3]{Likelihood}}, constructed via the
#'  \code{\link{stochastic_mediation_likelihood}} helper.
#'
#' @importFrom tmle3 define_lf LF_emp LF_fit Likelihood
#'
#' @keywords internal
stochastic_mediation_likelihood <- function(tmle_task, learner_list) {
  # covariates
  W_factor <- tmle3::define_lf(tmle3::LF_emp, "W")

  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continuous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  # treatment
  A_factor <- tmle3::define_lf(tmle3::LF_fit, "A",
    learner = learner_list[["A"]],
    bound = A_bound
  )

  # outcome
  Y_factor <- tmle3::define_lf(tmle3::LF_fit, "Y",
    learner = learner_list[["Y"]],
    type = "mean"
  )

  # construct and train likelihood
  factor_list <- list(W_factor, A_factor, Y_factor)
  likelihood_def <- tmle3::Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}

###############################################################################

#' Make Task for Derived Likelihood Factor e(A,W)
#'
#' @param tmle_task A \code{[tmle3]{tmle3_Task}} specifying the data and NPSEM
#'  for use in constructing components required for TML estimation.
#' @param likelihood A trained \code{[tmle3]{Likelihood}}, constructed via the
#'  \code{\link{stochastic_mediation_likelihood}} helper.
#'
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
make_e_task <- function(tmle_task, likelihood) {
  e_task <- sl3::sl3_Task$new(
    data = tmle_task$internal_data,
    outcome = tmle_task$npsem[["A"]]$variables,
    covariates = c(
      tmle_task$npsem[["Z"]]$variables,
      tmle_task$npsem[["W"]]$variables
    ),
    folds = tmle_task$folds
  )
  return(e_task)
}

###############################################################################

#' Make Task for Derived Likelihood Factor phi(W)
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
#' @keywords internal
make_phi_task <- function(tmle_task, likelihood) {
  mediator_task_fun_factory(param_type = "medshift")(
    tmle_task, likelihood
  )
}
