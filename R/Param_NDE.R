#' Parameter for the natural direct effect
#'
#' Parameter definition class. See
#' <https://www.ncbi.nlm.nih.gov/pubmed/22499725>
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom tmle3 Param_base
#' @importFrom sl3 sl3_Task
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_NDE, observed_likelihood, ...,
#'                      outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link[tmle3]{Likelihood}}
#'           corresponding to the observed likelihood.}
#'     \item{\code{...}}{Not currently used.}
#'     \item{\code{outcome_node}}{A \code{character}, giving the name of the
#'           node that should be treated as the outcome.}
#'   }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{The counterfactual likelihood for
#'           the treatment.}
#'     \item{\code{cf_likelihood_control}}{The counterfactual likelihood for
#'           the control.}
#'     \item{\code{treatment_task}}{\code{\link[tmle3]{tmle3_Task}} created by
#'           setting the intervention to the treatment condition: do(A = 1).}
#'     \item{\code{control_task}}{\code{\link[tmle3]{tmle3_Task}} created by
#'           setting the intervention to the control condition: do(A = 0).}
#' }
#'
#' @export
Param_NDE <- R6::R6Class(
  classname = "Param_NDE",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::Param_base,
  public = list(
    initialize = function(observed_likelihood,
                          ...,
                          outcome_node = "Y") {
      # copied from standard parameter definitions
      super$initialize(observed_likelihood, list(...),
        outcome_node = outcome_node
      )
      tmle_task <- observed_likelihood$training_task
      # counterfactual tasks
      treatment_task <-
        tmle_task$generate_counterfactual_task(
          uuid = uuid::UUIDgenerate(),
          new_data = data.table(A = 1)
        )
      control_task <-
        tmle_task$generate_counterfactual_task(
          uuid = uuid::UUIDgenerate(),
          new_data = data.table(A = 0)
        )
      # counterfactual likelihoods
      treatment_lf <- define_lf(LF_static, "A", value = 1)
      control_lf <- define_lf(LF_static, "A", value = 0)
      cf_likelihood_treatment <- CF_Likelihood$new(
        observed_likelihood, treatment_lf
      )
      cf_likelihood_control <- CF_Likelihood$new(
        observed_likelihood, control_lf
      )

      # store components
      private$.treatment_task <- treatment_task
      private$.control_task <- control_task
      private$.cf_likelihood_treatment <- cf_likelihood_treatment
      private$.cf_likelihood_control <- cf_likelihood_control
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      likelihood <- self$observed_likelihood
      treatment_task <- self$treatment_task
      control_task <- self$control_task
      cf_likelihood_treatment <- self$cf_likelihood_treatment
      cf_likelihood_control <- self$cf_likelihood_control

      # compute/extract g(1|W) and g(0|W)
      g1_est <- likelihood$get_likelihood(treatment_task, "A", fold_number)
      g0_est <- likelihood$get_likelihood(control_task, "A", fold_number)

      # treatment/control indicators
      cf_pA_treatment <- cf_likelihood_treatment$get_likelihoods(
        tmle_task, "A", fold_number
      )
      cf_pA_control <- cf_likelihood_control$get_likelihoods(
        tmle_task, "A", fold_number
      )

      # compute e(1|W,Z) & e(0|W,Z)
      e_est <- likelihood$get_likelihood(tmle_task, "E", fold_number)
      e1_est <- cf_pA_treatment * e_est + cf_pA_control * (1 - e_est)
      e0_est <- 1 - e1_est

      # clever covariates
      HY <- cf_pA_treatment * (e0_est / e1_est) * g0_est / g1_est^2 -
        cf_pA_control / g0_est
      HZ <- 1 / g0_est

      # output clever covariates
      return(list(Y = HY, psi_Z = HZ))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # get observed likelihood
      likelihood <- self$observed_likelihood
      treatment_task <- self$treatment_task
      control_task <- self$control_task
      cf_likelihood_control <- self$cf_likelihood_control

      # extract various likelihood components
      y <- tmle_task$get_tmle_node(self$outcome_node)
      m_est <- likelihood$get_likelihood(tmle_task, "Y", fold_number)

      # Y estimates under treatment/control
      m1_est <- likelihood$get_likelihood(treatment_task, "Y", fold_number)
      m0_est <- likelihood$get_likelihood(control_task, "Y", fold_number)

      # predict psi_Z for full dataset
      psi_Z_data <- data.table::as.data.table(list(
        tmle_task$get_tmle_node("W")
      ))
      full_psi_Z_task <- sl3::sl3_Task$new(
        psi_Z_data,
        covariates = tmle_task$npsem[["W"]]$variables,
        outcome_type = "continuous"
      )
      psi_Z_est <- likelihood$factor_list$psi_Z$learner$predict_fold(
        full_psi_Z_task, fold_number
      )

      # clever_covariates happens here but this is repeated computation
      HY <- self$clever_covariates(
        tmle_task,
        fold_number
      )[[self$outcome_node]]
      HZ <- self$clever_covariates(
        tmle_task,
        fold_number
      )[["psi_Z"]]

      cf_pA_control <- cf_likelihood_control$get_likelihoods(
        tmle_task, "A", fold_number
      )

      # compute individual scores for DY, DA, DW
      D_Y <- HY * (y - m_est)
      D_Z <- cf_pA_control * HZ * (m1_est - m0_est - psi_Z_est)
      D_W <- psi_Z_est

      # parameter and influence function
      theta <- mean(psi_Z_est)
      eif <- D_Y + D_Z + D_W - theta

      # output
      result <- list(psi = theta, IC = eif)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf(
        "NDE[%s_{A=1} - %s_{A=0}]", self$outcome_node, self$outcome_node
      )
      return(param_form)
    },
    treatment_task = function() {
      return(private$.treatment_task)
    },
    control_task = function() {
      return(private$.control_task)
    },
    cf_likelihood_treatment = function() {
      return(private$.cf_likelihood_treatment)
    },
    cf_likelihood_control = function() {
      return(private$.cf_likelihood_control)
    },
    update_nodes = function() {
      return(c(self$outcome_node))
    }
  ),
  private = list(
    .type = "NDE",
    .treatment_task = NULL,
    .control_task = NULL,
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL
  )
)
