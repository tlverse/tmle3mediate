################################################################################

#' Factory function that returns an \code{[sl3]{sl3_Task}} generator function
#' used in fitting the psi_Z(W) derived likelihood factor for NDE/NIE.
#'
#' @param tmle_task A \code{[tmle3]{tmle3_Task}} specifying the data and NPSEM
#'   for use in constructing components required for TML estimation.
#' @param likelihood A trained \code{[tmle3]{Likelihood}}, constructed via the
#'   \code{\link{stochastic_mediation_likelihood}} helper.
#' @param subset_value A \code{numeric}, which determines the treatment value used to
#'   subset the \code{[tmle3]{tmle3_Task}}.
#' @param param_type A \code{character}, giving the target parameter type for
#'   which to construct the psi_Z \code{[tmle3]{tmle3_Task}}.
#'
#' @return A \code{function} that takes a \code{[tmle3]{tmle3_Task}} and a
#'   \code{[tmle3]{Likelihood}} and returns a new \code{[sl3]{sl3_Task}}.
#'
#' @importFrom data.table as.data.table data.table setnames
#' @importFrom uuid UUIDgenerate
#' @importFrom sl3 sl3_Task
#'
#' @name psi_Z_task_factory
#'
#' @keywords internal
psi_Z_task_factory <- function(subset_value, param_type = c("NDE", "NIE")) {
  # determine the target parameter type
  param_type <- match.arg(param_type)

  # create the task generator
  make_psi_Z_task <- function(tmle_task, likelihood) {
    # create treatment task for intervention conditions
    # TODO: remove node name hard-coding
    treatment_task <-
      tmle_task$generate_counterfactual_task(
        uuid = uuid::UUIDgenerate(),
        new_data = data.table::data.table(A = 1)
      )
    # create counterfactual outcomes and construct pseudo-outcome
    m_pseudo <- likelihood$get_likelihood(treatment_task, "Y")

    if (param_type == "NDE") {
      # create control tasks for intervention conditions
      control_task <-
        tmle_task$generate_counterfactual_task(
          uuid = uuid::UUIDgenerate(),
          new_data = data.table::data.table(A = 0)
        )
      m0 <- likelihood$get_likelihood(control_task, "Y")
      # subtract control counterfactual likelihood for NDE psi_Z factor
      m_pseudo <- m_pseudo - m0
    }

    # create regression task for pseudo-outcome and baseline covariates
    psi_Z_data <- data.table::as.data.table(list(
      m_pseudo = m_pseudo,
      tmle_task$get_tmle_node("W")
    ))
    data.table::setnames(psi_Z_data,
                         c("m_pseudo", tmle_task$npsem[["W"]]$variables))
    psi_Z_task <- sl3::sl3_Task$new(
      data = psi_Z_data,
      outcome = "m_pseudo",
      covariates = tmle_task$npsem[["W"]]$variables,
      outcome_type = "continuous",
      folds = tmle_task$folds
    )

    # subset data: observations with treatment equal to subset_value
    row_index <- tmle_task$get_tmle_node("A") == subset_value
    return(psi_Z_task$subset_task(row_index))
  }
  return(make_psi_Z_task)
}
