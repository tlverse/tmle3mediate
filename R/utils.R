################################################################################

#' Factory function that returns an \code{\link[sl3]{sl3_Task}} generator
#' function used in fitting the appropriate mediator derived conditional
#' likelihood factor for NDE/NIE/medshift.
#'
#' @param param_type A \code{character}, giving the target parameter type for
#'   which to construct the psi_Z \code{[tmle3]{tmle3_Task}}.
#' @param subset_value A \code{numeric}, which determines the treatment value
#'   used to subset the \code{[tmle3]{tmle3_Task}}. \code{NA} by default for no
#'   subsetting.
#'
#' @return A \code{function} that takes a \code{[tmle3]{tmle3_Task}} and a
#'   \code{[tmle3]{Likelihood}} and returns a new \code{[sl3]{sl3_Task}}.
#'
#' @importFrom data.table as.data.table data.table setnames
#' @importFrom uuid UUIDgenerate
#' @importFrom sl3 sl3_Task
#'
#' @name mediator_task_fun_factory
#'
#' @keywords internal
mediator_task_fun_factory <- function(param_type = c("NDE", "NIE", "medshift"),
                                      subset_value = NA) {
  # determine the target parameter type
  param_type <- match.arg(param_type)

  # create the task generator
  make_mediator_task <- function(tmle_task, likelihood) {
    # create treatment task for intervention conditions
    # TODO: remove node name hard-coding
    treatment_task <-
      tmle_task$generate_counterfactual_task(
        uuid = uuid::UUIDgenerate(),
        new_data = data.table::data.table(A = 1)
      )
    # create counterfactual outcomes and construct pseudo-outcome
    m_pseudo <- likelihood$get_likelihood(treatment_task, "Y")

    if (param_type %in% c("NDE", "medshift")) {
      # create control tasks for intervention conditions
      control_task <-
        tmle_task$generate_counterfactual_task(
          uuid = uuid::UUIDgenerate(),
          new_data = data.table::data.table(A = 0)
        )
      m0 <- likelihood$get_likelihood(control_task, "Y")
      # subtract control counterfactual likelihood for NDE mediator factor
      m_pseudo <- m_pseudo - m0
    }

    # create regression task for pseudo-outcome and baseline covariates
    mediator_data <- data.table::as.data.table(list(
      m_pseudo = m_pseudo,
      tmle_task$get_tmle_node("W")
    ))
    data.table::setnames(mediator_data,
                         c("m_pseudo", tmle_task$npsem[["W"]]$variables))
    mediator_task <- sl3::sl3_Task$new(
      data = mediator_data,
      outcome = "m_pseudo",
      covariates = tmle_task$npsem[["W"]]$variables,
      outcome_type = "continuous",
      folds = tmle_task$folds
    )

    if (is.na(subset_value)) {
      # no subsetting
      return(mediator_task)
    } else {
      # subset data: observations with treatment equal to subset_value
      row_index <- tmle_task$get_tmle_node("A") == subset_value
      return(mediator_task$subset_task(row_index))
    }
  }
  return(make_mediator_task)
}
