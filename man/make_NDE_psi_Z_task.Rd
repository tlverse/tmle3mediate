% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmle3_Spec_NDE.R
\name{make_NDE_psi_Z_task}
\alias{make_NDE_psi_Z_task}
\title{Make Task for Derived Likelihood Factor psi_Z(W) for Natural Direct Effect}
\usage{
make_NDE_psi_Z_task(tmle_task, likelihood)
}
\arguments{
\item{tmle_task}{A \code{[tmle3]{tmle3_Task}} specifying the data and NPSEM
for use in constructing components required for TML estimation.}

\item{likelihood}{A trained \code{[tmle3]{Likelihood}}, constructed via the
\code{\link{stochastic_mediation_likelihood}} helper.}
}
\description{
Make Task for Derived Likelihood Factor psi_Z(W) for Natural Direct Effect
}
\keyword{internal}
