% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmle3_Spec_medshift.R
\name{make_phi_task}
\alias{make_phi_task}
\title{Make Task for Derived Likelihood Factor phi(W)}
\usage{
make_phi_task(tmle_task, likelihood)
}
\arguments{
\item{tmle_task}{A \code{[tmle3]{tmle3_Task}} specifying the data and NPSEM
for use in constructing components required for TML estimation.}

\item{likelihood}{A trained \code{[tmle3]{Likelihood}}, constructed via the
\code{\link{stochastic_mediation_likelihood}} helper.}
}
\description{
Make Task for Derived Likelihood Factor phi(W)
}
\keyword{internal}
