
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`tmle3mediate`

[![Travis-CI Build
Status](https://travis-ci.org/tlverse/tmle3mediate.svg?branch=master)](https://travis-ci.org/tlverse/tmle3mediate)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/tlverse/tmle3mediate?branch=master&svg=true)](https://ci.appveyor.com/project/tlverse/tmle3mediate)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tlverse/tmle3mediate/master.svg)](https://codecov.io/github/tlverse/tmle3mediate?branch=master)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Targeted Learning for Causal Mediation Analysis

**Authors:** \[TO FILL IN\]

-----

## What’s `tmle3mediate`?

The `tmle3mediate` R package is designed to \[TO FILL IN\]

-----

## Installation

Install the most recent *stable release* from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("tlverse/tmle3mediate")
```

-----

## Example

To illustrate how `tmle3mediate` may be used to estimate the effect of
applying a stochastic intervention to the treatment (`A`) while keeping
the mediator(s) (`Z`) fixed, consider the following example:

``` r
library(data.table)
library(origami)
library(sl3)
library(tmle3)
library(tmle3mediate)

# produces a simple data set based on ca causal model with mediation
make_mediation_data <- function(n_obs = 1000) {
  # baseline covariate -- simple, binary
  W <- rbinom(n_obs, 1, prob = 0.50)

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = W / 4 + 0.1))

  # single mediator to affect the outcome
  z1_prob <- 1 - plogis((A^2 + W) / (A + W^3 + 0.5))
  Z <- rbinom(n_obs, 1, prob = z1_prob)

  # create outcome as a linear function of A, W + white noise
  Y <- Z + A - 0.1 * W + rnorm(n_obs, mean = 0, sd = 0.25)

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c("Y", "Z", "A", "W"))
  return(data)
}

# set seed and simulate example data
set.seed(75681)
example_data <- make_mediation_data(100)
node_list <- list(W = "W", A = "A", Z = "Z", Y = "Y")

# consider an incremental propensity score intervention that triples (i.e.,
# delta = 3) the individual-specific odds of receiving treatment
delta_ipsi <- 3

# make learners for nuisance parameters
g_learners <- e_learners <- m_learners <- phi_learners <-
  Lrnr_cv$new(Lrnr_glm$new(), full_fit = TRUE)
learner_list <- list(Y = m_learners, A = g_learners)

# compute one-step estimate for an incremental propensity score intervention
tmle_spec <- tmle_medshift(delta = delta_ipsi,
                           e_learners = e_learners,
                           phi_learners = phi_learners,
                           max_iter = 5)
tmle_out <- tmle3(tmle_spec, example_data, node_list, learner_list)
tmle_out
#> A tmle3_Fit that took 5 step(s)
#>    type         param  init_est  tmle_est        se     lower   upper
#> 1: PIDE E[Y_{A=NULL}] 0.7938906 0.7927064 0.2039544 0.3929632 1.19245
#>    psi_transformed lower_transformed upper_transformed
#> 1:       0.7927064         0.3929632           1.19245
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tlverse/tmle3mediate/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/tlverse/tmle3mediate/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `tmle3mediate` R package, please cite the following:

-----

## License

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

-----

## References
