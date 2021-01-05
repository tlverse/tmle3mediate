
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`tmle3mediate`

[![Travis-CI Build
Status](https://travis-ci.com/tlverse/tmle3mediate.svg?branch=master)](https://travis-ci.com/tlverse/tmle3mediate)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/tlverse/tmle3mediate?branch=master&svg=true)](https://ci.appveyor.com/project/tlverse/tmle3mediate)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tlverse/tmle3mediate/master.svg)](https://codecov.io/github/tlverse/tmle3mediate?branch=master)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

> Targeted Learning for Causal Mediation Analysis

**Authors:** [Nima Hejaz](https://nimahejazi.org), James Duncan, David
McCoy, and [Mark van der Laan](https://vanderlaan-lab.org)

-----

## What’s `tmle3mediate`?

`tmle3mediate` is an adapter/extension R package in the `tlverse`
ecosystem that provides support for *causal mediation analysis*, for a
range of target parameters applicable in settings with mediating
variables. Causal effects for which estimation machinery are provided
include the controlled direct effect (Petersen, Sinisi, and van der Laan
2006; Didelez, Dawid, and Geneletti 2012; VanderWeele 2015), the natural
(in)direct effects (Robins and Greenland 1992; Zheng and van der Laan
2012), and the population intervention (in)direct effects (Dı́az and
Hejazi 2020). By building on the core `tlverse` grammar exposed by the
`tmle3` R package, `tmle3mediate` accommodates targeted maximum
likelihood (or targeted minimum loss-based) estimation of these causal
effect parameters through a unified interface. For a general discussion
of the framework of targeted minimum loss-based estimation and its
relationship to statistical causal inference, the motivated reader may
consider consulting van der Laan and Rose (2011) and van der Laan and
Rose (2018). A practical and accessible introduction using the `tlverse`
software ecosystem is provided in van der Laan et al. (2021) (see
<https://tlverse.org/tlverse-handbook>).

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

``` 
    @software{hejazi2021tmle3mediate-rpkg,
      author = {Hejazi, Nima S and Duncan, James and McCoy, David and
        {van der Laan}, Mark J},
      title = {{tmle3mediate}: Targeted Learning for Causal Mediation
        Analysis},
      year  = {2021},
      doi = {},
      url = {https://github.com/tlverse/tmle3mediate},
      note = {R package version 0.0.2}
    }
```

-----

## Related

  - [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
    providing tools to estimate the causal effect of stochastic
    treatment regimes in the mediation setting, including classical
    (IPW) and doubly robust one-step estimators. This is an
    implementation of the methodology explored by Dı́az and Hejazi
    (2020).

  - [R/`medoutcon`](https://github.com/nhejazi/medoutcon) - An R package
    providing doubly robust estimators (one-step, TMLE) of the
    interventional (in)direct effects, which are defined by joint
    static-stochastic interventions applied to the exposure and
    mediators, respectively. These effect definitions are similar to but
    more general than the natural (in)direct effects. This is an
    implementation of the methodology explored by Dı́az et al. (2020).

-----

## Funding

The development of this software was supported in part through [UC
Berkeley’s Biomedical Big Data training
program](http://bbd.berkeley.edu/), made possible by grant [T32
LM012417](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=)
from the National Institutes of Health.

-----

## License

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

-----

## References

<div id="refs" class="references">

<div id="ref-didelez2012direct">

Didelez, Vanessa, Philip Dawid, and Sara Geneletti. 2012. “Direct and
Indirect Effects of Sequential Treatments.” *arXiv Preprint
arXiv:1206.6840*.

</div>

<div id="ref-diaz2020causal">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2020nonparametric">

Dı́az, Iván, Nima S Hejazi, Kara E Rudolph, and Mark J van der Laan.
2020. “Non-Parametric Efficient Causal Mediation with Intermediate
Confounders.” *Biometrika*. <https://doi.org/10.1093/biomet/asaa085>.

</div>

<div id="ref-petersen2006estimation">

Petersen, Maya L, Sandra E Sinisi, and Mark J van der Laan. 2006.
“Estimation of Direct Causal Effects.” *Epidemiology*, 276–84.

</div>

<div id="ref-robins1992identifiability">

Robins, James M, and Sander Greenland. 1992. “Identifiability and
Exchangeability for Direct and Indirect Effects.” *Epidemiology*,
143–55.

</div>

<div id="ref-vdl2021targeted">

van der Laan, Mark J, Jeremy R Coyle, Nima S Hejazi, Ivana Malenica,
Rachael V Phillips, and Alan E Hubbard. 2021. *Targeted Learning in `R`:
Causal Data Science with the `tlverse` Software Ecosystem*. CRC Press.
<https://tlverse.org/tlverse-handbook>.

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-vdl2018targeted">

———. 2018. *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Science & Business Media.

</div>

<div id="ref-vanderweele2015explanation">

VanderWeele, Tyler. 2015. *Explanation in Causal Inference: Methods for
Mediation and Interaction*. Oxford University Press.

</div>

<div id="ref-zheng2012targeted">

Zheng, Wenjing, and Mark J van der Laan. 2012. “Targeted Maximum
Likelihood Estimation of Natural Direct Effects.” *International Journal
of Biostatistics* 8 (1). <https://doi.org/10.2202/1557-4679.1361>.

</div>

</div>
