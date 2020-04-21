library(here)
library(tidyverse)
source(here("/tests/testthat/helper-dgp.R"))

# generate data based on DGP components
sim_data <- function(n_obs) {
    ## set up DGP
    dgp <- make_dgp()

    ## generate  baseline covariates
    W_1 <- rbinom(n_obs, 1, prob = 0.50)
    W_2 <- rbinom(n_obs, 1, prob = 0.65)
    W_3 <- rbinom(n_obs, 1, prob = 0.35)
    W <- list(W_1, W_2, W_3) %>%
      bind_cols() %>%
      set_names(paste0("W", seq_len(3)))

    ## generate treatment conditional on W
    A <- rbinom(n_obs, 1, prob = dgp$g_mech(w = W))

    # generate mediators conditional on {A, W}
    z_probs <- dgp$z_mech(a = A, w = W)
    Z <- lapply(z_probs, function(z_prob) {
      rbinom(n_obs, 1, z_prob)
    }) %>%
    bind_cols() %>%
    set_names(paste0("Z", seq_len(3)))

    ## generate outcome conditional on {Z, A, W}
    Y <- dgp$m_mech(z = Z, a = A, w = W) + rnorm(length(A), mean = 0, sd = 0.5)

    ## construct data for output
    dat <- as_tibble(cbind(W, A, Z, Y))
    return(dat)
}
