# use custom package library on Savio cluster
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/david_mccoy/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# packages
library(here)
library(foreach)
library(future)
library(doFuture)
library(doRNG)
library(data.table)
library(tidyverse)
library(hal9001)
library(origami)
library(sl3)
library(tmle3)
devtools::load_all('/Users/davidmccoy/Documents/PhD/BBD_capstone/tmle3mediate')

# load scripts, parallelization, PRNG
source(here("mediation/tmle3mediate/sandbox/01_setup_data.R"))
source(here("mediation/tmle3mediate/sandbox/02_fit_estimators.R"))
registerDoFuture()
plan(multiprocess)

# simulation parameters
set.seed(7259)
n_sim <- 1                                   # number of simulations
n_obs <- cumsum(rep(sqrt(100), 5))^2          # sample sizes at root-n scale
n_obs <- n_obs[c(1)]
# perform simulation across sample sizes
sim_results <- lapply(n_obs, function(sample_size) {
  # get results in parallel
  results <- foreach(this_iter = seq_len(n_sim),
                     .options.multicore = list(preschedule = FALSE),
                     .errorhandling = "remove") %do% {
    data_sim <- sim_data(n_obs = sample_size)
    est_out <- fit_estimators(data = data_sim)
    #browser()
    return(est_out)
  }

  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  return(results_out)
})

# save results to file
names(sim_results) <- paste("n", n_obs, sep = "_")
timestamp <- str_replace_all(Sys.time(), " ", "_")
saveRDS(object = sim_results,
        file = here("data", paste0("tmle3mediate_", timestamp, ".rds")))



