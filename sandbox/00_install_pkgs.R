# set user-specific package library
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# set CRAN mirror
options(repos = structure(c(CRAN = "https://cran.rstudio.com/")))

# lie to pkgbuild, as per Jeremy
pkgbuild:::cache_set("has_compiler", TRUE)

# from CRAN
install.packages(c("here", "tidyverse", "remotes", "future", "future.apply",
                   "doFuture", "foreach", "doRNG", "data.table", "devtools",
                   "Rsolnp", "nnls", "glmnet", "Rcpp", "origami", "hal9001"),
                 lib = "/global/scratch/nhejazi/R")

# use remotes to install from GitHub
remotes::install_github(c("tlverse/sl3@master",
                          "tlverse/tmle3mediate@master"),
                        lib = "/global/scratch/nhejazi/R")

# update all packages
update.packages(ask = FALSE, lib.loc = "/global/scratch/nhejazi/R")
