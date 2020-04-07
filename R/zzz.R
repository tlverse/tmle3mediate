.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "tmle3mediate v", utils::packageDescription("tmle3mediate")$Version,
    ": Targeted Learning for Causal Mediation Analysis"
  ))
}
