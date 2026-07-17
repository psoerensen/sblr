phase10b_phase10a_helper <- if (file.exists(file.path("fixtures","blr-phase10a-scheduled-reference.R")))
  file.path("fixtures","blr-phase10a-scheduled-reference.R") else
  file.path("tests","testthat","fixtures","blr-phase10a-scheduled-reference.R")
source(phase10b_phase10a_helper)

phase10b_starting_commit <- "fb9fd03"
phase10b_configs <- phase10a_configs
phase10b_run <- phase10a_run
phase10b_normalize <- phase10a_normalize
phase10b_metadata <- function(name,cfg) list(
  starting_commit=phase10b_starting_commit, model="scheduled CSR BayesC",
  configuration=name, marker_count=6L, trait_count=1L, seed=1001L,
  chain_seeds=cfg$seeds,nchains=cfg$nchains,ncores=cfg$ncores,
  scheduler_controls=cfg,nit=8L,nburn=2L,nthin=1L,keep_chains=FALSE,
  schema="stblr_raw_v1",schema_version=1L,
  rng_ownership_version="scheduled_chain_rng_v1",
  engine_owner="chain",distribution_owner="chain",
  reference_mode="fresh R process",R_version=R.version.string,
  compiler="Rtools44 GCC 13.2 C++17")
