st_bayesc_scheduled_st_bayesc_scheduled_base_helper <- if (file.exists(file.path("fixtures","st-bayesc-scheduled-reference-base.R")))
  file.path("fixtures","st-bayesc-scheduled-reference-base.R") else
  file.path("tests","testthat","fixtures","st-bayesc-scheduled-reference-base.R")
source(st_bayesc_scheduled_st_bayesc_scheduled_base_helper)

st_bayesc_scheduled_starting_commit <- "fb9fd03"
st_bayesc_scheduled_configs <- st_bayesc_scheduled_base_configs
st_bayesc_scheduled_run <- st_bayesc_scheduled_base_run
st_bayesc_scheduled_normalize <- st_bayesc_scheduled_base_normalize
st_bayesc_scheduled_metadata <- function(name,cfg) list(
  starting_commit=st_bayesc_scheduled_starting_commit, model="scheduled CSR BayesC",
  configuration=name, marker_count=6L, trait_count=1L, seed=1001L,
  chain_seeds=cfg$seeds,nchains=cfg$nchains,ncores=cfg$ncores,
  scheduler_controls=cfg,nit=8L,nburn=2L,nthin=1L,keep_chains=FALSE,
  schema="stblr_raw_v1",schema_version=1L,
  rng_ownership_version="scheduled_chain_rng_v1",
  engine_owner="chain",distribution_owner="chain",
  reference_mode="fresh R process",R_version=R.version.string,
  compiler="Rtools44 GCC 13.2 C++17")
