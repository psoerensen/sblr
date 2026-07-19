phase17c_mt_config <- phase17b_mt_config
phase17c_mt_capture <- phase17b_mt_capture
phase17c_mt_native_raw <- phase17b_mt_native_raw

phase17c_mt_retained_counts <- function(config) {
  iteration <- 0:(config$nit + config$nburn - 1L)
  post_burn <- iteration >= config$nburn
  list(
    marker_retained_count = sum(post_burn &
      ((iteration - config$nburn) %% config$nthin == 0L)),
    covb_retained_count = if (isTRUE(config$updateB)) config$nit else 0L,
    covg_retained_count = config$nit,
    cove_retained_count = if (isTRUE(config$updateE)) config$nit else 0L,
    pi_retained_count = if (isTRUE(config$updatePi)) config$nit else 0L
  )
}

phase17c_mt_metadata <- function(id = 1L) {
  config <- phase17c_mt_config(id)
  c(list(
    starting_commit = "19a64e0",
    correction_phase = "Phase 17C",
    reference_mode = "structure_exact_numeric_tolerance",
    numeric_tolerance = 1e-12,
    structure_exact = TRUE,
    schema = "legacy positional native / named public fit",
    seed_policy = "R-generated seed; one fit-local std::mt19937",
    samples = config$n,
    markers = lengths(config$Xy)[1],
    traits = length(config$Xy),
    iterations = config$nit,
    burnin = config$nburn,
    thinning = config$nthin,
    updateB = config$updateB,
    updateE = config$updateE,
    updatePi = config$updatePi
  ), phase17c_mt_retained_counts(config))
}
