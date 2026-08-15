phase3_public_csr_prefix <- function(m) {
  prefix <- tempfile("phase3_csr_ld_")
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), rep(0, m + 1L))
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA",
    paste0("n_variants=", m), "nnz=0", "triangle=upper",
    "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"), paste0(prefix, ".meta.txt"))
  prefix
}

phase3_public_stats <- function(traits) {
  markers <- paste0("m", 1:3)
  wy <- lapply(seq_along(traits), function(index) {
    stats::setNames(c(2, -1, .5) * index, markers)
  })
  ww <- lapply(traits, function(.) stats::setNames(rep(50, 3), markers))
  names(wy) <- names(ww) <- traits
  list(wy = wy, ww = ww, yy = stats::setNames(rep(50, length(traits)), traits),
       n = 50L, m = 3L, marker_names = markers, trait_names = traits)
}

test_that("Phase 3 seed contract matches every frozen R and native vector", {
  cases <- data.frame(
    user_seed = c(0, 0, 0, 0, 0, 17, 17, 17),
    mode = c("single_trait", "single_trait", "independent_traits",
             "independent_traits", "joint_multitrait", "single_trait",
             "independent_traits", "joint_multitrait"),
    trait = c(NA, NA, "traitA", "traitB", NA, NA, "traitB", NA),
    chain = c(0, 1, 0, 1, 0, 0, 1, 1),
    expected = c(830191578, 160141543, 226943096, 286956759,
                 3100589946, 3397578794, 1132619387, 3700933392),
    stringsAsFactors = FALSE)
  actual_r <- mapply(
    sblr:::.blr_seed_v1, cases$user_seed, cases$mode, cases$trait,
    cases$chain)
  actual_native <- mapply(function(seed, mode, trait, chain) {
    sblr:::blr_phase3_seed_v1_internal(
      seed, mode, if (is.na(trait)) NULL else trait, chain)
  }, cases$user_seed, cases$mode, cases$trait, cases$chain)
  expect_identical(unname(actual_r), cases$expected)
  expect_identical(unname(actual_native), cases$expected)
  expect_true(any(actual_r > .Machine$integer.max))
  boundary <- sblr:::.blr_seed_v1(
    4294967295, "single_trait", NA_character_, 0L)
  expect_true(is.double(boundary) && boundary >= 0 && boundary <= 4294967295)
})

test_that("logical tasks are canonical and stable under trait reordering", {
  plan <- sblr:::.blr_logical_task_plan(
    "independent_traits", c("traitB", "traitA"), 2L)
  expect_identical(plan$task_id, c(
    "trait:traitB|chain:0", "trait:traitB|chain:1",
    "trait:traitA|chain:0", "trait:traitA|chain:1"))
  first <- sblr:::.blr_task_seeds_v1(
    0, "independent_traits", c("traitA", "traitB"), 2L)
  reordered <- sblr:::.blr_task_seeds_v1(
    0, "independent_traits", c("traitB", "traitA"), 2L)
  expect_identical(first["traitA", ], reordered["traitA", ])
  expect_identical(first["traitB", ], reordered["traitB", ])
  expect_identical(
    sblr:::.blr_logical_task_plan("joint_multitrait", c("T1", "T2"), 2L)$task_id,
    c("chain:0", "chain:1"))
})

test_that("explicit Phase 3 task seeds preserve uint32 values and identities", {
  chain <- sblr:::.blr_chain_controls(
    nit = 4L, nburn = 2L, nthin = 2L, seed = 0,
    nchains = 2L, ncores = 2L,
    chain_seeds = c(0, 4294967295))
  spec <- sblr:::resolve_blr_spec_from_wrapper(
    "sbayesc", "csr", c("traitA", "traitB"), c("m1", "m2"),
    chain, c(20, 21))
  expect_true(sblr:::validate_blr_resolved_spec(spec))
  contract <- sblr:::.blr_native_execution_contract(spec)
  expect_identical(contract$task_ids, c(
    "trait:traitA|chain:0", "trait:traitA|chain:1",
    "trait:traitB|chain:0", "trait:traitB|chain:1"))
  expect_true(all(contract$task_seeds >= 0 & contract$task_seeds <= 4294967295))
  bad <- spec
  bad$mcmc$task_seeds[1, 1] <- 4294967296
  expect_error(sblr:::validate_blr_resolved_spec(bad), "uint32")
})

test_that("retention and convergence plans have distinct exact indices", {
  examples <- list(
    c(1L, 1L), c(5L, 2L), c(4L, 4L), c(3L, 4L))
  expected <- list(1L, c(2L, 4L), 4L, integer())
  for (index in seq_along(examples)) {
    value <- sblr:::.blr_retention_plan(
      7L, examples[[index]][1], examples[[index]][2],
      retained_requested = FALSE)
    expect_identical(value$post_burn, expected[[index]])
    expect_identical(value$absolute_transition, 7L + expected[[index]])
  }
  expect_error(sblr:::.blr_retention_plan(
    0L, 3L, 4L, retained_requested = TRUE), "no draws")
  convergence <- sblr:::.blr_convergence_iteration_plan(7L, 5L)
  expect_identical(convergence$post_burn, 1:5)
  expect_identical(convergence$absolute_transition, 8:12)
  expect_identical(convergence$rng_draws, 0L)
})

test_that("legacy retention remains an explicit qualification boundary", {
  expect_identical(
    sblr:::.blr_retention_plan(0L, 5L, 2L, 0L)$post_burn,
    c(1L, 3L, 5L))
  expect_identical(
    sblr:::.blr_retention_plan(0L, 5L, 2L, 1L)$post_burn,
    c(2L, 4L))
})

test_that("the unified scheduler is neutral for identical seeds and retention", {
  fixture <- new.env(parent = globalenv())
  sys.source(test_path("fixtures", "st-bayesc-csr-reference.R"), envir = fixture)
  config <- fixture$st_bayesc_csr_reference_configurations$
    one_trait_one_chain_one_core
  prefix <- fixture$st_bayesc_csr_reference_write_csr()
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  inputs <- fixture$st_bayesc_csr_reference_inputs(config$traits)
  legacy_seed <- as.numeric(sblr:::blr_scalar_seeds_cpp(
    31L, 1L, 1L, integer())[1L])
  if (legacy_seed < 0) legacy_seed <- legacy_seed + 4294967296
  contract <- list(
    seed_contract_version = 1L,
    retention_contract_version = 1L,
    scheduler_version = 1L,
    task_seeds = legacy_seed,
    task_ids = "trait:trait1|chain:0",
    retained_transition_indices = 1:8)
  legacy <- fixture$st_bayesc_csr_reference_native(config, prefix, inputs)
  unified <- fixture$st_bayesc_csr_reference_native(
    config, prefix, inputs, execution_contract = contract)
  legacy$diagnostics$seconds_mean <- unified$diagnostics$seconds_mean <- 0
  legacy$diagnostics$seconds_max <- unified$diagnostics$seconds_max <- 0
  legacy$diagnostics["workers"] <- unified$diagnostics["workers"] <- list(NULL)
  expect_identical(legacy, unified)
})

test_that("eligible CSR fits execute and report Phase 3 contracts", {
  prefix <- phase3_public_csr_prefix(3L)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  run <- function(cores) stblr_csr(
    phase3_public_stats(c("traitA", "traitB")), ld_prefix = prefix,
    method = "sbayesc", nit = 5L, nburn = 1L, nthin = 2L,
    ncores = cores, seed = 0, updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, convergence = "none")
  serial <- run(1L)
  parallel <- run(2L)
  expect_identical(serial$bm, parallel$bm)
  expect_identical(serial$dm, parallel$dm)
  expect_identical(serial$vbs, parallel$vbs)
  expect_identical(serial$input$task_seeds_resolved,
                   parallel$input$task_seeds_resolved)
  expect_identical(serial$input$retained_transition_indices, c(2L, 4L))
  expect_identical(serial$input$convergence_iteration_indices, 1:5)
  expect_identical(serial$input$seed_contract_version, 1L)
  expect_identical(serial$input$retention_contract_version, 1L)
  expect_identical(serial$input$scheduler_version, 1L)
  workers <- parallel$diagnostics$native$workers
  expect_identical(workers$requested_cores, 2L)
  expect_identical(workers$scheduler_version, 1L)
  expect_identical(workers$diagnostics_rng_draws, 0L)
  expect_identical(workers$logical_task_order,
                   c("trait:traitA|chain:0", "trait:traitB|chain:0"))
  expect_true(workers$actual_team_size >= 1L)
  expect_identical(length(workers$task_worker_ids), 2L)
  capability <- sblr:::sparseLD_thread_info(2L)
  if (isTRUE(capability$openmp) &&
      capability$actual_threads_requested_region > 1L) {
    expect_gt(workers$actual_team_size, 1L)
    expect_gt(length(unique(workers$task_worker_ids)), 1L)
  }
})

test_that("formal convergence capture is observational and unthinned", {
  prefix <- phase3_public_csr_prefix(3L)
  on.exit(unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt"))), add = TRUE)
  run <- function(keep) stblr_csr(
    phase3_public_stats("traitA"), ld_prefix = prefix,
    method = "sbayesc", nit = 5L, nburn = 1L, nthin = 2L,
    ncores = 1L, seed = 17, updateB = FALSE, updateE = FALSE,
    updatePi = FALSE, convergence = "core",
    convergence_control = list(keep_traces = keep, warn = FALSE))
  discarded <- run(FALSE)
  retained <- run(TRUE)
  expect_identical(discarded$bm, retained$bm)
  expect_identical(discarded$dm, retained$dm)
  expect_identical(discarded$vbs, retained$vbs)
  expect_null(discarded$convergence_traces)
  expect_identical(dim(retained$convergence_traces$values)[1L], 5L)
  expect_identical(retained$input$convergence_iteration_indices, 1:5)
  expect_identical(retained$input$retained_transition_indices, c(2L, 4L))
})
