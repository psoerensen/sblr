args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else
  "research/blr_framework_contract/test_contract_fixtures.R"
source(file.path(dirname(normalizePath(script)), "blr_contract_fixtures.R"))

expectations <- 0L
expect_true <- function(x, label) {
  expectations <<- expectations + 1L
  if (!isTRUE(x)) stop("FAILED: ", label, call. = FALSE)
}
expect_identical <- function(x, y, label) {
  expectations <<- expectations + 1L
  if (!identical(x, y)) stop("FAILED: ", label, call. = FALSE)
}
expect_error <- function(expr, pattern, label) {
  expectations <<- expectations + 1L
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  if (is.null(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop("FAILED: ", label, if (is.null(message)) " (no error)" else paste0(" (", message, ")"), call. = FALSE)
  }
}

# Retention-contract version 1.
r <- blr_retained_indices_v1(5, 7, 3)
expect_identical(r$post_burn, c(3L, 6L), "post-burn retained indices")
expect_identical(r$absolute_transition, c(8, 11), "absolute retained transitions")
expect_identical(r$retained_draws, 2L, "retained draw count")
expect_error(blr_retained_indices_v1(0, 1, 2), "keeps no draws", "zero retained draws")
expect_error(blr_retained_indices_v1(0.5, 2, 1), "finite integer", "fractional burn-in")
expect_error(blr_retained_indices_v1(0, 0, 1), "finite integer", "zero sampling iterations")

# Analysis/execution policy.
expect_true(blr_valid_execution("single_trait", "serial", "none"), "single serial")
expect_true(blr_valid_execution("single_trait", "parallel", "chains"), "single parallel chains")
expect_true(blr_valid_execution("independent_traits", "parallel", "trait_chains"), "independent trait-chains")
expect_true(blr_valid_execution("joint_multitrait", "parallel", "chains"), "joint parallel chains")
expect_true(!blr_valid_execution("joint_multitrait", "parallel", "traits"), "joint traits rejected")
expect_true(!blr_valid_execution("single_trait", "serial", "chains"), "serial requires none")

# Fixed seed-contract vectors. Values are populated and independently checked
# against the exact algorithm during Phase 0 closeout.
reference <- blr_seed_reference_vectors_v1()
derived <- mapply(.blr_seed_v1, reference$user_seed, reference$analysis_mode,
                  reference$trait_id, reference$chain_index)
expect_true(all(is.finite(derived) & derived >= 0 & derived <= 4294967295), "native seeds are uint32")
expect_identical(unname(derived), reference$native_seed, "fixed seed-contract reference vectors")
seed_a <- blr_task_seeds_v1(0, "independent_traits", c("traitA", "traitB"), 2)
seed_b <- blr_task_seeds_v1(0, "independent_traits", c("traitB", "traitA"), 2)
expect_identical(seed_a["traitA", ], seed_b["traitA", ], "trait reorder preserves stable task seeds")
expect_identical(seed_a["traitB", ], seed_b["traitB", ], "second trait reorder invariant")
expect_identical(blr_task_seeds_v1(17, "single_trait", "traitA", 2),
                 blr_task_seeds_v1(17, "single_trait", "renamed", 2),
                 "single-trait sentinel is stable")
expect_identical(dim(seed_a), c(2L, 2L), "independent task seed shape")
expect_identical(dim(blr_task_seeds_v1(0, "joint_multitrait", c("t1", "t2"), 2)), NULL,
                 "joint task seed shape")
expect_identical(blr_task_seeds_v1(0, "independent_traits", c("traitA", "traitB"), 2),
                 blr_task_seeds_v1(0, "independent_traits", c("traitA", "traitB"), 2),
                 "execution mode does not enter task-seed derivation")

# Resolved specifications and provider/resource ownership.
single <- blr_spec_fixture("single_trait", "fixed_full")
independent <- blr_spec_fixture("independent_traits", "fixed_full")
joint_fixed <- blr_spec_fixture("joint_multitrait", "fixed_full")
joint_sampled <- blr_spec_fixture("joint_multitrait", "sampled_full")
summary_spec <- blr_summary_spec_fixture()
for (spec in list(single, independent, joint_fixed, joint_sampled, summary_spec)) {
  expect_true(blr_validate_spec_fixture(spec), paste("valid spec", spec$data$analysis_mode, spec$model$residual_policy))
}
parallel_independent <- independent
parallel_independent$compute$execution_mode <- "parallel"
parallel_independent$compute$parallelization <- "trait_chains"
parallel_independent$compute$cores <- 2L
expect_true(blr_validate_spec_fixture(parallel_independent), "parallel independent-trait spec")
expect_identical(parallel_independent$mcmc$task_seeds, independent$mcmc$task_seeds,
                 "serial and parallel specs use identical logical task seeds")
expect_identical(names(single), c("schema", "data", "model", "prior", "mcmc", "compute", "output"), "exact resolved envelope names")
expect_identical(names(single$schema), c("name", "version", "compatibility_id", "seed_contract_version", "retention_contract_version", "dimension_contract_version"), "exact resolved schema names")
expect_identical(names(single$data$operator_resources$bed_shared), c("resource_id", "operator_type", "marker_ids", "alleles", "genotype_coding", "centering", "standardization", "operator_scale", "storage", "block_eigen", "approximation", "provenance"), "exact operator-resource names")
expect_identical(names(single$data$providers$p1), c("provider_id", "trait_ids", "operator_resource_id", "local_to_global", "sufficient_statistics", "sample_size", "likelihood_regime", "residual_contract", "population", "effect_scale", "overlap_group", "provenance"), "exact provider names")
expect_identical(length(independent$data$operator_resources), 1L, "shared BED resource is not duplicated")
expect_true(all(vapply(independent$data$providers, `[[`, character(1), "operator_resource_id") == "bed_shared"), "independent providers share BED")
expect_identical(length(joint_fixed$data$providers), 1L, "joint common-sample provider is not factorized")
expect_identical(joint_fixed$data$providers$p_joint$trait_ids, c("trait1", "trait2"), "joint provider owns both traits")
expect_true(!identical(summary_spec$data$providers$p1$sample_size, summary_spec$data$providers$p2$sample_size), "summary provider sample sizes differ")
expect_true(!identical(names(summary_spec$data$providers$p1$local_to_global), names(summary_spec$data$providers$p2$local_to_global)), "summary marker coverage differs")
expect_true(is.null(joint_fixed$prior$residual_covariance$degrees_of_freedom) && !is.null(joint_fixed$prior$residual_covariance$fixed_value), "fixed full Ve contract")
expect_true(isTRUE(joint_sampled$prior$residual_covariance$sampled) && is.null(joint_sampled$prior$residual_covariance$fixed_value), "sampled full Ve contract")

bad_execution <- independent; bad_execution$compute$execution_mode <- "serial"; bad_execution$compute$parallelization <- "traits"
expect_error(blr_validate_spec_fixture(bad_execution), "invalid analysis/execution", "invalid execution fails")
bad_seed <- independent; bad_seed$mcmc$task_seeds <- c(1, 2)
expect_error(blr_validate_spec_fixture(bad_seed), "invalid shape", "malformed task seed array")
bad_seed_names <- independent; dimnames(bad_seed_names$mcmc$task_seeds)[[1]] <- c("wrong1", "wrong2")
expect_error(blr_validate_spec_fixture(bad_seed_names), "invalid trait or chain dimnames", "malformed task seed dimnames")
bad_seed_value <- independent; bad_seed_value$mcmc$task_seeds[1, 1] <- 4294967296
expect_error(blr_validate_spec_fixture(bad_seed_value), "invalid shape or uint32", "out-of-range task seed")
zero_seed <- independent; zero_seed$mcmc$task_seeds[1, 1] <- 0
expect_true(blr_validate_spec_fixture(zero_seed), "explicit native seed zero is valid")
bad_reference <- independent; bad_reference$data$providers$p1$operator_resource_id <- "missing"
expect_error(blr_validate_spec_fixture(bad_reference), "unknown operator", "unknown resource reference")
duplicate_resource <- independent
duplicate_resource$data$operator_resources <- list(
  bed_shared = duplicate_resource$data$operator_resources$bed_shared,
  bed_duplicate = duplicate_resource$data$operator_resources$bed_shared)
expect_error(blr_validate_spec_fixture(duplicate_resource), "resource IDs must be unique", "duplicate resource IDs")
factorized_joint <- joint_fixed
factorized_joint$data$providers <- list(
  p1 = blr_provider("p1", "trait1", "bed_shared", factorized_joint$data$global_markers, factorized_joint$data$global_markers, "common_sample", 100),
  p2 = blr_provider("p2", "trait2", "bed_shared", factorized_joint$data$global_markers, factorized_joint$data$global_markers, "common_sample", 100)
)
names(factorized_joint$data$providers) <- c("p1", "p2")
expect_error(blr_validate_spec_fixture(factorized_joint), "non-factorized", "joint provider cannot be factorized")

# Raw schema version 2 and non-dropping axes.
raw_single <- blr_raw_fixture("single_trait")
raw_independent <- blr_raw_fixture("independent_traits")
raw_joint <- blr_raw_fixture("joint_multitrait")
for (raw in list(raw_single, raw_independent, raw_joint)) {
  expect_true(blr_validate_raw_fixture(raw), paste("valid raw", raw$input$data$analysis_mode))
}
expect_identical(names(raw_single), c("schema", "model", "input", "posterior", "draws", "final", "derived", "diagnostics", "provenance"), "exact raw envelope names")
expect_identical(dim(raw_single$draws$realised_effects), c(1L, 1L, 3L, 1L), "single raw effect axes retained")
expect_identical(attr(raw_single$draws$realised_effects, "dim_axis_names"), c("draw", "chain", "marker", "trait"), "effect axis names")
expect_identical(dim(raw_single$draws$marker_variance), c(1L, 1L, 1L), "scalar variance axes retained")
expect_true(is.null(raw_independent$draws$joint_states), "independent joint states are present NULL")
expect_true(!is.null(raw_independent$draws$independent_trait_states), "independent traitwise states present")
expect_true(!is.null(raw_joint$draws$joint_states), "joint states present")
expect_true(is.null(raw_joint$draws$independent_trait_states), "joint traitwise states present NULL")
expect_true(all(c("latent_effects", "scaled_effects", "regional_marker_covariance") %in% names(raw_single$draws)), "required-present NULL draw fields")
expect_true(!any(c("pi", "pis", "pim", "state_probabilities", "pattern_probabilities") %in% names(raw_joint$posterior)), "probability names unambiguous")

cat("PASS: blr framework contract fixtures (", expectations, " expectations)\n", sep = "")
