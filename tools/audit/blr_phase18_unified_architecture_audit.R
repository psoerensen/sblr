root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
source_files <- list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)
source_text <- paste(vapply(source_files, function(x)
  paste(readLines(x, warn = FALSE), collapse = "\n"), character(1)),
  collapse = "\n")
ns <- read("NAMESPACE")
mt_cpp <- read("src/mtblr.cpp")
unified_r <- read("R/blr-unified.R")
convergence_r <- read("R/mtblr-convergence.R")
mt_default_core <- read("src/blr_mt_default_core_impl.h")
mt_bed_core <- read("src/blr_mt_bed_core_impl.h")
scheduled_types <- read("src/blr_csr_scheduled_bayesc_types.h")
scheduled_binding <- read("src/st_cpg_omp_csr_scheduled.cpp")
unified_tests <- read("tests/testthat/test-blr-unified-convergence.R")
reduction_tests <- read("tests/testthat/test-blr-operator-reductions.R")
reduction_source <- paste(read("tests/testthat/helper-blr-unified.R"),
                          reduction_tests)
slice <- function(text, begin, end) {
  start <- regexpr(begin, text, fixed = TRUE)[1L]
  stopifnot(start > 0L)
  tail <- substring(text, start)
  finish <- regexpr(end, tail, fixed = TRUE)[1L]
  if (finish < 1L) tail else substring(tail, 1L, finish - 1L)
}
csr_dispatch <- slice(mt_cpp, "Rcpp::List mtblr_csr_chains_raw_internal(",
                      "// Named schema adapter for the public R mtblr_block_eigen()")
block_dispatch <- slice(mt_cpp,
                        "Rcpp::List mtblr_block_eigen_chains_raw_internal(",
                        "namespace {\n\nstd::vector<std::string> mt_bed_character_vector")
exports <- sub("^export\\((.*)\\)$", "\\1",
               grep("^export\\(", strsplit(ns, "\n", fixed = TRUE)[[1]],
                    value = TRUE))
canonical <- c("stblr_csr", "stblr_block_eigen", "stblr_bed",
               "stblr_csr_annot", "mtblr_csr", "mtblr_block_eigen",
               "mtblr_bed")
redundant <- c("sblr", "stblr_bed_marker", "stblr_csr_bayesr",
               "stblr_csr_prior_annot", "stblr_csr_learn_annot",
               "stblr_csr_group_annot", "stblr_csr_sbayesrc_generic")
values <- c(
  canonical_public_fitting_exports = sum(canonical %in% exports),
  model_operator_capability_matrix = file.exists(
    file.path(root, "docs/dev/blr_unified_architecture.md")),
  canonical_model_spellings = all(c("bayesc", "sbayesc", "bayesr",
                                     "sbayesr", "bayesrc", "sbayesrc") %in%
    strsplit(read("docs/dev/blr_naming_conventions.md"), "[^a-z_]+")[[1]]),
  canonical_operator_identifiers = all(c("csr", "block_eigen", "packed_bed") %in%
    strsplit(read("docs/dev/blr_naming_conventions.md"), "[^a-z_]+")[[1]]),
  redundant_wrapper_count = sum(redundant %in% exports),
  experimental_fitting_export_count = sum("stblr_bed_marker" %in% exports),
  st_formatter_count = length(gregexpr("\\.as_stblr_fit <- function", source_text)[[1]]),
  mt_formatter_count = length(gregexpr("\\.as_mtblr_fit <- function", source_text)[[1]]),
  shared_convergence_engine_count = length(gregexpr(
    "\\.blr_convergence_scalar <- function", source_text)[[1]]),
  legacy_competing_convergence_engine_count = sum(grepl(
    "geweke.diag", strsplit(source_text, "\n", fixed = TRUE)[[1]],
    fixed = TRUE)),
  mt_csr_multichain_activation = grepl("nchains", read("R/mtblr-csr.R"), fixed = TRUE),
  mt_block_eigen_multichain_activation = grepl("nchains", read("R/mtblr-block-eigen.R"), fixed = TRUE),
  st_core_convergence_route_coverage = grepl(".blr_finalize_st_public", source_text, fixed = TRUE),
  st_block_eigen_public_activation = "stblr_block_eigen" %in% exports,
  common_fit_field_coverage = all(c("convergence_traces", "memory_estimate") %in%
    strsplit(read("R/blr-unified.R"), "[^A-Za-z0-9_]+")[[1]]),
  operator_reduction_owners = file.exists(file.path(root,
    "tests/testthat/test-blr-operator-reductions.R")),
  raw_schema_versions = !grepl('version[^\n]*=[^\n]*2', source_text),
  mt_csr_prepares_operator_once =
    lengths(regmatches(csr_dispatch, gregexpr("read_and_build_st_ld_csr(",
      csr_dispatch, fixed = TRUE))) == 1L &&
    grepl('Rcpp::_ ["operator_preparations"]=1', mt_cpp, fixed = TRUE),
  mt_block_eigen_prepares_operator_once =
    lengths(regmatches(block_dispatch, gregexpr("build_block_eigen(",
      block_dispatch, fixed = TRUE))) == 1L &&
    grepl('Rcpp::_ ["operator_preparations"]=1', mt_cpp, fixed = TRUE),
  mt_native_logical_chain_dispatch =
    grepl("schedule(static)", mt_cpp, fixed = TRUE) &&
    grepl("std::min(ncores,nchains)", mt_cpp, fixed = TRUE),
  worker_independent_seed_mapping =
    grepl("mt_summary_resolve_chain_seeds", mt_cpp, fixed = TRUE) &&
    grepl("seeds[static_cast<std::size_t>(chain)]", mt_cpp, fixed = TRUE),
  deterministic_mt_aggregation =
    grepl(".mtblr_summary_pool_raw", read("R/mtblr-summary-chains.R"),
          fixed = TRUE),
  scheduled_st_task_traces_before_aggregation =
    all(c("task_vbs", "task_vgs", "task_ves", "task_vle", "task_vld") %in%
      strsplit(scheduled_types, "[^A-Za-z0-9_]+")[[1]]) &&
    grepl('Rcpp::Named("trace")', scheduled_binding, fixed = TRUE),
  scheduled_trace_retention_independent =
    grepl("expect_identical(base$convergence_traces, retained$convergence_traces)",
          unified_tests, fixed = TRUE) &&
    grepl("expect_identical(base$convergence_traces, thinned$convergence_traces)",
          unified_tests, fixed = TRUE),
  all_st_public_routes_use_modern_convergence = all(vapply(
    c("R/stblr-public.R", "R/stblr-block-eigen.R", "R/stblr-csr-annot.R"),
    function(path) grepl(".blr_finalize_st_public", read(path), fixed = TRUE),
    logical(1))),
  executable_operator_reductions = all(vapply(
    c("stblr_csr", "stblr_block_eigen", "mtblr_csr",
      "mtblr_block_eigen", "stblr_bed"),
    grepl, logical(1), x = reduction_source, fixed = TRUE)) &&
    grepl("expect_equal", reduction_tests, fixed = TRUE),
  common_memory_ownership_categories = all(c(
    "shared_immutable_operator_data_bytes",
    "private_sampler_state_per_worker_bytes",
    "result_state_per_logical_chain_bytes",
    "convergence_trace_capture_bytes") %in%
    strsplit(unified_r, "[^A-Za-z0-9_]+")[[1]]),
  six_scientific_models = all(c(
    '"bayesc"', '"sbayesc"', '"bayesr"', '"sbayesr"',
    '"bayesrc"', '"sbayesrc"') %in%
    strsplit(unified_r, "[^A-Za-z0-9_\"]+")[[1]]),
  annotation_policy_inventory = all(c(
    "global", "fixed_marker", "group", "learned_logistic",
    "annotation_probit_stick") %in%
    strsplit(unified_r, "[^A-Za-z0-9_]+")[[1]]),
  s_models_reuse_kernels = grepl(
    'sbayesc = "bayesc", sbayesr = "bayesr"', unified_r, fixed = TRUE),
  exact_capability_matrix = grepl(
    ".blr_model_capability_matrix <- function", unified_r, fixed = TRUE) &&
    grepl("unsupported scientific model/operator combinations fail early",
          read("tests/testthat/test-blr-unified-public-contract.R"),
          fixed = TRUE),
  common_five_trace_contract = all(c(
    '"vbs"', '"vgs"', '"ves"', '"vle"', '"vld"') %in%
    strsplit(convergence_r, "[^A-Za-z0-9_\"]+")[[1]]),
  mt_vle_vld_chain_private = all(vapply(c(
    "result.vle", "result.vld", "diagonal_contribution"), grepl,
    logical(1), x = paste(mt_default_core, mt_bed_core), fixed = TRUE)),
  vld_identity = grepl("vld[t][it]=vgs[t][it]-vle[t][it]",
                       gsub("[[:space:]]", "", mt_default_core), fixed = TRUE) &&
    grepl("result.vld[trait][iteration]=result.vgs[trait][iteration]-result.vle[trait][iteration]",
          gsub("[[:space:]]", "", mt_bed_core), fixed = TRUE),
  common_trace_orientation = grepl(
    "iteration × trait", read("docs/dev/blr_unified_architecture.md"),
    fixed = TRUE),
  common_convergence_quantity_names = !grepl(
    "B_diag|G_diag|E_diag", convergence_r) &&
    grepl('c("vbs", "vgs", "ves", "vle", "vld")', convergence_r,
          fixed = TRUE),
  separate_mt_covariance_outputs = all(c(
    "cov_b_mean", "cov_g_mean", "cov_e_mean", "cov_b_final",
    "cov_g_final", "cov_e_final") %in%
    strsplit(unified_r, "[^A-Za-z0-9_]+")[[1]]),
  common_status_vocabulary = all(c(
    "computed", "computed_fewer_than_four_chains", "computed_partial",
    "not_updated", "not_applicable", "structural_zero", "constant",
    "constant_chain_mismatch", "unavailable_single_chain",
    "insufficient_draws", "nonfinite", "not_requested") %in%
    strsplit(convergence_r, "[^A-Za-z0-9_]+")[[1]]),
  explicit_scale_metadata = all(c(
    "genotype_scale", "effect_scale", "phenotype_scale", "ld_scale") %in%
    strsplit(unified_r, "[^A-Za-z0-9_]+")[[1]]),
  explicit_sample_size_metadata = all(c(
    "n_total", "n_used", "n_by_trait") %in%
    strsplit(unified_r, "[^A-Za-z0-9_]+")[[1]]),
  retention_independence = grepl(
    "convergence_control$keep_traces", read("docs/dev/blr_naming_conventions.md"),
    fixed = TRUE))
for (name in names(values)) cat(name, "=", values[[name]], "\n", sep = "")
guards <- c(values[["canonical_public_fitting_exports"]] == length(canonical),
            values[["redundant_wrapper_count"]] == 0,
            values[["experimental_fitting_export_count"]] == 0,
            values[["st_formatter_count"]] == 1,
            values[["mt_formatter_count"]] == 1,
            values[["shared_convergence_engine_count"]] == 1,
            values[["legacy_competing_convergence_engine_count"]] == 0,
            as.logical(values[c(2:4, 11:length(values))]))
cat("PHASE18_ARCHITECTURE_GUARDS_PASS=", all(guards), "\n", sep = "")
stopifnot(all(guards))
