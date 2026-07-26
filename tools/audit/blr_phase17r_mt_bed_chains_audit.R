root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
binding <- read("src/mtblr.cpp")
execution <- read("src/blr_mt_bed_chains_execution_impl.h")
aggregate <- read("src/blr_mt_bed_chains_aggregate_impl.h")
public <- read("R/mtblr-bed.R")
route <- sub(".*Rcpp::List mtblr_bed_chains_internal", "", binding)
route <- sub("// INTERNAL RESEARCH ONLY.*", "", route)
prepared <- sub(".*MtBedPreparedAdapter prepare_mt_bed_adapter", "", binding)
prepared <- sub("Rcpp::List mt_bed_marker_kernel_to_list.*", "", prepared)
count_fixed <- function(pattern, text) lengths(regmatches(text, gregexpr(pattern, text, fixed = TRUE)))
guards <- c(
  packed_owner_count = count_fixed("prepared.owner.reset(new PackedBedMatrix", binding) == 1L,
  borrowed_view_count = grepl("genotype_view() const", binding, fixed = TRUE),
  BED_read_count = count_fixed("read_bedfiles_to_packed_matrix(", prepared) == 1L,
  task_count_formula = grepl("MtBedChainTask> tasks(", route, fixed = TRUE),
  result_slot_formula = grepl("results[chain]", execution, fixed = TRUE),
  shared_phenotype_count = grepl("prepared.phenotype", route, fixed = TRUE),
  shared_marker_map_count = grepl("prepared.marker_maps", route, fixed = TRUE),
  shared_wy_count = grepl("prepared.marker_wy", route, fixed = TRUE),
  shared_order_count = grepl("prepared.marker_order", route, fixed = TRUE),
  private_residual_count = grepl("run_mt_bed_bayesc_core", execution, fixed = TRUE),
  private_workspace_count = grepl("run_mt_bed_bayesc_core", execution, fixed = TRUE),
  private_rng_count = grepl("execution.seed=task.seed", execution, fixed = TRUE),
  OpenMP_loop_count = count_fixed("#pragma omp parallel for", execution) == 1L,
  OpenMP_schedule = grepl("schedule(static)", execution, fixed = TRUE),
  R_API_inside_workers = !grepl("R::|Rf_|Rcpp::", execution),
  Rcpp_inside_workers = !grepl("Rcpp::", execution, fixed = TRUE),
  worker_print_count = !grepl("cout|cerr|Rprintf|Rcout", execution),
  failure_capture_count = grepl("catch (const std::exception&", execution, fixed = TRUE) &&
    grepl("catch (...)", execution, fixed = TRUE),
  aggregation_count = count_fixed("aggregate_mt_bed_chains(", route) == 1L,
  finalizer_count = count_fixed("finalize_mt_bed_chains_result(", route) == 1L &&
    count_fixed("finalize_mt_default_result(", aggregate) == 1L,
  legacy_adapter_count = count_fixed("make_mt_default_legacy_result(", route) == 1L,
  raw_converter_count = count_fixed("mtblr_legacy_to_raw(", route) == 1L,
  public_serial_route = grepl("mtblr_bed_internal(", public, fixed = TRUE) &&
    !grepl("mtblr_bed_chains_internal", public, fixed = TRUE),
  pooled_counts = grepl("marker_retained_count+=", aggregate, fixed = TRUE),
  sample_sd = grepl("summaries.size()-1", aggregate, fixed = TRUE)
)
for (name in names(guards)) cat(name, "=", as.integer(guards[[name]]), "\n", sep = "")
cat("task_count=nchains\nresult_slot=chain\nprivate_residual=nchains_logically\n")
if (!all(guards)) stop("Phase 17R ownership/threading audit failed")
