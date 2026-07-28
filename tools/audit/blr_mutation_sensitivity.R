guard <- function(spec) {
  identical(spec$fitters, c("stblr_csr", "stblr_csr_annot", "stblr_block_eigen",
    "stblr_bed", "mtblr_csr", "mtblr_block_eigen", "mtblr_bed")) &&
  setequal(spec$models, c("bayesc", "sbayesc", "bayesr", "sbayesr", "bayesrc", "sbayesrc")) &&
  identical(spec$s_prefix, "summary_statistics") &&
  identical(spec$maf_effect_s, "independent") &&
  identical(spec$raw_schema, 1L) && identical(spec$model_semantics, 2L) &&
  identical(spec$convergence_engine_count, 1L) &&
  identical(spec$selected_b, "effective") && identical(spec$selected_d, "binary") &&
  identical(spec$all_marker_default, FALSE) && identical(spec$package_check, TRUE) &&
  identical(spec$probability_fields, c("pi_final", "pi_mean", "pi_trace")) &&
  identical(spec$operator_preparations, 1L) &&
  identical(spec$annotation_preparations, 1L) &&
  identical(spec$worker_seed_policy, "logical_task") &&
  identical(spec$trace_checkpoint, "post_burn_every_iteration") &&
  identical(spec$memory_guard, "hard_pre_execution") &&
  identical(spec$annotation_alignment, "exact_marker_id") &&
  identical(spec$generated_interfaces, "source_owned")
}
base <- list(
  fitters = c("stblr_csr", "stblr_csr_annot", "stblr_block_eigen", "stblr_bed",
    "mtblr_csr", "mtblr_block_eigen", "mtblr_bed"),
  models = c("bayesc", "sbayesc", "bayesr", "sbayesr", "bayesrc", "sbayesrc"),
  s_prefix = "summary_statistics", maf_effect_s = "independent",
  raw_schema = 1L, model_semantics = 2L, convergence_engine_count = 1L,
  selected_b = "effective", selected_d = "binary", all_marker_default = FALSE,
  package_check = TRUE,
  probability_fields = c("pi_final", "pi_mean", "pi_trace"),
  operator_preparations = 1L, annotation_preparations = 1L,
  worker_seed_policy = "logical_task",
  trace_checkpoint = "post_burn_every_iteration",
  memory_guard = "hard_pre_execution",
  annotation_alignment = "exact_marker_id",
  generated_interfaces = "source_owned")
mutations <- list(
  removed_fitter = function(x) { x$fitters <- x$fitters[-1L]; x },
  alias_fitter = function(x) { x$fitters[[1L]] <- "sblr"; x },
  mixed_case_model = function(x) { x$models[[1L]] <- "BayesC"; x },
  s_means_maf = function(x) { x$s_prefix <- "maf_s"; x },
  implicit_maf_effect_s = function(x) { x$maf_effect_s <- "model_name"; x },
  raw_schema_bump = function(x) { x$raw_schema <- 2L; x },
  semantic_version_removed = function(x) { x$model_semantics <- NULL; x },
  second_convergence_engine = function(x) { x$convergence_engine_count <- 2L; x },
  selected_b_latent = function(x) { x$selected_b <- "latent"; x },
  selected_d_component = function(x) { x$selected_d <- "component"; x },
  all_marker_default = function(x) { x$all_marker_default <- TRUE; x },
  probability_alias = function(x) { x$probability_fields[[1L]] <- "pi"; x },
  operator_rebuilt_per_chain = function(x) { x$operator_preparations <- 2L; x },
  annotations_rebuilt_per_chain = function(x) { x$annotation_preparations <- 2L; x },
  worker_dependent_seed = function(x) { x$worker_seed_policy <- "worker"; x },
  thinned_diagnostic_trace = function(x) { x$trace_checkpoint <- "retained_only"; x },
  soft_memory_guard = function(x) { x$memory_guard <- "truncate"; x },
  positional_annotation_alignment = function(x) { x$annotation_alignment <- "row_order"; x },
  stale_generated_registration = function(x) { x$generated_interfaces <- "stale"; x },
  package_check_removed = function(x) { x$package_check <- FALSE; x })
if (!guard(base)) stop("base architecture rejected", call. = FALSE)
detected <- vapply(mutations, function(mutate) !guard(mutate(base)), logical(1))
if (!all(detected)) stop("undetected mutations: ",
  paste(names(detected)[!detected], collapse = ", "), call. = FALSE)
cat(sprintf("mutation_sensitivity=%d/%d PASS\n", sum(detected), length(detected)))
