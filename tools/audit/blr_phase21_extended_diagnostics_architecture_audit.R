root <- normalizePath(if (file.exists("DESCRIPTION")) "." else "../..",
                      winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
r_unified <- read("R/blr-unified.R")
r_extended <- read("R/blr-extended-convergence.R")
r_convergence <- read("R/mtblr-convergence.R")
cpp_default <- read("src/blr_mt_default_core_impl.h")
cpp_bed <- read("src/blr_mt_bed_core_impl.h")
cpp_st_probability <- paste(read("src/blr_csr_bayesr_core_impl.h"),
  read("src/blr_bed_bayesr_core_impl.h"), collapse = "\n")
cpp_st_selected <- paste(vapply(c(
  "src/blr_csr_bayesc_core_impl.h", "src/blr_csr_prior_bayesc_core_impl.h",
  "src/blr_csr_scheduled_bayesc_core_impl.h",
  "src/blr_bed_scheduled_bayesc_core_impl.h",
  "src/blr_csr_bayesr_core_impl.h", "src/blr_bed_bayesr_core_impl.h",
  "src/blr_csr_sbayesrc_core_impl.h", "src/blr_bed_bayesrc_core_impl.h"),
  read, character(1)), collapse = "\n")
cpp_st_annotations <- paste(read("src/blr_csr_group_bayesc_core_impl.h"),
  read("src/blr_csr_learned_annotation_bayesc_core_impl.h"),
  read("src/blr_csr_sbayesrc_core_impl.h"), collapse = "\n")
public <- paste(vapply(c("R/stblr-public.R", "R/stblr-block-eigen.R",
  "R/stblr-csr-annot.R", "R/mtblr-csr.R", "R/mtblr-block-eigen.R",
  "R/mtblr-bed.R"), read, character(1)), collapse = "\n")
count_definition <- function(text, pattern) length(gregexpr(
  pattern, text, fixed = TRUE)[[1L]][gregexpr(pattern, text, fixed = TRUE)[[1L]] > 0L])
guards <- c(
  ONE_RHAT_IMPLEMENTATION = count_definition(r_convergence,
    ".blr_convergence_rhat_basic <- function") == 1L,
  ONE_ESS_IMPLEMENTATION = count_definition(r_convergence,
    ".blr_convergence_ess <- function") == 1L,
  ONE_MCSE_MEAN_IMPLEMENTATION = count_definition(r_convergence,
    "mcse_mean <- posterior_sd / sqrt(ess_mean)") == 1L,
  ONE_GENERIC_BUNDLE = grepl("blr_convergence_trace_bundle", r_unified,
                             fixed = TRUE),
  AUTO_CORE_ONLY = grepl("identical(convergence, \"extended\")", r_convergence,
                         fixed = TRUE),
  EXTENDED_MODE_ALL_PUBLIC = length(gregexpr(
    'c("auto", "none", "core", "extended")', public, fixed = TRUE)[[1L]]) >= 7L,
  EXTENDED_GROUPS_EXPLICIT = grepl(".blr_extended_group_values", r_extended,
                                   fixed = TRUE),
  STRICT_LOWER_ORDER = grepl("for (int col=0;col<nt;++col)", cpp_default,
                             fixed = TRUE) &&
    grepl("for (int row=col+1;row<nt;++row", cpp_default, fixed = TRUE),
  NO_COVARIANCE_DIAGONAL_DUPLICATION = grepl("nt*(nt-1)/2", cpp_default,
                                             fixed = TRUE),
  BINARY_COMPLEMENT_DEDUPLICATED = grepl("pi_active", r_extended, fixed = TRUE),
  FULL_JOINT_OPT_IN = grepl("full_probability_states", cpp_default,
                            fixed = TRUE),
  SELECTED_MARKERS_EXPLICIT = grepl("all-marker shortcuts are not supported",
                                    r_extended, fixed = TRUE),
  NO_ALL_MARKER_DEFAULT = !grepl('selected_markers = "all"', public,
                                 fixed = TRUE),
  DIRECT_INDEXED_CAPTURE = grepl("trace_spec->selected_markers", cpp_default,
                                 fixed = TRUE) &&
    grepl("convergence->selected_markers", cpp_bed, fixed = TRUE),
  CHAIN_PRIVATE_CAPTURE = grepl("MtExtendedTraceResult extended", cpp_default,
                                fixed = TRUE),
  RNG_NEUTRAL_CAPTURE = !grepl("sample_uniform", paste(
    regmatches(cpp_default, gregexpr("extended[^;]*", cpp_default))[[1L]],
    collapse = "\n"), fixed = TRUE),
  KEEP_CHAINS_INDEPENDENT = !grepl("keep_chains", cpp_default, fixed = TRUE),
  NTHIN_INDEPENDENT = grepl("extended.selected_b", cpp_default, fixed = TRUE),
  HARD_MEMORY_GUARD = grepl("allow_large_traces", r_extended, fixed = TRUE),
  NO_SILENT_TRUNCATION = grepl("proceed without truncation", r_extended,
                               fixed = TRUE),
  RAW_SCHEMA_VERSION_1 = grepl("version = 1L", r_unified, fixed = TRUE),
  MODEL_SEMANTICS_VERSION_2 = grepl("model_semantics_version = 2L", r_unified,
                                    fixed = TRUE),
  NO_NEW_MODEL_OPERATOR = !grepl("phase21.*kernel", cpp_default,
                                 ignore.case = TRUE),
  MT_COVARIANCE_CAPTURE = grepl("extended.cov_b", cpp_default, fixed = TRUE),
  MT_PROBABILITY_CAPTURE = grepl("extended.component_pi", cpp_default,
                                 fixed = TRUE),
  MT_ANNOTATION_CAPTURE = grepl("extended.annotation_alpha", cpp_default,
                                fixed = TRUE),
  MT_SELECTED_CAPTURE = grepl("extended.selected_component", cpp_default,
                              fixed = TRUE),
  ST_PI_CAPTURE_ADAPTER = grepl("pi_active", r_extended, fixed = TRUE),
  ST_SELECTION_S_ADAPTER = grepl("selection_s", r_extended, fixed = TRUE),
  ST_COMPONENT_PI_NATIVE_CAPTURE = grepl("convergence_pi_task", cpp_st_probability,
    fixed = TRUE) && grepl("component_pi", r_extended, fixed = TRUE),
  ST_SELECTED_NATIVE_CAPTURE = grepl("convergence_markers", cpp_st_selected,
    fixed = TRUE) && grepl("selected_component", r_extended, fixed = TRUE),
  ST_GROUP_LEARNED_NATIVE_CAPTURE =
    grepl("convergence_group_pi", cpp_st_annotations, fixed = TRUE) &&
    grepl("convergence_eta_pi", cpp_st_annotations, fixed = TRUE) &&
    grepl("convergence_alpha_task", cpp_st_annotations, fixed = TRUE))
for (name in names(guards)) cat(name, "=", guards[[name]], "\n", sep = "")
cat("PHASE21_ARCHITECTURE_GUARDS_TRUE=", sum(guards), "\n", sep = "")
cat("PHASE21_ARCHITECTURE_GUARDS_TOTAL=", length(guards), "\n", sep = "")
cat("PHASE21_IMPLEMENTATION_READY=", all(guards), "\n", sep = "")
if (!all(guards)) {
  stop("A completed Phase 21 architecture guard regressed.", call. = FALSE)
}
