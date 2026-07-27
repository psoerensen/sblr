root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
has <- function(text, ...) all(vapply(list(...), grepl, logical(1), x = text,
                                      fixed = TRUE))
count <- function(text, pattern) {
  hit <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  if (identical(hit, -1L)) 0L else length(hit)
}
bayesrc <- read("R/mtblr-bayesrc.R")
resolver <- read("R/mtblr-bayesr.R")
csr <- read("R/mtblr-csr.R")
block <- read("R/mtblr-block-eigen.R")
bed <- read("R/mtblr-bed.R")
types <- read("src/blr_mt_bayesrc_types.h")
summary_core <- read("src/blr_mt_default_core_impl.h")
bed_core <- read("src/blr_mt_bed_core_impl.h")
binding <- read("src/mtblr.cpp")
convergence <- read("R/mtblr-convergence.R")
tests <- paste(read("tests/testthat/test-mtblr-bayesrc-model.R"),
               read("tests/testthat/test-mtblr-bayesrc-operators.R"),
               read("tests/testthat/test-mtblr-bayesrc-annotations.R"))
guards <- c(
  ONE_ANNOTATION_COMPONENT_MODEL = count(types, "struct MtBayesRCSpec") == 1L,
  ONE_PROBIT_STICK_CONVENTION = has(types, "st_bayesrc_annotation_prior.h") &&
    count(bayesrc, ".mtblr_bayesrc_prior_probabilities <- function") == 1L,
  ONE_UNIQUE_NULL_STATE = has(types, "prior[0]=theta[0]", "state=1" ) ||
    has(types, "prior[0]=theta[0]", "for (std::size_t state=1"),
  ONE_PATTERN_VECTOR = has(types, "std::vector<double> pattern_prior",
                            "const std::vector<double>& omega"),
  PHASE19_STATE_DESCRIPTOR = has(types, "MtJointStateSpec") &&
    has(summary_core, "sample_mt_bayesr_marker") &&
    has(bed_core, "mt_joint_marker_kernel"),
  ONE_BASE_B_SCALING = has(summary_core, "mt_bayesr_base_effects") &&
    has(bed_core, "mt_bayesr_base_effects"),
  NO_ANNOTATION_B_MATRICES = !grepl("annotation_B", paste(types, summary_core,
    bed_core), fixed = TRUE),
  NO_PATTERN_ANNOTATION_COEFFICIENTS = !grepl("pattern_alpha", paste(types,
    summary_core, bed_core), fixed = TRUE),
  BED_BAYESRC = has(resolver, 'packed_bed = c("bayesc", "bayesr", "bayesrc")'),
  CSR_SBAYESRC = has(resolver, 'csr = c("sbayesc", "sbayesr", "sbayesrc")'),
  BLOCK_SBAYESRC = has(resolver,
    'block_eigen = c("sbayesc", "sbayesr", "sbayesrc")'),
  SHARED_ANNOTATION_PREPARATION = count(bayesrc,
    ".stblr_align_bed_bayesrc_annotations(") == 1L,
  SHARED_OPERATOR_PREPARATION = has(binding, "operator_preparations",
    "prepare_mt_bed_adapter"),
  NATIVE_DETERMINISTIC_CHAINS = has(binding, "schedule(static)",
    "std::min(ncores,nchains)"),
  WORKER_INDEPENDENT_SEEDS = has(binding, "mt_summary_resolve_chain_seeds",
    "seeds[static_cast<std::size_t>(chain)]"),
  RAW_SCHEMA_VERSION_ONE = has(csr, "mtblr_raw", "version"),
  MODEL_SEMANTICS_VERSION_TWO = has(csr, "model_semantics_version = 2L") &&
    has(block, "model_semantics_version = 2L") &&
    has(bed, "model_semantics_version = 2L"),
  SELECTION_S_INDEPENDENT = has(resolver,
    "selection_s_active <- !is.null(selection_s)") &&
    !grepl("selection_s <- 0", paste(resolver, bayesrc), fixed = TRUE),
  CORE_CONVERGENCE_UNCHANGED = count(convergence,
    ".blr_convergence_rhat_basic <- function") == 1L &&
    count(convergence, ".blr_convergence_ess <- function") == 1L,
  EXPLICIT_MEMORY_OWNERSHIP = has(bayesrc,
    "bayesrc_shared_annotation_bytes",
    "bayesrc_private_annotation_bytes_per_worker",
    "bayesrc_chain_annotation_result_bytes",
    "bayesrc_formatted_annotation_output_bytes"),
  PERMANENT_REDUCTIONS = has(tests, "reduces exactly to Phase 19 SBayesR",
    "CSR and exact block eigen reduce numerically"),
  EXPLICIT_ALIGNMENT = has(bayesrc,
    "external row order is never trusted", "exact_marker_id_match"),
  PRIOR_POSTERIOR_PROBABILITIES_DISTINCT = has(bayesrc,
    "prior_component_probabilities", "prior_active_probabilities") &&
    has(csr, "component_probabilities"),
  PI_NOT_REDEFINED = has(bayesrc, 'fit["pi_final"] <- list(NULL)',
    "pattern_pi_final"),
  ANNOTATION_NOT_IN_CORE_CONVERGENCE = !grepl("annotation_coefficients",
    convergence, fixed = TRUE)
)
for (name in names(guards)) cat(name, "=", guards[[name]], "\n", sep = "")
if (!all(guards)) stop("Phase 20 architecture audit failed: ",
  paste(names(guards)[!guards], collapse = ", "), call. = FALSE)
cat("PHASE20_ARCHITECTURE_AUDIT=PASS\n")
