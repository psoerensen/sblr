root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
bayesr_r <- read("R/mtblr-bayesr.R")
unified_r <- read("R/blr-unified.R")
csr_r <- read("R/mtblr-csr.R")
block_r <- read("R/mtblr-block-eigen.R")
bed_r <- read("R/mtblr-bed.R")
kernel <- read("src/blr_mt_bayesr_kernel_impl.h")
types <- read("src/blr_mt_bayesr_types.h")
summary_core <- read("src/blr_mt_default_core_impl.h")
bed_core <- read("src/blr_mt_bed_core_impl.h")
binding <- read("src/mtblr.cpp")
tests <- paste(read("tests/testthat/test-mtblr-bayesr-model.R"),
               read("tests/testthat/test-mtblr-bayesr-operators.R"),
               read("tests/testthat/test-blr-model-semantics.R"))
all_r <- paste(vapply(list.files(file.path(root, "R"), "[.]R$", full.names = TRUE),
  function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
  character(1)), collapse = "\n")
has <- function(text, ...) all(vapply(list(...), grepl, logical(1), x = text,
                                      fixed = TRUE))
count <- function(text, pattern) {
  hit <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  if (identical(hit, -1L)) 0L else length(hit)
}
guards <- c(
  ONE_JOINT_STATE_SPECIFICATION = count(types, "struct MtJointStateSpec") == 1L,
  DETERMINISTIC_STATE_ORDERING = has(bayesr_r, 'state_names <- c("null"',
                                      "rep(active, each = length(positive))"),
  ONE_UNIQUE_NULL_STATE = has(bayesr_r, "length(null) != 1L",
                               "state_component <- c(0L"),
  ONE_COMPONENT_MULTIPLIER_SPEC = has(bayesr_r, "one leading zero",
                                      "unique ascending positive"),
  ONE_PROBABILITY_NORMALIZATION_UTILITY = count(kernel,
    "joint probabilities cannot be normalized") == 2L,
  ONE_BASE_B_SCALING_CONTRACT = count(kernel,
    "mt_bayesr_base_effects") == 1L,
  NO_COMPONENT_SPECIFIC_B_MATRICES = !grepl("component_B", paste(kernel,
    summary_core, bed_core), fixed = TRUE),
  ALL_MT_OPERATORS_BAYESR = has(bayesr_r,
    'packed_bed = c("bayesc", "bayesr")',
    'csr = c("sbayesc", "sbayesr")',
    'block_eigen = c("sbayesc", "sbayesr")'),
  BAYESC_RETAINED = has(csr_r, '"sbayesc"') &&
    has(block_r, '"sbayesc"') && has(bed_r, '"bayesc"'),
  NO_BAYESRC_PHASE19 = all(!vapply(list(csr_r, block_r, bed_r),
    function(x) grepl('method %in% c("bayesc", "bayesr", "sbayesr", "bayesrc"',
                      x, fixed = TRUE), logical(1))),
  SHARED_OPERATOR_PREPARATION = has(binding, "operator_preparations", "schedule(static)",
    "std::min(ncores,nchains)"),
  WORKER_INDEPENDENT_SEEDS = has(binding, "mt_summary_resolve_chain_seeds",
    "seeds[static_cast<std::size_t>(chain)]"),
  COMMON_OUTPUT_FIELDS = has(csr_r, "component_final", "component_probabilities",
    "pi_trace", "model_parameters"),
  COMMON_CONVERGENCE_FIELDS = has(tests, "vbs[T1]", "vld[T1]",
    "convergence_traces"),
  RAW_SCHEMA_VERSION_ONE = has(csr_r,
    '.is_mtblr_raw <- function(raw)',
    'identical(as.integer(raw$schema$version), 1L)'),
  MEMORY_OWNERSHIP_CATEGORIES = has(bayesr_r,
    "bayesr_shared_state_descriptors_bytes",
    "bayesr_private_worker_state_bytes", "bayesr_chain_result_bytes",
    "bayesr_component_output_bytes"),
  PERMANENT_REDUCTION_OWNERS = has(tests, "exact block eigen",
    "unit scale in every operator"),
  S_PREFIX_DATA_LEVEL = has(bayesr_r,
    'csr = c("sbayesc", "sbayesr")',
    'packed_bed = c("bayesc", "bayesr")'),
  PRIOR_KERNEL_SHARED = has(unified_r,
    'sbayesc = "bayesc", sbayesr = "bayesr"'),
  SELECTION_S_INDEPENDENT = has(bayesr_r,
    "selection_s_active <- !is.null(selection_s)") &&
    !grepl('selection_s <- 0', bayesr_r, fixed = TRUE),
  MODEL_SEMANTICS_VERSION_TWO = has(csr_r,
    "model_semantics_version", "s_prefix_means_summary_statistics"),
  SELECTION_MAF_PROVENANCE = has(bayesr_r,
    "selection_maf_source", "selection_maf_alignment_status",
    "selection_maf_fallback_used"),
  ONE_RHAT_ENGINE = count(all_r, ".blr_convergence_rhat_basic <- function") == 1L,
  ONE_ESS_ENGINE = count(all_r, ".blr_convergence_ess <- function") == 1L)
for (name in names(guards)) cat(name, "=", guards[[name]], "\n", sep = "")
if (!all(guards)) stop("Phase 19 architecture audit failed: ",
                       paste(names(guards)[!guards], collapse = ", "), call. = FALSE)
cat("PHASE19_ARCHITECTURE_AUDIT=PASS\n")
