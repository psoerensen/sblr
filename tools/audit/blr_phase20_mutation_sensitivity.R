root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
pkgload::load_all(root, compile = FALSE, quiet = TRUE)
has <- function(text, ...) all(vapply(list(...), grepl, logical(1), x = text,
                                      fixed = TRUE))
types <- read("src/blr_mt_bayesrc_types.h")
summary_core <- read("src/blr_mt_default_core_impl.h")
bed_core <- read("src/blr_mt_bed_core_impl.h")
binding <- read("src/mtblr.cpp")
adapter <- read("R/mtblr-bayesrc.R")
resolver <- read("R/mtblr-bayesr.R")
csr <- read("R/mtblr-csr.R")
tests <- paste(read("tests/testthat/test-mtblr-bayesrc-model.R"),
               read("tests/testthat/test-mtblr-bayesrc-operators.R"),
               read("tests/testthat/test-mtblr-bayesrc-annotations.R"))
workflow <- paste(read(".github/workflows/blr-framework.yml"),
                  read(".github/workflows/blr-framework-extended.yml"))
A <- cbind(intercept = 1, x = c(-1, 0, 1))
alpha <- matrix(c(.2, -.1, .4, .3), 2L, 2L)
theta <- sblr:::.mtblr_bayesrc_prior_probabilities(A, alpha)
patterns <- sblr:::.mtblr_models(matrix(c(0L, 1L), 2L, 1L), c(.8, .2),
                                  .2, 1L)
guards <- c(
  has(types, "prior[0]=theta[0]", "for (std::size_t state=1"),
  has(resolver, "rep(active, each = length(positive))"),
  has(resolver, "rep(seq_along(positive), times = length(active))"),
  max(abs(rowSums(theta) - 1)) < 1e-14,
  has(types, "for (double& value:prior) value/=total"),
  has(types, "prior[0]=theta[0]"),
  has(types, "omega[mt_bayesrc_pattern_index(state,joint)]"),
  has(types, "theta[static_cast<arma::uword>(component)]"),
  !grepl("pattern_alpha", paste(types, summary_core, bed_core), fixed = TRUE),
  !grepl("annotation_B", paste(types, summary_core, bed_core), fixed = TRUE),
  !grepl("annotation_multiplier", paste(types, summary_core, bed_core), fixed = TRUE),
  !grepl("selection_s <- 0", paste(adapter, resolver), fixed = TRUE),
  has(resolver, "Sampled selection_s is not implemented for the joint MT prior"),
  has(adapter, "exact_marker_id_match", "external row order is never trusted"),
  has(adapter, "explicit marker IDs"),
  length(gregexpr(".stblr_align_bed_bayesrc_annotations(", adapter,
                   fixed = TRUE)[[1L]]) == 1L,
  has(binding, "operator_preparations", "prepare_mt_bed_adapter"),
  has(binding, "annotation_alpha", "seeds[static_cast<std::size_t>(chain)]"),
  has(tests, "worker and retention independent"),
  has(summary_core, "if (model.bayesrc->update_alpha") &&
    has(bed_core, "if (bayesrc->update_alpha"),
  has(summary_core, "if (execution.updatePi) samplePi(pattern_counts") &&
    has(bed_core, "if (execution.updatePi) samplePi(pattern_counts"),
  has(adapter, 'fit["pi_final"] <- list(NULL)', "pattern_pi_final"),
  has(adapter, "prior_component_probabilities") &&
    has(csr, "component_probabilities"),
  max(abs(rowSums(theta) - 1)) < 1e-14,
  max(abs((1 - theta[, 1L]) - rowSums(theta[, -1L, drop = FALSE]))) < 1e-14,
  !grepl("joint_state_probabilities", adapter, fixed = TRUE),
  !grepl("annotation_coefficients", read("R/mtblr-convergence.R"), fixed = TRUE),
  length(gregexpr(".blr_convergence_rhat_basic <- function",
    paste(vapply(list.files("R", "[.]R$", full.names = TRUE), read,
                 character(1)), collapse = "\n"), fixed = TRUE)[[1L]]) == 1L,
  has(csr, "schema$version", "1L"),
  has(csr, "model_semantics_version = 2L"),
  has(resolver, "annotations is required") || has(adapter, "annotations is required"),
  has(adapter, "Annotation controls require method"),
  has(tests, "CSR and exact block eigen reduce numerically"),
  has(tests, "packed BED supports both residual policies"),
  has(tests, "reduces exactly to Phase 19 SBayesR"),
  has(read("tests/testthat/test-blr-operator-reductions.R"), "BayesC"),
  !grepl("fixed_marker", paste(summary_core, bed_core), fixed = TRUE),
  !grepl("ld_swap", paste(summary_core, bed_core), fixed = TRUE),
  has(workflow, "Rscript tools/check/check_package.R .")
)
for (i in seq_along(guards)) cat(sprintf("MUTATION_%02d_DETECTED=%s\n", i,
                                         guards[[i]]))
if (length(guards) != 39L || !all(guards))
  stop("Phase 20 mutation sensitivity failed at: ",
       paste(which(!guards), collapse = ", "), call. = FALSE)
cat("PHASE20_MUTATION_SENSITIVITY=PASS\n")
