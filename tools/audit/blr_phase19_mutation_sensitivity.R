root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
pkgload::load_all(root, compile = FALSE, quiet = TRUE)
patterns <- sblr:::.mtblr_models(matrix(c(0L, 1L), 2L, 1L), c(.8, .2),
                                  .2, 1L)
fails <- function(expr, pattern = NULL) {
  value <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  !is.null(value) && (is.null(pattern) || grepl(pattern, value, fixed = TRUE))
}
valid <- sblr:::.mtblr_bayesr_spec("bayesr", patterns, NULL, 2L,
                                    c(0, .1, 1))
semantics_tests <- read("tests/testthat/test-blr-model-semantics.R")
kernel <- read("src/blr_mt_bayesr_kernel_impl.h")
summary_core <- read("src/blr_mt_default_core_impl.h")
bed_core <- read("src/blr_mt_bed_core_impl.h")
binding <- read("src/mtblr.cpp")
bayesr_r <- read("R/mtblr-bayesr.R")
csr_r <- read("R/mtblr-csr.R")
block_r <- read("R/mtblr-block-eigen.R")
bed_r <- read("R/mtblr-bed.R")
workflow <- paste(read(".github/workflows/blr-framework.yml"),
                  read(".github/workflows/blr-framework-extended.yml"))
all_r <- paste(vapply(list.files("R", "[.]R$", full.names = TRUE), read,
                      character(1)), collapse = "\n")
has <- function(text, ...) all(vapply(list(...), grepl, logical(1), x = text,
                                      fixed = TRUE))
guards <- c(
  fails(sblr:::.mtblr_bayesr_spec("bayesr", patterns, NULL, 2L, c(0, 0, 1))),
  fails(sblr:::.mtblr_bayesr_spec("bayesr", patterns, NULL, 2L, c(.1, 1))),
  fails(sblr:::.mtblr_bayesr_spec("bayesr", patterns, NULL, 2L, c(0, 1, .1))),
  identical(valid$joint_names, c("null", "1__component_1", "1__component_2")),
  identical(valid$joint_component, c(0L, 1L, 2L)),
  identical(valid$joint_multiplier, c(0, .1, 1)),
  has(summary_core, "selected>0 ?", "effective") &&
    has(bed_core, "selected>0", "effective"),
  has(bayesr_r, "component must be zero exactly for null marker patterns"),
  has(bayesr_r, "positive multiplier"),
  has(kernel, "spec.multiplier[state]*marker_scale"),
  has(summary_core, "mt_bayesr_base_effects") &&
    has(bed_core, "mt_bayesr_base_effects"),
  has(bayesr_r, "selection_s + 1"),
  identical(sblr:::.mtblr_bayesr_spec(
    "bayesr", patterns, c(.2, .4), 2L, c(0, 1),
    selection_s = -1)$marker_scale, c(1, 1)),
  fails(sblr:::.mtblr_bayesr_spec("sbayesr", patterns, NULL, 2L, c(0, 1))),
  identical(sblr:::.mtblr_bayesr_spec("bayesr", patterns, c(.2, .4), 2L,
    c(0, 1), selection_s = -1)$marker_scale, c(1, 1)),
  abs(sum(valid$patterns$probabilities) - 1) < 1e-15,
  has(summary_core, "cmodel=*model.pi_prior", "samplePi(cmodel, pi"),
  !grepl("component_B", paste(kernel, summary_core, bed_core), fixed = TRUE),
  has(bayesr_r, "rowSums(raw$marker$component_probabilities)" ) ||
    has(csr_r, "rowSums(raw$marker$component_probabilities)"),
  has(summary_core, "dm[t][i]") && has(bed_core, "result.dm[trait][marker]"),
  !grepl("joint_state_probabilities", paste(csr_r, block_r, bed_r), fixed = TRUE),
  has(binding, "operator_preparations", "storage_owners"),
  has(binding, "build_block_eigen(", "operator_preparations"),
  has(binding, "prepare_mt_bed_adapter", "dispatch_mt_bed_chain_tasks"),
  has(binding, "seeds[static_cast<std::size_t>(chain)]"),
  !grepl("rerun", paste(summary_core, bed_core), fixed = TRUE),
  lengths(regmatches(all_r, gregexpr(".blr_convergence_rhat_basic <- function",
    all_r, fixed = TRUE))) == 1L,
  !grepl("pi_pattern", read("R/mtblr-convergence.R"), fixed = TRUE),
  has(csr_r, '.is_mtblr_raw <- function(raw)',
      'identical(as.integer(raw$schema$version), 1L)'),
  has(read("tests/testthat/test-blr-operator-reductions.R"), "BayesC"),
  !grepl('method %in% c("bayesc", "bayesr", "sbayesr", "bayesrc"',
         paste(csr_r, block_r, bed_r), fixed = TRUE),
  has(workflow, "Rscript tools/check/check_package.R ."),
  identical(sblr:::.blr_model_semantics(
    "sbayesr", "csr")$effect_scale_policy, "component"),
  identical(sblr:::.blr_model_semantics(
    "sbayesc", "csr")$effect_scale_policy, "unit"),
  identical(sblr:::.blr_model_semantics(
    "sbayesrc", "csr")$effect_scale_policy, "component"),
  fails(sblr:::.mtblr_resolve_public_method("bayesr", "csr")),
  fails(sblr:::.mtblr_resolve_public_method("sbayesr", "packed_bed")),
  fails(sblr:::.mtblr_resolve_public_method("bayesr", "block_eigen")),
  has(bayesr_r, "selection_s_active <- !is.null(selection_s)"),
  !grepl("selection_s <- 0", bayesr_r, fixed = TRUE),
  has(all_r, "data_level", "prior_kernel", "effect_scale_policy"),
  has(all_r, "model_semantics_version = 2L",
      "s_prefix_means_summary_statistics"),
  has(all_r, "selection_maf_source", "selection_maf_alignment_status",
      "selection_maf_fallback_used"),
  has(all_r, "allow_reference_maf_for_selection_s = FALSE"),
  has(semantics_tests, "cannot be silently reinterpreted"),
  has(semantics_tests, "prior_kernel", "data_level"),
  has(semantics_tests, "selection_s = -1", "selection_s_active"),
  !grepl('method = c("bayesrc"', paste(csr_r, block_r, bed_r), fixed = TRUE),
  !grepl('method = c("sbayesrc"', paste(csr_r, block_r, bed_r), fixed = TRUE))
for (i in seq_along(guards)) cat(sprintf("MUTATION_%02d_DETECTED=%s\n", i,
                                         guards[[i]]))
if (length(guards) != 49L || !all(guards))
  stop("Phase 19 mutation sensitivity failed at: ",
       paste(which(!guards), collapse = ", "), call. = FALSE)
cat("PHASE19_MUTATION_SENSITIVITY=PASS\n")
