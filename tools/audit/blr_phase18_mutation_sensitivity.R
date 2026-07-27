root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
all_r <- paste(vapply(list.files("R", "[.]R$", full.names = TRUE), read,
                      character(1)), collapse = "\n")
tests <- paste(vapply(list.files("tests/testthat", "[.]R$", full.names = TRUE),
                      read, character(1)), collapse = "\n")
names_doc <- read("docs/dev/blr_naming_conventions.md")
arch <- read("docs/dev/blr_unified_architecture.md")
workflow <- paste(read(".github/workflows/blr-framework.yml"),
                  read(".github/workflows/blr-framework-extended.yml"))
ns <- read("NAMESPACE")
mt_cpp <- read("src/mtblr.cpp")
scheduled_core <- read("src/blr_csr_scheduled_bayesc_core_impl.h")
scheduled_binding <- read("src/st_cpg_omp_csr_scheduled.cpp")
unified_convergence <- read("tests/testthat/test-blr-unified-convergence.R")
reductions <- read("tests/testthat/test-blr-operator-reductions.R")
reduction_source <- paste(read("tests/testthat/helper-blr-unified.R"), reductions)
unified_r <- read("R/blr-unified.R")
convergence_r <- read("R/mtblr-convergence.R")
mt_default_core <- read("src/blr_mt_default_core_impl.h")
mt_bed_core <- read("src/blr_mt_bed_core_impl.h")
public_contract <- read("tests/testthat/test-blr-unified-public-contract.R")
has <- function(text, ...) all(vapply(list(...), grepl, logical(1), x = text,
                                      fixed = TRUE))
guards <- c(
  has(names_doc, "Methods accept one exact lowercase spelling"),
  has(names_doc, "`cpg` is an internal historical kernel", "never a public model"),
  !has(ns, "export(stblr_csr_bayesr)"),
  !has(ns, "export(stblr_bed_marker)"),
  has(names_doc, "`ncores` always means requested concurrent logical MCMC tasks"),
  has(arch, "Seeds attach to logical tasks before dispatch"),
  has(arch, "unpooled", "post-burn", "unthinned"),
  has(arch, "Compact-chain retention", "independent"),
  has(arch, "unthinned"),
  length(gregexpr(".blr_convergence_scalar <- function", all_r, fixed = TRUE)[[1]]) == 1L,
  length(gregexpr(".blr_convergence_rhat_basic <- function", all_r, fixed = TRUE)[[1]]) == 1L,
  length(gregexpr(".blr_convergence_ess <- function", all_r, fixed = TRUE)[[1]]) == 1L,
  !grepl("geweke[.]diag", all_r),
  has(tests, "fit$convergence"),
  has(tests, "convergence_traces"),
  has(tests, "auto", "nchains"),
  has(arch, "no per-quantity warning"),
  has(names_doc, "analytical memory", "never described as", "measured RSS"),
  has(arch, "shared immutable state"),
  has(arch, "deterministic task/result order"),
  has(tests, 'expect_false(any(c("pi", "pim", "pis", "covb", "vb", "bm_sd")'),
  !grepl("fit$comp_prob <-", all_r, fixed = TRUE),
  has(names_doc, "Backend details live in `diagnostics`"),
  has(names_doc, "cov_b_final", "cov_b_mean"),
  has(arch, "float-packed reconstructed dense blocks", "not low-rank"),
  has(tests, "unsupported", "block_eigen"),
  has(tests, "operator reduction"),
  has(tests, "one-trait", "ST", "MT"),
  has(all_r, '.is_mtblr_raw <- function(raw)',
      'identical(as.integer(raw$schema$version), 1L)'),
  !grepl('"bayesrc"', paste(read("R/mtblr-csr.R"),
    read("R/mtblr-block-eigen.R"), read("R/mtblr-bed.R")), fixed = TRUE),
  !grepl("marker_b_values", all_r, fixed = TRUE),
  has(workflow, "Rscript tools/check/check_package.R ."),
  has(mt_cpp, "mtblr_csr_chains_raw_internal", "operator_preparations"),
  has(mt_cpp, "mtblr_block_eigen_chains_raw_internal", "build_block_eigen("),
  has(mt_cpp, "schedule(static)", "std::min(ncores,nchains)"),
  has(mt_cpp, "mt_summary_resolve_chain_seeds", "seeds[static_cast<std::size_t>(chain)]"),
  has(mt_cpp, "std::vector<sblr::core::SparseLdCsrStorage> storage_owners",
              "std::vector<sblr::core::BlockEigenStorage> storage_owners"),
  has(scheduled_core, "result.task_vbs=std::move(vbs_task)",
                      "result.task_vld=std::move(vlds_task)"),
  has(scheduled_binding, 'Rcpp::Named("trace")',
                         'Rcpp::Named("chains")=chains'),
  has(unified_convergence,
      "expect_identical(base$convergence_traces, retained$convergence_traces)",
      "expect_identical(base$convergence_traces, thinned$convergence_traces)"),
  has(reduction_source, "stblr_csr", "stblr_block_eigen",
                  "mtblr_csr", "mtblr_block_eigen"),
  has(reduction_source, "stblr_bed", "one-trait MT and ST"),
  has(all_r, "shared_immutable_operator_data_bytes",
             "private_sampler_state_per_worker_bytes",
             "result_state_per_logical_chain_bytes"),
  !grepl("operator_preparations\"]=0", mt_cpp, fixed = TRUE),
  has(names_doc, "`bayesc`, `sbayesc`, `bayesr`, `sbayesr`, `bayesrc`, and `sbayesrc`"),
  has(unified_r, 'sbayesc = "bayesc", sbayesr = "bayesr"'),
  has(unified_r, "global", "fixed_marker", "group", "learned_logistic",
      "annotation_probit_stick"),
  has(unified_r, ".blr_model_capability_matrix <- function"),
  has(public_contract, "unsupported scientific model/operator combinations fail early"),
  has(mt_default_core, "diagonal_contribution", "result.vle", "result.vld"),
  has(mt_bed_core, "data.marker_maps[marker].xx", "result.vle", "result.vld"),
  has(mt_default_core, "vld[t][it]=vgs[t][it]-vle[t][it]"),
  has(mt_bed_core, "result.vgs[trait][iteration] - result.vle[trait][iteration]"),
  !grepl("B_diag|G_diag|E_diag", convergence_r),
  has(convergence_r, 'c("vbs", "vgs", "ves", "vle", "vld")'),
  has(unified_r, "cov_b_mean", "cov_b_final"),
  has(convergence_r, "not_applicable", "structural_zero", "not_updated"),
  has(unified_r, "genotype_scale", "effect_scale", "phenotype_scale", "ld_scale"),
  has(unified_r, "n_total", "n_used", "n_by_trait"),
  has(names_doc, "`keep_chains` retains compact logical-chain records",
      "`convergence_control$keep_traces` independently retains convergence arrays"))
for (i in seq_along(guards)) cat(sprintf("MUTATION_%02d=%s\n", i, guards[i]))
cat("ALL_60_CRITICAL_MUTATIONS_DETECTED=", all(guards), "\n", sep = "")
stopifnot(length(guards) == 60L, all(guards))
