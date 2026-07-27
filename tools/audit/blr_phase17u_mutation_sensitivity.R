root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
engine <- read("R/mtblr-convergence.R")
cpp <- read("src/mtblr.cpp")
extractor <- read("src/blr_mt_bed_convergence_trace_impl.h")
execution <- read("src/blr_mt_bed_chains_execution_impl.h")
public <- read("R/mtblr-bed.R")
contract <- read("docs/dev/blr_mt_bed_convergence_contract.md")
description <- read("DESCRIPTION")
namespace <- read("NAMESPACE")
workflow <- read(".github/workflows/blr-framework.yml")
has <- function(text, ...) {
  all(vapply(c(...), function(token) grepl(token, text, fixed = TRUE),
             logical(1)))
}
unchanged <- function(path) {
  identical(system2("git", c("diff", "--quiet", "--", path)), 0L)
}
guards <- c(
  has(extractor, "nburn+iteration"),
  has(cpp, "build_mt_bed_convergence_trace_bundle(\n   results"),
  has(engine, "additional_thinning = FALSE"),
  has(engine, "nrow(x) - half + 1L"),
  has(engine, 'ties.method = "average"'),
  has(engine, "ranks - 3 / 8"),
  has(engine, "rhat_folded"),
  has(engine, "max(rhat_rank, rhat_folded)"),
  has(engine, ".blr_convergence_split_chains(draws)"),
  has(engine, "variance * (n - 1) / n"),
  has(engine, "stats::var(colMeans(x))"),
  has(engine, "rho_even + rho_odd > 0"),
  has(engine, "rho[lag + 1L] + rho[lag + 2L] >"),
  has(engine, "ess_q05", "ess_q95", "min(ess_q05, ess_q95)"),
  !has(engine, "ess_tail <- ess_bulk"),
  has(engine, "mean_ess <- .blr_convergence_ess(split)"),
  has(engine, "posterior_sd / sqrt(ess_mean)"),
  has(engine, '"constant"'),
  has(engine, '"constant_chain_mismatch"'),
  has(engine, 'if (!updated) "not_updated"'),
  has(engine, "E_diag", "updateE"),
  has(engine, "nchains < 2L"),
  has(engine, '"computed_partial"'),
  has(cpp, "capture_convergence_traces"),
  has(cpp, "results, nburn, nit, updateB, updateE"),
  !grepl("B_lower|G_lower|E_lower", extractor),
  !grepl("marker_id|marker_effect", extractor),
  !grepl("pi_trace|model_probability", extractor),
  has(engine, "burnin_included = FALSE"),
  !grepl("nthin", extractor, fixed = TRUE),
  lengths(regmatches(cpp, gregexpr(
    "dispatch_mt_bed_chain_tasks\\(", cpp))) == 1L,
  lengths(regmatches(cpp, gregexpr(
    "prepare_mt_bed_adapter\\(", cpp))) == 3L,
  !grepl("rhat|ess|mcse", execution, ignore.case = TRUE),
  !grepl("convergence =", public, fixed = TRUE),
  !grepl("mtblr_bed_convergence_trace_internal", public, fixed = TRUE),
  !grepl("fit\\$convergence", public),
  !grepl("convergence warning", public, ignore.case = TRUE),
  !grepl("posterior::", engine, fixed = TRUE),
  !grepl("posterior", description, ignore.case = TRUE),
  unchanged("NAMESPACE"),
  has(contract, "raw version 1"),
  unchanged("src/blr_mt_bed_core_impl.h"),
  unchanged("src/blr_mt_bed_chains_aggregate_impl.h"),
  has(workflow, "Rscript tools/check/check_package.R .")
)
names(guards) <- sprintf("MUTATION_%02d", seq_along(guards))
for (name in names(guards)) cat(name, "=", guards[[name]], "\n", sep = "")
cat("ALL_44_CRITICAL_MUTATIONS_DETECTED=", all(guards), "\n", sep = "")
stopifnot(length(guards) == 44L, all(guards))
