root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
cpp <- read("src/mtblr.cpp")
types <- read("src/blr_mt_bed_convergence_types.h")
extractor <- read("src/blr_mt_bed_convergence_trace_impl.h")
engine <- read("R/mtblr-convergence.R")
public <- read("R/mtblr-bed.R")
description <- read("DESCRIPTION")
workflow <- read(".github/workflows/blr-framework.yml")

diff_clean <- function(path) {
  identical(system2("git", c("diff", "--quiet", "--", path)), 0L)
}
count_fixed <- function(text, token) {
  hit <- gregexpr(token, text, fixed = TRUE)[[1L]]
  if (length(hit) == 1L && hit[1L] < 0L) 0L else length(hit)
}

values <- c(
  internal_convergence_route_count =
    count_fixed(cpp, "Rcpp::List mtblr_bed_convergence_trace_internal("),
  trace_bundle_type_count =
    count_fixed(types, "struct MtBedConvergenceTraceBundle"),
  trace_extractor_count =
    count_fixed(extractor, "build_mt_bed_convergence_trace_bundle("),
  production_rhat_engine_count =
    count_fixed(engine, ".blr_convergence_rhat_basic <- function"),
  production_ess_engine_count =
    count_fixed(engine, ".blr_convergence_ess <- function"),
  production_mcse_engine_count =
    count_fixed(engine, "mcse_mean <- posterior_sd / sqrt(ess_mean)"),
  public_convergence_formal_count =
    count_fixed(public, "convergence ="),
  public_convergence_call_count =
    count_fixed(public, "mtblr_bed_convergence_trace_internal"),
  postburn_only =
    grepl("nburn+iteration", extractor, fixed = TRUE) &&
      !grepl("bundle.values.push_back(trace[iteration])", extractor,
             fixed = TRUE),
  pooled_trace_rejected =
    grepl("results, nburn, nit", cpp, fixed = TRUE) &&
      !grepl("aggregate.pooled, nburn, nit", cpp, fixed = TRUE),
  keep_chains_independent =
    grepl("capture_convergence_traces", cpp, fixed = TRUE) &&
      grepl("results, nburn, nit, updateB, updateE", cpp, fixed = TRUE),
  tier1_only =
    all(vapply(c("B_diag", "G_diag", "E_diag"),
               grepl, logical(1), x = extractor, fixed = TRUE)) &&
      !grepl("marker_effect|model_probability|lower_triangle", extractor),
  dependency_free =
    !grepl("posterior::|coda::|matrixStats::|checkmate::", engine),
  posterior_runtime_dependency =
    grepl("posterior", description, ignore.case = TRUE),
  schema_version =
    grepl("mtblr_convergence_trace_bundle", cpp, fixed = TRUE) &&
      grepl("Rcpp::_\\[\"version\"\\]=1", cpp),
  native_core_hash_unchanged = diff_clean("src/blr_mt_bed_core_impl.h"),
  public_adapter_hash_unchanged = diff_clean("R/mtblr-bed.R"),
  package_check_retained =
    grepl("Rscript tools/check/check_package.R .", workflow, fixed = TRUE)
)
for (name in names(values)) cat(name, "=", values[[name]], "\n", sep = "")

guards <- c(
  ONE_INTERNAL_CONVERGENCE_ROUTE =
    values[["internal_convergence_route_count"]] == 1,
  ONE_TRACE_BUNDLE = values[["trace_bundle_type_count"]] == 1,
  ONE_RHAT_ENGINE = values[["production_rhat_engine_count"]] == 1,
  ONE_ESS_ENGINE = values[["production_ess_engine_count"]] == 1,
  ONE_MCSE_ENGINE = values[["production_mcse_engine_count"]] == 1,
  ZERO_PUBLIC_CONVERGENCE_FORMALS =
    values[["public_convergence_formal_count"]] == 0,
  ZERO_PUBLIC_CONVERGENCE_CALLS =
    values[["public_convergence_call_count"]] == 0,
  KEEP_CHAINS_INDEPENDENT = as.logical(values[["keep_chains_independent"]]),
  TIER1_ONLY = as.logical(values[["tier1_only"]]),
  NO_POSTERIOR_RUNTIME_DEPENDENCY =
    !as.logical(values[["posterior_runtime_dependency"]]),
  PUBLIC_ADAPTER_UNCHANGED =
    as.logical(values[["public_adapter_hash_unchanged"]]),
  NUMERICAL_CORE_UNCHANGED =
    as.logical(values[["native_core_hash_unchanged"]])
)
for (name in names(guards)) cat(name, "=", guards[[name]], "\n", sep = "")
stopifnot(all(guards), values[["trace_extractor_count"]] == 1,
          as.logical(values[["dependency_free"]]),
          as.logical(values[["postburn_only"]]),
          as.logical(values[["pooled_trace_rejected"]]),
          as.logical(values[["schema_version"]]),
          as.logical(values[["package_check_retained"]]))
