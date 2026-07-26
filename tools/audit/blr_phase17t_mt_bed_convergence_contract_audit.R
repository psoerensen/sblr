root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
contract <- read("docs/dev/blr_mt_bed_convergence_contract.md")
public <- read("R/mtblr-bed.R")
engine <- read("R/mtblr-convergence.R")
native <- paste(vapply(c("src/blr_mt_bed_core_impl.h",
  "src/blr_mt_bed_chains_types.h", "src/blr_mt_bed_chains_execution_impl.h",
  "src/blr_mt_bed_chains_aggregate_impl.h", "src/mtblr.cpp"), read,
  character(1)), collapse = "\n")
has <- function(pattern) grepl(pattern, contract, fixed = TRUE)
values <- c(
  CURRENT_RHAT_IMPLEMENTATION_COUNT = sum(grepl("rhat_rank|rhat_folded", c(public, native, engine))),
  CURRENT_ESS_IMPLEMENTATION_COUNT = sum(grepl("ess_bulk|ess_tail", c(public, native, engine))),
  CURRENT_MCSE_IMPLEMENTATION_COUNT = sum(grepl("mcse_mean", c(public, native, engine))),
  CURRENT_PUBLIC_CONVERGENCE_FORMAL = grepl("convergence =", public, fixed = TRUE),
  CURRENT_NATIVE_DIAGNOSTIC_TRACE_BUNDLE = grepl("MtBedConvergenceTraceBundle", native, fixed = TRUE),
  CURRENT_CHAIN_VBS_TRACE = grepl("summary.vbs=core.vbs", native, fixed = TRUE),
  CURRENT_CHAIN_VGS_TRACE = grepl("summary.vgs=core.vgs", native, fixed = TRUE),
  CURRENT_CHAIN_VES_TRACE = grepl("summary.ves=core.ves", native, fixed = TRUE),
  CURRENT_MARKER_ITERATION_TRACES = grepl("marker_iteration_trace", native, fixed = TRUE),
  CURRENT_FULL_COVARIANCE_ITERATION_TRACES = grepl("covariance_iteration_trace", native, fixed = TRUE),
  CURRENT_PI_ITERATION_TRACES = grepl("pi_iteration_trace", native, fixed = TRUE),
  POSTBURN_POLICY_EXPLICIT = has("nburn + 1") && has("nburn + nit"),
  SPLIT_POLICY_EXPLICIT = has("floor(N/2)") && has("discard the central draw"),
  RANK_NORMALIZATION_EXPLICIT = has("qnorm((r - 3/8)/(S + 1/4))"),
  FOLDED_RHAT_EXPLICIT = has("rhat=max(rhat_rank,rhat_folded)"),
  BULK_ESS_EXPLICIT = has("positive sequence") && has("monotone sequence"),
  TAIL_ESS_EXPLICIT = has("draw<=q05") && has("draw<=q95"),
  MCSE_EXPLICIT = has("mcse_mean=posterior_sd/sqrt(ess_mean)"),
  CONSTANT_POLICY_EXPLICIT = has("status `constant`") && has("NA R-hat, ESS, and MCSE"),
  WARNING_POLICY_EXPLICIT = has("at most one main-thread warning per fit"),
  KEEP_CHAINS_INDEPENDENCE_EXPLICIT = has("must work with `keep_chains=FALSE`"),
  MEMORY_POLICY_EXPLICIT = has("tier1_trace_bytes") && has("selected_d_trace_bytes"),
  STAGED_IMPLEMENTATION_EXPLICIT = has("## 33. Phase 17U implementation") && has("## 34. Phase 17V plan")
)
for (name in names(values)) cat(name, "=", values[[name]], "\n", sep = "")
stopifnot(values[["CURRENT_RHAT_IMPLEMENTATION_COUNT"]] == 1,
          values[["CURRENT_ESS_IMPLEMENTATION_COUNT"]] == 1,
          values[["CURRENT_MCSE_IMPLEMENTATION_COUNT"]] == 1,
          !as.logical(values[["CURRENT_PUBLIC_CONVERGENCE_FORMAL"]]),
          as.logical(values[["CURRENT_NATIVE_DIAGNOSTIC_TRACE_BUNDLE"]]),
          all(as.logical(values[c(
            "CURRENT_CHAIN_VBS_TRACE", "CURRENT_CHAIN_VGS_TRACE",
            "CURRENT_CHAIN_VES_TRACE")])),
          !any(as.logical(values[c(
            "CURRENT_MARKER_ITERATION_TRACES",
            "CURRENT_FULL_COVARIANCE_ITERATION_TRACES",
            "CURRENT_PI_ITERATION_TRACES")])),
          all(as.logical(values[12:length(values)])))
