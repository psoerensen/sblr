root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
contract <- read("docs/dev/blr_mt_bed_extended_convergence_contract.md")
public <- read("R/mtblr-bed.R")
native <- paste(vapply(c("src/mtblr.cpp", "src/blr_mt_bed_core_impl.h",
  "src/blr_mt_bed_convergence_types.h",
  "src/blr_mt_bed_convergence_trace_impl.h"), read, character(1)),
  collapse = "\n")
has <- function(...) all(vapply(list(...), grepl, logical(1), x = contract,
                                fixed = TRUE))

values <- c(
  CURRENT_TIER2_TRACE_IMPLEMENTATION_COUNT =
    sum(grepl("MtBedExtendedTrace", strsplit(native, "\n", fixed = TRUE)[[1]],
              fixed = TRUE)),
  CURRENT_TIER3_TRACE_IMPLEMENTATION_COUNT =
    sum(grepl("marker_b_values|marker_d_values", native)),
  CURRENT_PUBLIC_EXTENDED_MODE =
    grepl('c("auto", "none", "core", "extended")', public, fixed = TRUE),
  CURRENT_EXTENDED_NATIVE_ROUTE =
    grepl("extended_convergence", read("src/mtblr.cpp"), fixed = TRUE),
  TIER2_COVARIANCE_CONTRACT_EXPLICIT =
    has("Tier 2A", "B/G/E off-diagonal", "raw covariance"),
  TIER2_PROBABILITY_CONTRACT_EXPLICIT =
    has("Tier 2B", "pi_active", "Pattern diagnostics are marginal"),
  TIER3_MARKER_CONTRACT_EXPLICIT =
    has("Tier 3", "effective b", "binary inclusion state"),
  CANONICAL_CHECKPOINT_EXPLICIT =
    has("end of every", "completed Gibbs", "after the marker sweep",
        "accumulator block", "next iteration"),
  LOWER_TRIANGLE_ORDER_EXPLICIT =
    has("strict lower triangle", "column-major order", "Qoff=T*(T-1)/2"),
  NO_DIAGONAL_DUPLICATION =
    has("Tier 2A adds only off-diagonals", "Tier 1 diagonals are never copied"),
  NULL_ACTIVE_DEDUPLICATION_EXPLICIT =
    has("pi_mass:null_active", "not a second summary row", "overview count"),
  PATTERN_SELECTION_EXPLICIT =
    has("public model names", "one-based model-pattern indices",
        "order is preserved"),
  MARKER_SELECTION_EXPLICIT =
    has("Exactly one of `marker_ids` and `marker_indices`",
        "positions in that final selected marker order"),
  BINARY_DIAGNOSTIC_SEMANTICS_EXPLICIT =
    has("all-zero/all-one", "constant_chain_mismatch", "Tail ESS may be unavailable"),
  TRACE_OWNERSHIP_EXPLICIT =
    has("Every live chain exclusively owns", "no shared trace buffer",
        "no R/Rcpp activity"),
  MEMORY_FORMULAS_EXPLICIT =
    has("B = updateB * 8*C*N*Qoff", "d = 4*C*N*K*T", "O(C*N)"),
  LARGE_REQUEST_POLICY_EXPLICIT =
    has("allow_large_traces=TRUE", "never truncated"),
  OUTPUT_COMPATIBILITY_EXPLICIT =
    has("Existing `fit$convergence_traces$scope`", "remain present"),
  WARNING_POLICY_EXPLICIT =
    has("At most one advisory warning", "diagnostic_key"),
  STAGED_IMPLEMENTATION_EXPLICIT =
    has("Phase 17X", "Phase 17Y", "separate contract/implementation phase"))

for (name in names(values)) cat(name, "=", values[[name]], "\n", sep = "")
guards <- c(
  values[["CURRENT_TIER2_TRACE_IMPLEMENTATION_COUNT"]] == 0,
  values[["CURRENT_TIER3_TRACE_IMPLEMENTATION_COUNT"]] == 0,
  !values[["CURRENT_PUBLIC_EXTENDED_MODE"]],
  !values[["CURRENT_EXTENDED_NATIVE_ROUTE"]],
  as.logical(values[5:length(values)]))
cat("ALL_EXTENDED_CONTRACT_GUARDS_PASS=", all(guards), "\n", sep = "")
stopifnot(all(guards))
