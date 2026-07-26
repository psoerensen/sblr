root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
public <- read("R/mtblr-bed.R")
formatter <- read("R/mtblr-csr.R")
namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
native <- c(
  "src/blr_mt_bed_chains_types.h",
  "src/blr_mt_bed_chains_execution_impl.h",
  "src/blr_mt_bed_chains_aggregate_impl.h",
  "src/blr_mt_bed_types.h", "src/blr_mt_bed_access.h",
  "src/blr_mt_bed_core_impl.h", "src/blr_mt_covariance_rng.h",
  "src/mtblr.cpp", "R/RcppExports.R", "src/RcppExports.cpp")
native_unchanged <- system2("git", c("diff", "--quiet", "HEAD", "--", native),
                            stdout = FALSE, stderr = FALSE) == 0L
env <- new.env(parent = baseenv())
sys.source(file.path(root, "R", "mtblr-bed.R"), env)
f <- formals(env$mtblr_bed)
count_fixed <- function(pattern, text) {
  lengths(regmatches(text, gregexpr(pattern, text, fixed = TRUE)))
}
report <- c(
  public_mt_bed_export_count = sum(namespace == "export(mtblr_bed)"),
  public_nchains_formal = "nchains" %in% names(f),
  public_ncores_formal = "ncores" %in% names(f),
  public_chain_seeds_formal = "chain_seeds" %in% names(f),
  public_keep_chains_formal = "keep_chains" %in% names(f),
  public_default_nchains = identical(f$nchains, 1L),
  public_default_ncores = identical(f$ncores, 1L),
  public_default_chain_seeds = is.null(f$chain_seeds),
  public_default_keep_chains = identical(f$keep_chains, FALSE),
  public_chains_native_call_count = count_fixed(
    "raw <- mtblr_bed_chains_internal(", public),
  public_serial_native_call_count = count_fixed(
    "raw <- mtblr_bed_internal(", public),
  memory_helper_count = count_fixed(
    ".mtblr_bed_memory_estimate <- function", public),
  general_mt_formatter_count = count_fixed(
    ".as_mtblr_fit <- function", formatter),
  raw_schema_version = grepl("mtblr_raw version 1", public, fixed = TRUE) ||
    grepl("raw_schema_version", formatter, fixed = TRUE),
  backend = grepl('backend = "mt_bed_bayesc"', public, fixed = TRUE),
  data_level = grepl('data_level = "individual"', public, fixed = TRUE),
  native_phase17r_hashes_unchanged = native_unchanged)
for (name in names(report)) cat(name, "=", report[[name]], "\n", sep = "")
guards <- c(
  ONE_PUBLIC_MT_BED_EXPORT = report[["public_mt_bed_export_count"]] == 1,
  PUBLIC_MULTICHAIN_FORMALS_ACTIVE = all(as.logical(report[c(
    "public_nchains_formal", "public_ncores_formal",
    "public_chain_seeds_formal", "public_keep_chains_formal",
    "public_default_nchains", "public_default_ncores",
    "public_default_chain_seeds", "public_default_keep_chains")])),
  ONE_CHAINS_NATIVE_CALL = report[["public_chains_native_call_count"]] == 1,
  ZERO_SERIAL_NATIVE_CALLS_FROM_PUBLIC =
    report[["public_serial_native_call_count"]] == 0,
  ONE_MEMORY_HELPER = report[["memory_helper_count"]] == 1,
  ONE_GENERAL_MT_FORMATTER = report[["general_mt_formatter_count"]] == 1,
  PHASE17R_NATIVE_UNCHANGED = native_unchanged,
  RAW_VERSION_ONE = isTRUE(as.logical(report[["raw_schema_version"]])))
for (name in names(guards)) cat(name, "=", toupper(guards[[name]]), "\n", sep = "")
if (!all(guards)) stop("Phase 17S public activation audit failed")

