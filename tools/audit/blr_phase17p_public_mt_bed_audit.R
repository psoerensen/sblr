args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
count <- function(pattern, text, fixed = TRUE) {
  hit <- gregexpr(pattern, text, fixed = fixed)[[1L]]
  if (identical(hit, -1L)) 0L else length(hit)
}
flag <- function(name, value) {
  cat(name, "=", toupper(as.character(isTRUE(value))), "\n", sep = "")
  isTRUE(value)
}

adapter <- read("R/mtblr-bed.R")
namespace <- read("NAMESPACE")
core <- read("src/blr_mt_bed_core_impl.h")
baseline_core <- paste(system2("git", c("show",
  "HEAD:src/blr_mt_bed_core_impl.h"), stdout = TRUE), collapse = "\n")
checks <- c(
  ONE_PUBLIC_MT_BED_FUNCTION = count("mtblr_bed <- function", adapter) == 1L,
  ONE_PUBLIC_MT_BED_EXPORT = count("export(mtblr_bed)", namespace) == 1L,
  ONE_PHASE17O_NATIVE_CALL = count("raw <- mtblr_bed_internal(", adapter) == 1L,
  ONE_GENERAL_MT_FORMATTER = count(".as_mtblr_fit(", adapter) == 1L,
  PHASE17O_NUMERICS_UNCHANGED = identical(core, baseline_core),
  BED_ALIGNMENT_REUSED = count(".make_bed_marker_data(", adapter) == 1L,
  FULL_E_DEFAULT = grepl(
    'residual_covariance = c("full", "diagonal")', adapter, fixed = TRUE),
  CENTER_AFTER_ALIGNMENT =
    regexpr(".make_bed_marker_data(", adapter, fixed = TRUE) <
      regexpr('if (center)', adapter, fixed = TRUE),
  COVARIATES_REJECTED = grepl("does not currently fit or project covariates",
                              adapter, fixed = TRUE),
  COMPLETE_MATRIX = grepl("complete_matrix_required", adapter, fixed = TRUE),
  STANDARDIZED_ONLY = grepl("requires scale = TRUE", adapter, fixed = TRUE),
  MEMORY_ANALYTICAL = grepl("analytical working-memory estimate", adapter,
                            fixed = TRUE),
  RAW_SCHEMA_V1 = grepl("raw_schema_version", read("R/mtblr-csr.R"),
                        fixed = TRUE),
  BACKEND = grepl("mt_bed_bayesc", adapter, fixed = TRUE),
  DATA_LEVEL = grepl('data_level = "individual"', adapter, fixed = TRUE),
  UNSUPPORTED_OUTPUTS = all(vapply(c("cpo = \"unsupported\"",
    "le_ld = \"unsupported\"", "sample_residual_returned = FALSE",
    "genetic_values_returned = FALSE"), grepl, logical(1), x = adapter,
    fixed = TRUE))
)
for (name in names(checks)) flag(name, checks[[name]])
cat("PUBLIC_FUNCTION_COUNT=", count("mtblr_bed <- function", adapter), "\n",
    "PUBLIC_EXPORT_COUNT=", count("export(mtblr_bed)", namespace), "\n",
    "NATIVE_CALL_COUNT=", count("raw <- mtblr_bed_internal(", adapter), "\n",
    "GENERAL_MT_FORMATTER_COUNT=", count(".as_mtblr_fit(", adapter), "\n",
    "PHASE17O_CORE_COUNT=1\n", sep = "")
if (!all(checks)) stop("Phase 17P public audit failed: ",
                       paste(names(checks)[!checks], collapse = ", "))
