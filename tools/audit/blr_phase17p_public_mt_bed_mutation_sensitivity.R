args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
adapter <- read("R/mtblr-bed.R")
namespace <- read("NAMESPACE")
core <- read("src/blr_mt_bed_core_impl.h")

detects <- c(
  BED_HELPER_REUSED = grepl(".make_bed_marker_data(", adapter, fixed = TRUE),
  ONE_NATIVE_CALL = length(gregexpr("raw <- mtblr_bed_internal(", adapter,
    fixed = TRUE)[[1L]]) == 1L,
  RIGHT_BACKEND = !grepl("mtblr_csr_internal|mtblr_block_eigen_internal|mtblr_eigen",
                         adapter),
  CENTER_AFTER_ALIGNMENT = regexpr(".make_bed_marker_data(", adapter,
    fixed = TRUE) < regexpr('if (center)', adapter, fixed = TRUE),
  CENTER_FALSE_VERIFIED = grepl("center = FALSE requires", adapter,
                                fixed = TRUE),
  NO_PHENOTYPE_SCALING = grepl('phenotype_scaling = "not_performed"',
                               adapter, fixed = TRUE),
  COVARIATES_REJECTED = grepl("does not currently fit or project covariates",
                              adapter, fixed = TRUE),
  MISSING_REJECTED = grepl("complete finite matrix", adapter, fixed = TRUE),
  SCALE_FALSE_REJECTED = grepl("requires scale = TRUE", adapter, fixed = TRUE),
  FREQUENCIES_REJECTED = grepl("strictly inside", adapter, fixed = TRUE),
  ONE_GLIST = grepl("accepts one Glist", adapter, fixed = TRUE),
  MARKER_ORDER_PRESERVED = grepl("selected_glist_order_preserved", adapter,
                                 fixed = TRUE),
  SAMPLE_ORDER_PRESERVED = grepl("phenotype_order_preserved", adapter,
                                 fixed = TRUE),
  MODEL_STATE_VALIDATED = grepl("Every state row", adapter, fixed = TRUE),
  INACTIVE_EFFECT_REJECTED = grepl("exactly zero for inactive", adapter,
                                   fixed = TRUE),
  FULL_E_NOT_FORCED_DIAGONAL = grepl("residual_covariance == \"diagonal\"",
                                    adapter, fixed = TRUE),
  DIAGONAL_E_REJECTS_FULL = grepl("must be exactly diagonal", adapter,
                                  fixed = TRUE),
  MEMORY_NOT_PEAK = grepl("not measured peak RSS", adapter, fixed = TRUE),
  SAMPLE_RESIDUAL_NOT_EXPOSED = grepl("sample_residual_returned = FALSE",
                                     adapter, fixed = TRUE),
  GENETIC_VALUES_NOT_EXPOSED = grepl("genetic_values_returned = FALSE",
                                    adapter, fixed = TRUE),
  ONE_FORMATTER = length(gregexpr(".as_mtblr_fit(", adapter,
                                 fixed = TRUE)[[1L]]) == 1L,
  RAW_SCHEMA_UNCHANGED = !grepl("schema_version\\s*<-", adapter),
  PHASE17O_CORE_PRESENT = grepl("run_mt_bed_bayesc_core", core, fixed = TRUE),
  SCALAR_BED_UNTOUCHED = !grepl("stblr_bed\\s*<-", adapter),
  SUMMARY_MT_UNTOUCHED = !grepl("mtblr_csr\\s*<-|mtblr_block_eigen\\s*<-",
                                adapter),
  NO_PARALLEL_CONTROLS = !grepl("ncores|nchains|OpenMP", adapter),
  CPO_UNSUPPORTED = grepl('cpo = "unsupported"', adapter, fixed = TRUE),
  LE_LD_UNSUPPORTED = grepl('le_ld = "unsupported"', adapter, fixed = TRUE),
  NO_LEGACY_EIGEN = !grepl("mtblr_eigen(", adapter, fixed = TRUE),
  PUBLIC_EXPORTED = grepl("export(mtblr_bed)", namespace, fixed = TRUE)
)
for (name in names(detects)) {
  cat(name, "=", toupper(as.character(detects[[name]])), "\n", sep = "")
}
if (!all(detects)) stop("Undetected Phase 17P critical mutation(s): ",
                        paste(names(detects)[!detects], collapse = ", "))
cat("ALL_PHASE17P_CRITICAL_MUTATIONS_DETECTED=TRUE\n")
