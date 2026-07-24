root <- normalizePath(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ".", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
contract <- read("docs/dev/blr_mt_bed_internal_contract.md")
helper <- read("tests/testthat/helper-mtblr-bed-contract.R")
tests <- read("tests/testthat/test-blr-framework-phase17n.R")
namespace <- read("NAMESPACE")
packed <- read("src/packed_bed.h")
bayesc <- read("src/st_cpg_omp_individual_scheduled_chains.cpp")
family <- read("src/blr_bed_family_types.h")

probes <- c(
  BED_CODE_MAPPING = grepl("c(0L, 1L, 2L, 3L)", tests, fixed = TRUE),
  STANDARDIZED_MISSING_ZERO = grepl("out[is.na(out[, j]), j] <- 0", helper, fixed = TRUE),
  RAW_MISSING_2P = grepl("out[is.na(out[, j]), j] <- 2 * af[j]", helper, fixed = TRUE),
  ROW_ORDER = grepl("rows <- c(7L, 2L, 5L, 1L, 6L)", helper, fixed = TRUE),
  MARKER_ORDER = grepl("cls <- list(c(4L, 2L, 1L), c(3L, 1L))", helper, fixed = TRUE),
  PARTIAL_BYTE = grepl("jbase + 3 < n", bayesc, fixed = TRUE) &&
    grepl("partial-byte", contract, fixed = TRUE),
  NO_PER_CHAIN_PACKED_COPY = grepl("zero per-marker or per-chain genotype copies", contract, fixed = TRUE),
  NO_MCMC_BED_READ = grepl("zero\\s+MCMC-time BED reads", contract, perl = TRUE),
  SCALAR_NOT_JOINT = grepl("do not\\s+fit a joint multivariate likelihood", contract, perl = TRUE),
  COVARIATES_UNSUPPORTED = grepl("pre-adjusted phenotypes", contract, fixed = TRUE),
  MISSING_PHENOTYPES_UNSUPPORTED = grepl("complete finite phenotype matrix only", contract, fixed = TRUE),
  SAMPLE_MARKER_RESIDUAL_DISTINCT = grepl("r_jt = x_j'R_t", contract, fixed = TRUE) &&
    grepl("R = Y - X B_eff", contract, fixed = TRUE),
  FULL_E_PRESENT = grepl("C_k   = P + w_j D_k Omega D_k", contract, fixed = TRUE),
  OUTPUT_SCHEMA_DECIDED = grepl("reuses `mtblr_raw` version 1", contract, fixed = TRUE),
  NO_INDIVIDUAL_MT_SAMPLER = !file.exists(file.path(root, "src/blr_mt_bed_core_impl.h")),
  NO_PUBLIC_MTBLR_BED = !grepl("export(mtblr_bed)", namespace, fixed = TRUE),
  NO_LEGACY_EIGEN_BASE = grepl("must not use\\s+legacy `mtblr_eigen\\(\\)`", contract, perl = TRUE),
  PHASE15B_VIEW_RETAINED = grepl("reuses `BedPackedGenotypeView` unchanged", contract, fixed = TRUE) &&
    grepl("struct BedPackedGenotypeView", family, fixed = TRUE)
)
for (name in names(probes)) {
  cat(name, "=", toupper(as.character(probes[[name]])), "\n", sep = "")
}
if (!all(probes)) {
  stop("Undetected Phase 17N critical mutation(s): ",
       paste(names(probes)[!probes], collapse = ", "))
}
cat("ALL_PHASE17N_CRITICAL_MUTATIONS_DETECTED=TRUE\n")
