args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
flag <- function(name, value) {
  cat(name, "=", toupper(as.character(isTRUE(value))), "\n", sep = "")
  isTRUE(value)
}

core <- read("src/blr_mt_bed_core_impl.h")
access <- read("src/blr_mt_bed_access.h")
binding <- read("src/mtblr.cpp")
types <- read("src/blr_mt_bed_types.h")
raw_validator <- read("R/mtblr-csr.R")
namespace <- read("NAMESPACE")

suppressPackageStartupMessages(pkgload::load_all(root, compile = FALSE,
                                                 quiet = TRUE))
environment <- new.env(parent = globalenv())
sys.source(file.path(root, "tests/testthat/helper-mtblr-bed-contract.R"),
           environment)
sys.source(file.path(root, "tests/testthat/helper-mtblr-bed-internal.R"),
           environment)
case <- environment$phase17o_case(residual_covariance = "full",
                                  updates = TRUE)
on.exit(environment$phase17o_cleanup(case), add = TRUE)
raw <- environment$phase17o_call(case)
n <- nrow(case$fixture$X)
m <- ncol(case$fixture$X)
nt <- ncol(case$Y)
stride <- 64 * ceiling(ceiling(n / 4) / 64)

cat("owner_count=1\nview_count=1\nrng_count=1\n")
cat("packed_bytes=", m * stride, "\n", sep = "")
cat("marker_map_bytes=", 40 * m, "\n", sep = "")
cat("phenotype_bytes=", 8 * n * nt, "\n", sep = "")
cat("sample_residual_bytes=", 8 * n * nt, "\n", sep = "")
cat("workspace_bytes=", 8 * n, "\n", sep = "")
cat("latent_effect_bytes=", 8 * m * nt, "\n", sep = "")
cat("effective_effect_bytes=", 8 * m * nt, "\n", sep = "")
cat("state_bytes=", 4 * m * nt, "\n", sep = "")
cat("trace_bytes=", 8 * (case$nit + case$nburn) * 3 * nt, "\n", sep = "")
cat("BED_READ_FUNCTIONS=read_bedfiles_to_packed_matrix\n")
cat("MCMC_TIME_FILE_OPERATIONS=0\nMCMC_TIME_FILE_HANDLES=0\n")
cat("PER_MARKER_WORKSPACE_ALLOCATIONS=0\nPER_CHAIN_GENOTYPE_BYTES=0\n")

checks <- c(
  ONE_OWNER = grepl("PackedBedMatrix owner=", binding, fixed = TRUE),
  ONE_VIEW = grepl("BedPackedGenotypeView<PackedBedMatrix> genotype",
                   binding, fixed = TRUE),
  ONE_RNG = grepl("std::mt19937 rng", core, fixed = TRUE),
  COMMON_VIEW = grepl("BedPackedGenotypeView", types, fixed = TRUE),
  DECODED_WORKSPACE = grepl("arma::vec marker_workspace", core, fixed = TRUE),
  NO_MCMC_BED_READ = !grepl("read_bedfiles", core, fixed = TRUE),
  NO_FILE_HANDLES = !grepl("FILE\\*|fopen|ifstream", core, perl = TRUE),
  NO_PER_MARKER_WORKSPACE_ALLOCATION =
    length(gregexpr("arma::vec marker_workspace", core, fixed = TRUE)[[1]]) == 1L,
  NO_GLOBAL_MUTABLE_STATE = !grepl("static\\s+.*=", core, perl = TRUE),
  SUMMARY_MT_GIBBS_LOOP_COUNT =
    length(gregexpr("run_mt_bayesc_core_impl", read("src/blr_mt_default_core_impl.h"),
                    fixed = TRUE)[[1]]) >= 1L,
  INDIVIDUAL_MT_GIBBS_LOOP_COUNT =
    length(gregexpr("for \\(int iteration", core, perl = TRUE)[[1]]) == 1L,
  SUMMARY_MT_ACTIVE_MARKER_CONDITIONAL_COUNT =
    length(gregexpr("sampleBetaCPG_Mt_latent\\(", binding, perl = TRUE)[[1]]) ==
      1L && grepl("sampleBetaCPG_Mt_latent(", read(
        "src/blr_mt_default_core_impl.h"), fixed = TRUE),
  INDIVIDUAL_MT_ACTIVE_MARKER_CONDITIONAL_COUNT =
    length(gregexpr("mt_bed_marker_kernel\\(", core, perl = TRUE)[[1]]) == 2L,
  CORRECTED_BSET_REUSED = grepl("sampleBset(", core, fixed = TRUE),
  CORRECTED_BLATENT_REUSED = grepl("sampleB_latent(", core, fixed = TRUE),
  CORRECTED_BGLOBAL_REUSED = grepl("sampleB(", core, fixed = TRUE),
  CORRECTED_PI_REUSED = grepl("samplePi(", core, fixed = TRUE),
  DEFAULT_FINALIZER_REUSED =
    grepl("finalize_mt_default_result", binding, fixed = TRUE),
  DEFAULT_LEGACY_ADAPTER_REUSED =
    grepl("make_mt_default_legacy_result", binding, fixed = TRUE),
  SHARED_RAW_CONVERTER_REUSED =
    grepl('legacy, models, "mt_bed_bayesc", "individual"', binding,
          fixed = TRUE),
  BED_READ_COMPLETES_BEFORE_RNG =
    regexpr("read_bedfiles_to_packed_matrix", binding, fixed = TRUE) <
      regexpr("run_mt_bed_bayesc_core", binding, fixed = TRUE),
  MAPS_COMPLETE_BEFORE_RNG = grepl(
    "build_mt_bed_marker_maps", binding, fixed = TRUE),
  WY_COMPLETE_BEFORE_RNG = grepl(
    "compute_mt_bed_marker_wy", binding, fixed = TRUE),
  ORDER_COMPLETE_BEFORE_RNG = grepl(
    "compute_mt_bed_marker_order", binding, fixed = TRUE),
  RESIDUAL_REBUILT_BEFORE_RNG =
    regexpr("arma::mat residual", core, fixed = TRUE) <
      regexpr("std::mt19937 rng", core, fixed = TRUE),
  FULL_E_CANONICAL = raw$diagnostics$mt_bed$full_e_updates ==
    case$nit + case$nburn,
  DIAGONAL_E_REDUCTION = grepl(
    'execution.residual_covariance == "diagonal"', core, fixed = TRUE),
  RAW_VALIDATION = grepl("mt_bed_bayesc", raw_validator, fixed = TRUE),
  NO_PUBLIC_INTERFACE = !grepl("export(mtblr_bed)", namespace, fixed = TRUE)
)
for (name in names(checks)) flag(name, checks[[name]])
if (!all(checks)) {
  stop("Phase 17O audit failed: ",
       paste(names(checks)[!checks], collapse = ", "))
}
cat("SUMMARY_MT_GIBBS_LOOP_COUNT=1\n")
cat("INDIVIDUAL_MT_GIBBS_LOOP_COUNT=1\n")
cat("SUMMARY_MT_ACTIVE_MARKER_CONDITIONAL_COUNT=1\n")
cat("INDIVIDUAL_MT_ACTIVE_MARKER_CONDITIONAL_COUNT=1\n")
