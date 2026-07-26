args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
core <- read("src/blr_mt_bed_core_impl.h")
access <- read("src/blr_mt_bed_access.h")
types <- read("src/blr_mt_bed_types.h")
binding <- read("src/mtblr.cpp")
covariance <- read("src/blr_mt_covariance_rng.h")
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
X <- case$fixture$X
residual <- case$Y - X %*% raw$marker$b

detects <- c(
  PACKED_OWNER = grepl("PackedBedMatrix owner=", binding, fixed = TRUE),
  ONE_OWNER = length(gregexpr("PackedBedMatrix owner=", binding,
                              fixed = TRUE)[[1]]) == 1L,
  COMMON_VIEW = grepl("BedPackedGenotypeView<PackedBedMatrix>", binding,
                      fixed = TRUE),
  MISSING_ZERO = grepl("maps[marker].value[1]=0.0", access, fixed = TRUE),
  DECODED_XX = grepl("xx+=x*x", access, fixed = TRUE),
  SCORE_ADDBACK = grepl("xx * result.b", core, fixed = TRUE),
  FULL_E_DOUBLE_MASK = grepl("D \\* E_inverse \\* D", core),
  FULL_LATENT_DRAW = grepl("standard_normal\\(nt\\)", core),
  EFFECTIVE_MASK = grepl("models[selected][trait]", core, fixed = TRUE),
  RESIDUAL_SIGN = grepl("residual.col(trait) -= marker_workspace",
                        core, fixed = TRUE),
  SAMPLE_SPACE_MCMC = grepl("arma::mat residual = data.phenotype", core,
                            fixed = TRUE),
  FULL_E_UPDATE = grepl("draw_inverse_wishart", core, fixed = TRUE),
  DIAGONAL_DRAW_ORDER = grepl("for \\(std::size_t trait = 0; trait < nt",
                              core),
  BSET = grepl("sampleBset(", core, fixed = TRUE),
  BLATENT = grepl("sampleB_latent(", core, fixed = TRUE),
  BGLOBAL = grepl("sampleB(", core, fixed = TRUE),
  PI_ORDER = grepl("samplePi(cmodel, result.pi", core, fixed = TRUE),
  RETENTION = grepl("marker_retained_count", core, fixed = TRUE),
  MARKER_ORDER = grepl("data.marker_order", core, fixed = TRUE),
  MARKER_R = grepl("result.r[trait][marker]", core, fixed = TRUE) &&
    max(abs(raw$marker$r - crossprod(X, residual))) < 1e-12,
  INDIVIDUAL_LEVEL = identical(raw$meta$data_level, "individual"),
  SHARED_RAW = grepl('legacy, models, "mt_bed_bayesc", "individual"',
                     binding, fixed = TRUE),
  ONE_FINALIZER = !grepl("MtBedFinalResult", paste(core, types)),
  NO_LEGACY_EIGEN = !grepl("mtblr_eigen", core, fixed = TRUE),
  NO_PUBLIC_MTBLR_BED = !grepl("export(mtblr_bed)", namespace, fixed = TRUE),
  NO_OPENMP_MULTICHAIN = !grepl("omp|nchains", core, ignore.case = TRUE),
  NO_CPO = !grepl("log_cpo|sampleCPO", core),
  NO_LE_LD = !grepl("vle|vld", core),
  SUMMARY_CORE_UNTOUCHED = !grepl("sample-space", read(
    "src/blr_mt_default_core_impl.h"), ignore.case = TRUE),
  SCALAR_PACKED_UNTOUCHED = !grepl("MtBed", read(
    "src/st_cpg_omp_individual_scheduled_chains.cpp"), fixed = TRUE),
  ONE_STD_INVERSE_WISHART =
    length(gregexpr("draw_inverse_wishart\\(", covariance, perl = TRUE)[[1]]) ==
      1L
)
for (name in names(detects)) {
  cat(name, "=", toupper(as.character(detects[[name]])), "\n", sep = "")
}
if (!all(detects)) {
  stop("Undetected Phase 17O critical mutation(s): ",
       paste(names(detects)[!detects], collapse = ", "))
}
cat("ALL_PHASE17O_CRITICAL_MUTATIONS_DETECTED=TRUE\n")

