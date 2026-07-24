root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
r <- read("R/mtblr-block-eigen.R")
common <- read("R/mtblr-csr.R")
native <- read("src/mtblr.cpp")
diagnostic_adapter <- read("src/st_block_eigen_rcpp.h")
namespace <- read("NAMESPACE")

required_formals <- c(
  "stats", "Glist", "block_start", "operator_sharing", "eigen_filter",
  "eigen_tau", "eigen_eta", "summary_reference", "trait_metadata",
  "marker_policy", "sample_overlap", "method", "n", "sets", "b", "h2",
  "pi", "models", "pimodels", "vg", "vb", "ve", "ssb_prior",
  "sse_prior", "updateB", "updateE", "updatePi", "nub", "nue", "nit",
  "nburn", "nthin", "seed", "verbose")
helper_names <- c(
  ".mtblr_marker_metadata", ".mtblr_normalize_stats",
  ".mtblr_normalize_scale", ".mtblr_align", ".mtblr_models", ".mtblr_sets",
  ".mtblr_cov", ".is_mtblr_raw", ".validate_mtblr_raw", ".as_mtblr_fit")

pkgload::load_all(root, compile = FALSE, quiet = TRUE)
actual_formals <- names(formals(mtblr_block_eigen))
checks <- c(
  PUBLIC_FUNCTION_EXPORTED =
    grepl("export(mtblr_block_eigen)", namespace, fixed = TRUE),
  PUBLIC_FORMALS = identical(actual_formals, required_formals),
  SHARED_HELPERS_PRESENT =
    all(vapply(helper_names, exists, logical(1), envir = asNamespace("sblr"),
               inherits = FALSE)),
  ONE_MT_GIBBS_LOOP =
    length(gregexpr("run_mt_bayesc_core_impl<", native, fixed = TRUE)[[1L]]) == 1L,
  ONE_ACTIVE_MT_MARKER_CONDITIONAL =
    length(gregexpr("sampleBetaCPG_Mt_latent(", native, fixed = TRUE)[[1L]]) >= 1L,
  ONE_MT_LEGACY_TO_RAW_CONVERTER =
    length(gregexpr("Rcpp::List mtblr_legacy_to_raw(", native,
                    fixed = TRUE)[[1L]]) == 1L,
  ONE_GENERAL_MT_FIT_FORMATTER =
    length(gregexpr(".as_mtblr_fit <- function(", common,
                    fixed = TRUE)[[1L]]) == 1L,
  PUBLIC_BLOCK_EIGEN_USES_PHASE17L =
    grepl("mtblr_block_eigen_raw_internal(", r, fixed = TRUE) &&
    grepl("run_mt_block_eigen_adapter(", native, fixed = TRUE),
  EXTERNAL_SUMMARY_REJECTED =
    grepl("External GWAS/reference-panel projection is not yet supported", r,
          fixed = TRUE),
  PROVENANCE_FIELD_COVERAGE =
    all(vapply(c("bed_files", "n_bed", "cls", "rows", "af", "marker_ids",
                 "marker_metadata", "source", "scale", "sample_size"),
               grepl, logical(1), x = r, fixed = TRUE)),
  ALIGNMENT_STATUS_COVERAGE =
    all(vapply(c("bed_provenance_status", "row_provenance_status",
                 "column_provenance_status",
                 "allele_frequency_provenance_status",
                 "wy_transformation_status"),
               grepl, logical(1), x = r, fixed = TRUE)),
  BLOCK_DIAGNOSTIC_COVERAGE =
    all(vapply(c("owner_count", "trait_owner", "n_kept", "mu_min", "shrink"),
               grepl, logical(1), x = paste(r, native, diagnostic_adapter),
               fixed = TRUE)),
  SAMPLE_OVERLAP_POLICY = grepl("sample_overlap must be exactly", r,
                                fixed = TRUE),
  RESIDUAL_COVARIANCE_POLICY =
    grepl('residual_covariance_policy = "diagonal"', r, fixed = TRUE),
  SUMMARY_WW_POLICY =
    grepl("validated_by_construction_not_used_as_runtime_diagonal", r,
          fixed = TRUE))

cat("PUBLIC_FORMALS=", paste(actual_formals, collapse = ","), "\n", sep = "")
cat("SHARED_HELPER_COUNT=", sum(vapply(helper_names, exists, logical(1),
    envir = asNamespace("sblr"), inherits = FALSE)), "\n", sep = "")
cat("RAW_CONVERTER_COUNT=",
    length(gregexpr("Rcpp::List mtblr_legacy_to_raw(", native,
                    fixed = TRUE)[[1L]]), "\n", sep = "")
cat("GENERAL_FORMATTER_COUNT=",
    length(gregexpr(".as_mtblr_fit <- function(", common,
                    fixed = TRUE)[[1L]]), "\n", sep = "")
for (name in names(checks)) cat(name, "=", checks[[name]], "\n", sep = "")
if (!all(checks)) stop("Phase 17M public audit failed: ",
                       paste(names(checks)[!checks], collapse = ", "))
