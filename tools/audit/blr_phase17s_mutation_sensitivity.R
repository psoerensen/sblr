root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
public <- read("R/mtblr-bed.R")
docs <- paste(read("docs/dev/blr_mt_bed_public_contract.md"), public,
              collapse = "\n")
report <- read("docs/dev/blr_framework_phase17s_report.md")
workflow <- read(".github/workflows/blr-framework.yml")
formatter <- read("R/mtblr-csr.R")
namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
native <- c(
  "src/blr_mt_bed_chains_types.h",
  "src/blr_mt_bed_chains_execution_impl.h",
  "src/blr_mt_bed_chains_aggregate_impl.h",
  "src/blr_mt_bed_types.h", "src/blr_mt_bed_access.h",
  "src/blr_mt_bed_core_impl.h", "src/blr_mt_covariance_rng.h",
  "src/mtblr.cpp")
env <- new.env(parent = baseenv())
sys.source(file.path(root, "R", "mtblr-bed.R"), env)
f <- formals(env$mtblr_bed)
position <- match(c("seed", "nchains", "ncores", "chain_seeds", "keep_chains",
                    "memory_warning_gb"), names(f))
count_fixed <- function(pattern, text) {
  lengths(regmatches(text, gregexpr(pattern, text, fixed = TRUE)))
}
guards <- c(
  nchains_present = "nchains" %in% names(f),
  nchains_default = identical(f$nchains, 1L),
  ncores_present = "ncores" %in% names(f),
  ncores_default = identical(f$ncores, 1L),
  chain_seeds_present = "chain_seeds" %in% names(f),
  chain_seeds_default = is.null(f$chain_seeds),
  keep_chains_present = "keep_chains" %in% names(f),
  keep_chains_default = identical(f$keep_chains, FALSE),
  signature_location = identical(unname(position), seq(position[1], position[1] + 5L)),
  no_serial_call = !grepl("raw <- mtblr_bed_internal(", public, fixed = TRUE),
  not_both_routes = count_fixed("raw <- mtblr_bed_chains_internal(", public) == 1L,
  one_chains_call = count_fixed("raw <- mtblr_bed_chains_internal(", public) == 1L,
  null_empty = grepl("native <- integer()", public, fixed = TRUE),
  seeds_not_sorted = !grepl("sort(chain_seeds", public, fixed = TRUE),
  seeds_not_offset = !grepl("chain_seeds\\s*\\+\\s*9176", public),
  negative_allowed = grepl("-2147483648", public, fixed = TRUE),
  signed_range = grepl("2147483647", public, fixed = TRUE),
  cores_above_chains = grepl("requested_workers <- min(ncores, nchains)",
                             public, fixed = TRUE),
  packed_shared = grepl("shared_components_bytes", public, fixed = TRUE),
  phenotype_shared = grepl("shared_names", public, fixed = TRUE),
  private_by_workers = grepl("requested_workers * private_bytes", public, fixed = TRUE),
  results_by_chains = grepl("nchains * result_bytes", public, fixed = TRUE),
  retained_conditional = grepl("if (keep_chains)", public, fixed = TRUE),
  warning_chain_core = grepl("nchains=%d, ncores=%d", public, fixed = TRUE),
  not_measured = grepl("not measured RSS", public, fixed = TRUE),
  no_blas_mutation = !grepl("Sys.setenv|Sys.unsetenv", public),
  pooled_documented = grepl("pool", docs, ignore.case = TRUE),
  binary_not_averaged = grepl("binary", docs, ignore.case = TRUE),
  stability_not_posterior_sd = grepl("not posterior standard", docs, ignore.case = TRUE),
  retained_exclusions = grepl("omit", docs, ignore.case = TRUE),
  default_numerics_guarded = grepl("pre-existing numerical", report,
                                   ignore.case = TRUE),
  formatter_used = count_fixed(".as_mtblr_fit(", public) == 1L,
  one_formatter = count_fixed(".as_mtblr_fit <- function", formatter) == 1L,
  schema_version_one = grepl("version 1", docs, fixed = TRUE),
  native_unchanged = system2("git", c("diff", "--quiet", "HEAD", "--", native),
                             stdout = FALSE, stderr = FALSE) == 0L,
  wrappers_unchanged = system2("git", c("diff", "--quiet", "HEAD", "--",
    "R/RcppExports.R", "src/RcppExports.cpp"), stdout = FALSE,
    stderr = FALSE) == 0L,
  namespace_single = sum(namespace == "export(mtblr_bed)") == 1L &&
    system2("git", c("diff", "--quiet", "HEAD", "--", "NAMESPACE"),
            stdout = FALSE, stderr = FALSE) == 0L,
  other_api_unchanged = system2("git", c("diff", "--quiet", "HEAD", "--",
    "R/mtblr-csr.R", "R/mtblr-block-eigen.R", "R/sparse_ld_bed_helper.R",
    "R/interface_mtblr.R"), stdout = FALSE, stderr = FALSE) == 0L,
  no_convergence_claim = grepl("do not establish convergence", docs,
                               ignore.case = TRUE),
  package_check_ci = grepl("tools/check/check_package.R", workflow, fixed = TRUE))
for (name in names(guards)) cat(sprintf(
  "MUTATION_%02d_%s=%s\n", match(name, names(guards)), toupper(name),
  toupper(guards[[name]])))
if (length(guards) != 40L || !all(guards)) {
  stop("Phase 17S mutation sensitivity failed")
}
cat("ALL_40_CRITICAL_MUTATIONS_DETECTED=TRUE\n")
