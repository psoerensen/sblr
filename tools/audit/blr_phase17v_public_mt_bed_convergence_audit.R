root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pkgload::load_all(root, compile = FALSE, quiet = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
                             collapse = "\n")
public <- read("R/mtblr-bed.R")
public_lines <- readLines(file.path(root, "R/mtblr-bed.R"), warn = FALSE)
engine <- read("R/mtblr-convergence.R")
formatter <- read("R/mtblr-csr.R")
namespace <- read("NAMESPACE")
values <- c(
  public_convergence_formal = grepl("convergence = c(\"auto\", \"none\", \"core\")", public, fixed = TRUE),
  public_convergence_control_formal = grepl("convergence_control = NULL", public, fixed = TRUE),
  public_convergence_default = identical(formals(get(
    "mtblr_bed", asNamespace("sblr")))$convergence,
                                         quote(c("auto", "none", "core"))),
  public_convergence_control_default = is.null(formals(get(
    "mtblr_bed", asNamespace("sblr")))$convergence_control),
  public_chains_route_count = sum(trimws(public_lines) ==
    "mtblr_bed_chains_internal"),
  public_trace_route_count = sum(trimws(public_lines) ==
    "mtblr_bed_convergence_trace_internal"),
  conditional_route_count = grepl("native_route <- if", public, fixed = TRUE),
  public_convergence_formatter_count = grepl("fit$convergence <- raw$diagnostics$convergence", formatter, fixed = TRUE),
  public_convergence_warning_builder_count = grepl(".mtblr_convergence_warning_messages <-", engine, fixed = TRUE),
  public_convergence_warning_call_count = grepl("warning(raw$diagnostics$convergence$warning_messages[1L]", public, fixed = TRUE),
  convergence_memory_helper_count = grepl(".mtblr_convergence_memory_estimate <-", engine, fixed = TRUE),
  raw_schema_version = grepl("schema version 1", formatter, fixed = TRUE),
  native_phase17u_hashes_unchanged = system2("git", c("diff", "--quiet", "--", "src")) == 0,
  wrapper_hashes_unchanged = system2("git", c("diff", "--quiet", "--", "R/RcppExports.R", "src/RcppExports.cpp")) == 0,
  namespace_unchanged = system2("git", c("diff", "--quiet", "--", "NAMESPACE")) == 0)
for (name in names(values)) cat(name, "=", values[[name]], "\n", sep = "")
guards <- c(
  PUBLIC_CONVERGENCE_FORMALS_ACTIVE = all(values[c(
    "public_convergence_formal", "public_convergence_control_formal")]),
  PUBLIC_DEFAULT_AUTO = all(values[c(
    "public_convergence_default", "public_convergence_control_default")]),
  ONE_ROUTE_PER_FIT = values[["conditional_route_count"]],
  ONE_PUBLIC_CONVERGENCE_FORMAT_PATH = values[["public_convergence_formatter_count"]],
  ONE_AGGREGATED_WARNING_PATH = all(values[c(
    "public_convergence_warning_builder_count",
    "public_convergence_warning_call_count")]),
  ONE_CONVERGENCE_MEMORY_HELPER = values[["convergence_memory_helper_count"]],
  RAW_VERSION_ONE = values[["raw_schema_version"]],
  PHASE17U_NATIVE_UNCHANGED = values[["native_phase17u_hashes_unchanged"]],
  WRAPPERS_UNCHANGED = values[["wrapper_hashes_unchanged"]],
  NAMESPACE_UNCHANGED = values[["namespace_unchanged"]])
for (name in names(guards)) cat(name, "=", guards[[name]], "\n", sep = "")
stopifnot(all(guards))
