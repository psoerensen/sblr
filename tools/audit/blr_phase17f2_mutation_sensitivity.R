root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)))
  commandArgs(trailingOnly = TRUE)[[1L]] else ".", winslash = "/", mustWork = TRUE)
read_text <- function(path) paste(readLines(file.path(root, path), warn = FALSE),
  collapse = "\n")
must_change <- function(label, original, mutated, detector) {
  if (identical(original, mutated)) stop(label, ": mutation did not apply")
  if (!isTRUE(detector(mutated))) stop(label, ": owning contract did not detect mutation")
  cat(sprintf("PASS\t%s\n", label))
}

core <- read_text("src/blr_mt_default_core_impl.h")
must_change("remove updateB guard", core,
  sub("if (execution.updateB) {", "if (true) {", core, fixed = TRUE),
  function(x) !grepl("if (execution.updateB) {", x, fixed = TRUE))
must_change("strict post-burn", core,
  sub("it >= execution.nburn", "it > execution.nburn", core, fixed = TRUE),
  function(x) grepl("it > execution.nburn", x, fixed = TRUE))
must_change("absolute thinning", core,
  sub("(it - execution.nburn) % execution.nthin", "it % execution.nthin",
    core, fixed = TRUE),
  function(x) grepl("it % execution.nthin", x, fixed = TRUE))

finalizer <- read_text("src/blr_mt_default_finalize_impl.h")
must_change("wrong pi denominator", finalizer,
  sub("result.pi_retained_count", "result.pi_retained_count - 1.0",
    finalizer, fixed = TRUE),
  function(x) grepl("pi_retained_count - 1.0", x, fixed = TRUE))

seed <- read_text("src/blr_bed_family_types.h")
must_change("logical-chain seed mapping", seed,
  sub("9176", "9177", seed, fixed = TRUE),
  function(x) !grepl("9176", x, fixed = TRUE))

adapter <- read_text("src/mtblr.cpp")
must_change("legacy output position swap", adapter,
  sub("result[17][t][i] = pi_mean[i];", "result[16][t][i] = pi_mean[i];",
    adapter, fixed = TRUE),
  function(x) !grepl("result[17][t][i] = pi_mean[i];", x, fixed = TRUE))

formatter <- read_text("R/interface_mtblr.R")
orientation_token <- "fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)"
must_change("MT trait orientation", formatter,
  sub(orientation_token,
    "fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = FALSE)",
    formatter, fixed = TRUE),
  function(x) !grepl(orientation_token, x, fixed = TRUE))

types <- read_text("src/blr_mt_default_types.h")
must_change("Rcpp type in binding-neutral header", types,
  sub("namespace sblr", "Rcpp::NumericVector forbidden_binding;\nnamespace sblr",
    types, fixed = TRUE),
  function(x) grepl("Rcpp::", x, fixed = TRUE))

cat("MUTATION_SENSITIVITY=8/8\n")
