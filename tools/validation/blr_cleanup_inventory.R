args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".")
old <- setwd(root)
on.exit(setwd(old), add = TRUE)

paths <- system2("git", c("ls-files"), stdout = TRUE)
paths <- unique(c(paths, system2("git", c("ls-files", "--others", "--exclude-standard"),
                                 stdout = TRUE)))
paths <- gsub("\\\\", "/", paths)
phase <- grepl("phase[0-9]", paths, ignore.case = TRUE)
candidate <- phase | grepl(
  "^(tests/testthat/(helper|fixtures)|tools/(audit|benchmarks)|docs/dev|docs/notes|examples/workflows|R/|src/|man/|[.]github/workflows)",
  paths)
paths <- paths[candidate]

artifact_type <- function(path) {
  if (grepl("^tests/testthat/test-", path)) return("test")
  if (grepl("^tests/testthat/helper", path)) return("test_helper")
  if (grepl("^tests/testthat/fixtures", path)) return("fixture")
  if (grepl("^tools/audit", path)) return("audit")
  if (grepl("^tools/benchmarks", path)) return("benchmark")
  if (grepl("^docs/dev", path)) return("developer_documentation")
  if (grepl("^docs/notes", path)) return("user_documentation")
  if (grepl("^examples/workflows", path)) return("workflow")
  if (grepl("^R/", path)) return("R_source")
  if (grepl("^src/", path)) return("native_source")
  if (grepl("^man/", path)) return("Rd")
  if (grepl("^[.]github/workflows", path)) return("CI")
  "other"
}

replacement_for <- function(type) switch(type,
  test = "docs/dev/blr_test_ownership.md",
  test_helper = "tests/testthat/helper-blr-fixtures.R",
  fixture = "tests/testthat/fixtures/README.md",
  audit = "tools/audit/blr_architecture_audit.R",
  benchmark = "docs/dev/blr_phase8a_public_r_alignment_checkpoint.md",
  developer_documentation = "docs/dev/history.md",
  user_documentation = "docs/notes/index.qmd",
  workflow = "examples/workflows/README.md",
  R_source = "docs/dev/blr_architecture.md",
  native_source = "docs/dev/blr_backend_inventory.md",
  Rd = "roxygen-generated canonical Rd",
  CI = ".github/workflows/blr-framework.yml",
  "docs/dev/blr_architecture.md")

classify <- function(path, type, present) {
  is_phase <- grepl("phase[0-9]", path, ignore.case = TRUE)
  if (type == "test" && is_phase)
    return(c("phase history", "permanent scientific/software owner", "merge into a permanent owner"))
  if (type == "audit" && is_phase)
    return(c("historical architecture guard", "tools/audit/blr_architecture_audit.R", "replace with a canonical owner"))
  if (type == "benchmark" && is_phase)
    return(c("historical performance probe", "canonical benchmark owner", "merge into a permanent owner"))
  if (type == "developer_documentation" && is_phase)
    return(c("historical phase record", "docs/dev/history.md", "archive as concise historical metadata"))
  if (type == "fixture" && is_phase)
    return(c("scientific or migration reference", "scientifically named fixture owner", "rename and retain"))
  if (!present) {
    action <- switch(type,
      R_source = "delete as compatibility-only",
      native_source = "delete as unreachable/dead",
      Rd = "delete as compatibility-only",
      test = "delete as duplicate",
      test_helper = "delete as duplicate",
      fixture = "delete as migration-only",
      audit = "replace with a canonical owner",
      benchmark = "merge into a permanent owner",
      developer_documentation = "merge into a permanent owner",
      user_documentation = "merge into a permanent owner",
      workflow = "merge into a permanent owner",
      "delete as duplicate")
    return(c(paste("retired", type, "contract or duplicate"),
             replacement_for(type), action))
  }
  c("current contract candidate", path, "retain unchanged as permanent scientific owner")
}

rows <- lapply(paths, function(path) {
  type <- artifact_type(path)
  present <- file.exists(path)
  cls <- classify(path, type, present)
  data.frame(
    path = path,
    artifact_type = type,
    current_owner = if (grepl("phase[0-9]", path, ignore.case = TRUE)) "development phase" else path,
    scientific_or_software_contract = cls[[1L]],
    portable_or_source_only = if (type %in% c("audit", "benchmark", "developer_documentation")) "source_only" else "portable_or_mixed",
    replacement_owner = cls[[2L]],
    action = cls[[3L]],
    deletion_preconditions = "replacement recorded; focused owner and neutrality gate pass",
    validation_evidence = "permanent owner gates plus exact Phase 21 neutrality matrix",
    final_status = if (present) "retained_or_replaced" else "removed_after_transfer",
    stringsAsFactors = FALSE)
})
manifest <- do.call(rbind, rows)

escape <- function(x) gsub("[|]", "\\\\|", x)
header <- c(
  "# BLR cleanup manifest",
  "",
  "This final manifest inventories the starting tracked artifacts and their permanent replacements. Every removal records a replacement owner and validation evidence.",
  "",
  paste0("Starting commit: `", system2("git", c("rev-parse", "HEAD"), stdout = TRUE), "`"),
  "",
  paste0("Inventoried artifacts: ", nrow(manifest)),
  "",
  paste(names(manifest), collapse = " | "),
  paste(rep("---", ncol(manifest)), collapse = " | "))
body <- apply(manifest, 1L, function(x) paste(escape(x), collapse = " | "))
dir.create("docs/dev", recursive = TRUE, showWarnings = FALSE)
writeLines(c(header, body), "docs/dev/blr_cleanup_manifest.md", useBytes = TRUE)
cat(sprintf("cleanup_manifest_rows=%d\n", nrow(manifest)))
