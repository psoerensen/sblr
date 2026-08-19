# Study 12 is scientifically preserved but paused pending redesign. This
# entry point validates traceability metadata and never fits or writes.

study12_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        any(readLines(file.path(path, "DESCRIPTION"), warn = FALSE) ==
          "Package: sblrbench")) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Cannot locate sblrbench root.",
      call. = FALSE)
    path <- parent
  }
}

study12_validation_only <- function(root = study12_root()) {
  old <- setwd(root); on.exit(setwd(old), add = TRUE)
  base <- file.path("studies", "12_logvar_annotation_evaluation")
  spec1 <- source(file.path(base, "spec.R"), local = TRUE)$value
  spec2 <- source(file.path(base, "spec_aspect2.R"), local = TRUE)$value
  documents <- file.path(base, c(
    "README.md", "ASPECTS.md", "design.md", "design_aspect2.md",
    "decision_aspect2_addendum.md", "report.qmd", "report_aspect2.qmd",
    "synthesis.qmd"))
  checks <- c(
    paused = identical(spec1$execution_status, "paused_pending_redesign") &&
      identical(spec2$execution_status, "paused_pending_redesign"),
    aspects = identical(spec1$aspect,
      "aspect1_canonical_informative_annotation") &&
      identical(spec2$aspect, "aspect2_sbayesr_lv_truth_v1"),
    historical_dataset = identical(spec1$historical_data$dataset_id,
      "human_independent") && identical(spec2$historical_data$dataset_id,
      "human_independent"),
    historical_qgdata = identical(spec1$historical_data$qgdata_sha,
      "6cca5819e711d326cfb2614d7e9d9f34942612cd") &&
      identical(spec2$historical_data$qgdata_sha,
      "6cca5819e711d326cfb2614d7e9d9f34942612cd"),
    local_only = identical(spec1$outputs$promote_reference, FALSE) &&
      identical(spec2$outputs$promote_reference, FALSE),
    recovery_checkpoint = identical(spec1$historical_source$git_checkpoint,
      "7a5cbf947af9b3441538fc64d7966c693e1433f8") &&
      identical(spec2$historical_source$git_checkpoint,
      "7a5cbf947af9b3441538fc64d7966c693e1433f8"),
    documents = all(file.exists(documents)))
  if (!all(checks)) stop("Study 12 traceability validation failed: ",
    paste(names(checks)[!checks], collapse = ", "), call. = FALSE)
  invisible(list(status = spec1$execution_status,
    aspects = c(spec1$aspect, spec2$aspect),
    historical_data = spec1$historical_data,
    historical_source = spec1$historical_source,
    reports = documents[grepl("[.]qmd$", documents)],
    evidence = file.path(root, "results", "local",
      "12_logvar_annotation_evaluation"),
    promoted = FALSE, writes = FALSE, checks = checks))
}

study12_execute <- function(gate = "validate", root = study12_root()) {
  if (gate %in% c("validate", "status"))
    return(study12_validation_only(root))
  stop("Study 12 is paused pending redesign; scientific gate `", gate,
    "` is deliberately unavailable in Phase 8B. Historical executable ",
    "sources remain recoverable from Git checkpoint ",
    "7a5cbf947af9b3441538fc64d7966c693e1433f8.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
inline <- args[startsWith(args, "--gate=")]
at <- match("--gate", args)
gate <- if (length(inline)) {
  sub("--gate=", "", inline[[1L]], fixed = TRUE)
} else if (!is.na(at) && at < length(args)) {
  args[[at + 1L]]
} else {
  "validate"
}
invisible(study12_execute(gate))
