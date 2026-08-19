# Single entry point for SBayesRV Research Gate 1.
# Run from the sblr repository root with:
#   Rscript research/sbayesrv/analysis.R

research_dir <- file.path("research", "sbayesrv")
if (!file.exists(file.path(research_dir, "prototype.R"))) {
  stop("Run analysis.R from the sblr repository root.", call. = FALSE)
}

source(file.path(research_dir, "prototype.R"), local = FALSE)
source(file.path(research_dir, "qualification.R"), local = FALSE)

qualification <- run_sbayesrv_qualification(stop_on_failure = TRUE)
print(qualification, row.names = FALSE)

output_dir <- file.path(research_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  qualification,
  file.path(output_dir, "gate1_deterministic_qualification.csv"),
  row.names = FALSE
)

gradient_rows <- qualification[grepl("gradient$", qualification$gate), ]
summary_lines <- c(
  "# SBayesRV Research Gate 1 deterministic summary",
  "",
  paste0("Checks: ", nrow(qualification)),
  paste0("Passed: ", sum(qualification$pass)),
  paste0("Failed: ", sum(!qualification$pass)),
  paste0(
    "Maximum conditional-gradient discrepancy: ",
    format(max(qualification$value[
      qualification$gate == "conditional_gradient"]), digits = 16)
  ),
  paste0(
    "Maximum collapsed-gradient discrepancy: ",
    format(max(qualification$value[
      qualification$gate == "collapsed_gradient"]), digits = 16)
  ),
  paste0(
    "Maximum gradient tolerance: ",
    format(max(gradient_rows$tolerance), digits = 16)
  ),
  "",
  "This is deterministic research qualification, not predictive or benchmark evidence."
)
writeLines(summary_lines, file.path(output_dir, "gate1_summary.md"))

invisible(qualification)
