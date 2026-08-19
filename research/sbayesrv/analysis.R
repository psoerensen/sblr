# Deterministic SBayesRV Gate 1 qualification only.
# Run from the sblr repository root:
#   Rscript research/sbayesrv/analysis.R

research_dir <- file.path("research", "sbayesrv")
if (!file.exists(file.path(research_dir, "prototype.R")) ||
    !file.exists("DESCRIPTION")) {
  stop("Run analysis.R from the sblr repository root.", call. = FALSE)
}

source(file.path(research_dir, "prototype.R"), local = FALSE)
source(file.path(research_dir, "qualification.R"), local = FALSE)
qualification <- run_sbayesrv_qualification(stop_on_failure = TRUE)
print(qualification, row.names = FALSE)
cat(sprintf("\nSBayesRV deterministic qualification: %d/%d passed.\n",
            sum(qualification$pass), nrow(qualification)))
invisible(qualification)
