root <- normalizePath(if (file.exists("DESCRIPTION")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE)
bench <- data.frame(
  family = c("CSR BayesC", "scheduled CSR BayesC", "CSR BayesR", "CSR SBayesRC",
    "fixed-prior CSR", "group CSR", "learned-annotation CSR", "block-eigen",
    "legacy multivariate"),
  script = c("blr_phase4_csr_bayesc.R", "blr_phase10d_scheduled_csr.R",
    "blr_phase6_csr_bayesr.R", "blr_phase8_csr_sbayesrc.R",
    "blr_phase9c_csr_prior_bayesc.R", "blr_phase9e_csr_group_bayesc.R",
    "blr_phase9g_csr_learned_annotation_bayesc.R", NA, NA),
  baseline = c(rep("runtime/completed-fit RSS/I/O", 7),
    "tiny smoke tests only", "no canonical benchmark"),
  peak_rss = c(rep("opt-in", 7), "absent", "absent"),
  io = c(rep("CSR loaded once per fit; page-cache effects apply", 7),
    "CSR plus BED-derived block construction; page-cache effects apply",
    "in-memory dense/sparse summaries; CSR variant opens LD files per fit"),
  stringsAsFactors = FALSE)
bench$exists <- !is.na(bench$script) & file.exists(file.path(root, "tools", "benchmarks", bench$script))
cat("Phase 17A benchmark inventory; no cross-backend speed ranking.\n")
cat("Completed-fit RSS is not peak RSS; workloads are backend-specific.\n")
print(bench, row.names = FALSE)
