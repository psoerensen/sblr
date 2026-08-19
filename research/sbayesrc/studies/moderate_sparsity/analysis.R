spec <- source(file.path("studies", "06_annotation_models",
  "06d_moderate_sparsity", "spec.R"), local = TRUE)$value
evidence <- file.path("results", "local", "10_moderate_sparsity_sbayesrc")
invisible(list(spec = spec, evidence = evidence,
  evidence_available = dir.exists(evidence),
  status = "STUDY10-R2 — local, review-pending"))
