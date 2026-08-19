spec <- source(file.path("studies", "06_annotation_models",
  "06c_high_information", "spec.R"), local = TRUE)$value
evidence <- file.path("results", "local", "08_high_information_sbayesrc")
invisible(list(spec = spec, evidence = evidence,
  evidence_available = dir.exists(evidence),
  status = "STUDY08-R2 — local, review-pending"))
