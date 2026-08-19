spec <- source(file.path("studies", "06_annotation_models",
  "06b_joint_vs_em", "spec.R"), local = TRUE)$value
evidence <- file.path("results", "local", "07_joint_em_sbayesrc")
invisible(list(spec = spec, evidence = evidence,
  evidence_available = dir.exists(evidence),
  status = "STUDY07-R5 — EM demonstration failed; review required"))
