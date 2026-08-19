spec <- source(file.path("studies", "06_annotation_models",
  "06e_implementation_concordance", "spec.R"), local = TRUE)$value
evidence <- file.path("results", "local", "11_gctb_sbayesrc_concordance")
invisible(list(spec = spec, evidence = evidence,
  evidence_available = dir.exists(evidence),
  status = "A4B4 does not resolve mixing; production GCTB Phase B deferred"))
