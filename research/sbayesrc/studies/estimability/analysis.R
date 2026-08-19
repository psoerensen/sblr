# Study 06A is closed. This entry point validates accepted compact evidence;
# it does not reopen fitting or regenerate scientific results.
spec <- source(file.path("studies", "06_annotation_models", "06a_estimability",
  "spec.R"), local = TRUE)$value
source(file.path("studies", "06_annotation_models", "06a_estimability",
  "estimability.R"), local = TRUE)

capsule <- file.path("results", "reference", "06_annotation_models",
  "current-stop")
inventory <- utils::read.csv(file.path(capsule, "checksums.csv"),
  stringsAsFactors = FALSE)
paths <- file.path(capsule, inventory$file)
if (any(!file.exists(paths)) ||
    !identical(unname(as.numeric(file.info(paths)$size)),
      as.numeric(inventory$size_bytes)) ||
    !identical(unname(tools::md5sum(paths)), inventory$md5))
  stop("The historical Study 06A stop capsule failed its raw-MD5 inventory.",
    call. = FALSE)

accepted <- c(
  file.path("results", "reference", "06_annotation_models",
    "final_decision.json"),
  file.path("results", "reference", "06_annotation_models",
    "final_cross_implementation_comparison.csv"),
  file.path("results", "reference", "06_annotation_models",
    "final_hierarchy_of_evidence.csv"))
if (any(!file.exists(accepted)))
  stop("Accepted Study 06A decision evidence is incomplete.", call. = FALSE)

invisible(list(spec = spec, capsule = capsule, accepted = accepted,
  reproduction = study06_estimability_validate("."),
  status = "CLOSED — EST-R2"))
