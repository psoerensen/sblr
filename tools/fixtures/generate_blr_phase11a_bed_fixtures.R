# Manual maintenance tool. Run one model per fresh R process, for example:
# Rscript tools/fixtures/generate_blr_phase11a_bed_fixtures.R bayesr
args <- commandArgs(trailingOnly = TRUE)
model <- if (length(args)) args[[1L]] else stop("supply bayesc, bayesr, or bayesrc")
stopifnot(model %in% c("bayesc", "bayesr", "bayesrc"))
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11a-bed-reference.R"))
value <- phase11a_capture(model, ncores = 1L,
  nchains = if (model == "bayesc") 2L else 1L, seed = 71L)
value$metadata <- list(reference_mode = "fresh R process", model = model,
  source_commit = "2a3cddb", rng_audit = if (model == "bayesc")
    "worker-thread-persistent distribution risk" else "chain-local distributions")
dir.create(file.path("tests", "testthat", "fixtures", "blr_phase11a"),
  recursive = TRUE, showWarnings = FALSE)
saveRDS(value, file.path("tests", "testthat", "fixtures", "blr_phase11a",
  paste0(model, ".rds")), version = 3)
