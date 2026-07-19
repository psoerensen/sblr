compile_source <- !identical(Sys.getenv("SBLR_FIXTURE_USE_EXISTING_BUILD"), "true")
if (!"sblr" %in% loadedNamespaces())
  pkgload::load_all(".", compile = compile_source, quiet = TRUE)
source(paste0("tests/testthat/fixtures/blr_phase17b_mt_default/",
  "blr-phase17b-mt-default-reference.R"))
source(paste0("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/",
  "blr-phase17c-mt-default-corrected-reference.R"))

out <- "tests/testthat/fixtures/blr_phase17c_mt_default_corrected"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (id in 1:3) {
  reference <- list(
    metadata = phase17c_mt_metadata(id),
    raw = phase17c_mt_capture(id, FALSE),
    fit = phase17c_mt_capture(id, TRUE)
  )
  path <- file.path(out, sprintf("config-%d.rds", id))
  saveRDS(reference, path, version = 3)
  message("Wrote ", path)
}
