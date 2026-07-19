root <- normalizePath(if (file.exists("DESCRIPTION")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE)
setwd(root)
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("tests/testthat/fixtures/blr_phase17b_mt_default/blr-phase17b-mt-default-reference.R")
out <- "tests/testthat/fixtures/blr_phase17b_mt_default"
for (id in 1:3) {
  reference <- list(metadata = phase17b_mt_metadata(id),
    raw = phase17b_mt_capture(id, FALSE), fit = phase17b_mt_capture(id, TRUE))
  saveRDS(reference, file.path(out, sprintf("config-%d.rds", id)), version = 3)
}
message("Wrote three unchanged-public-path Phase 17B references.")
