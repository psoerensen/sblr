# Manual maintenance tool; never sourced by ordinary tests.
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11b-bed-bayesc-reference.R"))
out <- file.path("tests", "testthat", "fixtures", "blr_phase11b_bed_bayesc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
specs <- list(single_1x1 = c("single", 1L, 1L),
  multichain_2x1 = c("multichain", 2L, 1L),
  multichain_2x2 = c("multichain", 2L, 2L))
for (nm in names(specs)) {
  z <- specs[[nm]]
  value <- phase11b_capture(z[[1]], as.integer(z[[3]]), as.integer(z[[2]]), 71L)
  saveRDS(list(metadata = phase11b_metadata(z[[1]], as.integer(z[[2]]),
    as.integer(z[[3]])), raw = value$raw, fit = value$fit),
    file.path(out, paste0(nm, ".rds")), version = 3)
}
