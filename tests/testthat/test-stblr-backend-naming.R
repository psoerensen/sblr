source_sblr_test_file <- function(path) {
  candidates <- c(path, file.path("..", "..", path))
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) stop("Could not find ", path, call. = FALSE)
  source(path)
}

if (!exists(".fit_stblr_bed_bayesr", mode = "function")) {
  source_sblr_test_file("R/sparse_ld_bed_helper.R")
}

test_that("backend-specific internal helper names are available", {
  expect_true(exists(".format_stblr_csr_bayesc_fit", mode = "function"))
  expect_true(exists(".format_stblr_csr_bayesr_fit", mode = "function"))
  expect_true(exists(".format_stblr_bed_bayesc_fit", mode = "function"))
  expect_true(exists(".format_stblr_bed_bayesr_fit", mode = "function"))

  expect_true(exists(".fit_stblr_csr_bayesc", mode = "function"))
  expect_true(exists(".fit_stblr_csr_bayesr", mode = "function"))
  expect_true(exists(".fit_stblr_bed_bayesc", mode = "function"))
  expect_true(exists(".fit_stblr_bed_bayesr", mode = "function"))
})

test_that("compatibility aliases for old internal BayesR names remain", {
  expect_true(exists(".stblr_csr_bayesr_experimental", mode = "function"))
  expect_true(exists(".stblr_bed_marker_bayesr_experimental", mode = "function"))
  expect_equal(names(formals(.stblr_bed_marker_bayesr_experimental)), "...")
})

test_that("supported public CSR and BED docs do not label BayesR experimental", {
  docs <- c(
    "man/stblr_csr.Rd",
    "man/stblr_csr_bayesr.Rd",
    "man/stblr_bed.Rd",
    "docs/dev/stblr_csr_bayesr_design.md",
    "docs/dev/archive/stblr_bayesr_backend_design.md",
    "docs/dev/stblr_bayesr_ldswap_design.md",
    "docs/dev/stblr_backend_naming.md"
  )
  docs <- docs[file.exists(docs)]
  text <- unlist(lapply(docs, readLines, warn = FALSE), use.names = FALSE)
  bayesr_lines <- grep("BayesR|bayesr|stblr_csr_bayesr|stblr_bed", text, value = TRUE)
  active_lines <- bayesr_lines[!grepl("compatibility alias|experimental\\(\\)", bayesr_lines)]
  expect_false(any(grepl("\\b[Ee]xperimental\\b", active_lines)))
})
