phase9d1_path <- function(...) {
  path <- file.path(...)
  if (file.exists(path)) path else file.path("..", "..", ...)
}

source(phase9d1_path(
  "tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"
))

phase9d1_active <- function(lines, text) {
  sum(grepl(text, lines, fixed = TRUE) & !grepl("^\\s*//", lines))
}

test_that("Phase 9D1 group execution is extracted once at the approved seam", {
  source_lines <- readLines(phase9d1_path("src", "st_cpg_omp_csr_group.cpp"), warn = FALSE)
  core_lines <- readLines(phase9d1_path("src", "blr_csr_group_bayesc_core_impl.h"), warn = FALSE)
  source_text <- paste(source_lines, collapse = "\n")
  core_text <- paste(core_lines, collapse = "\n")

  expect_equal(phase9d1_active(core_lines, "for (int isort = 0; isort < m; ++isort)"), 1L)
  expect_equal(phase9d1_active(source_lines, "for (int isort = 0; isort < m; ++isort)"), 0L)
  expect_equal(phase9d1_active(source_lines, "#include \"blr_csr_group_bayesc_core_impl.h\""), 1L)
  expect_match(core_text, "#ifndef SBLR_BLR_CSR_GROUP_BAYESC_CORE_IMPL_H", fixed = TRUE)
  expect_match(core_text, "#ifndef SBLR_CSR_GROUP_BAYESC_CORE_IMPL_TRANSLATION_UNIT", fixed = TRUE)
  expect_match(core_text, "const int g = group(iu);", fixed = TRUE)
  expect_match(core_text, "group_pi_t(gu)", fixed = TRUE)
  expect_match(core_text, "group_vb_multiplier_t(gu)", fixed = TRUE)
  expect_match(core_text, "sampleGroupVbMultipliers_ST_csr_group(", fixed = TRUE)
  expect_match(core_text, "normalize_group_vb", fixed = TRUE)
  expect_false(grepl("use_old|use_new|old_path|new_path|fallback", source_text,
                     ignore.case = TRUE))
  expect_match(source_text, "static Rcpp::List stblr_csr_group_bayesc_result_to_raw(", fixed = TRUE)
  expect_match(source_text, "for (int chain = 0; chain < nchains; ++chain)", fixed = TRUE)

  cpp <- list.files(phase9d1_path("src"), pattern = "\\.(cpp|h)$", full.names = TRUE)
  include_count <- sum(vapply(cpp, function(path) {
    any(grepl("#include \"blr_csr_group_bayesc_core_impl.h\"",
              readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1)))
  expect_equal(include_count, 1L)
})

test_that("Phase 9D1 permanent group references remain exact", {
  expect_length(phase9a_configs$group, 3L)
  for (name in names(phase9a_configs$group)) {
    reference <- readRDS(phase9d1_path(
      "tests", "testthat", "fixtures", "blr_phase9a_group", paste0(name, ".rds")
    ))
    config <- phase9a_configs$group[[name]]
    expect_identical(phase9a_normalize(phase9a_run("group", config, TRUE)), reference$raw,
                     info = paste(name, "raw"))
    expect_identical(phase9a_normalize(phase9a_run("group", config, FALSE)), reference$fit,
                     info = paste(name, "formatted"))
  }
})

test_that("Phase 9D1 group execution remains reproducible across core order", {
  comparable <- function(value) {
    value <- phase9a_normalize(value)
    value$input$ncores <- 0L
    value
  }
  config <- phase9a_configs$group$group_chains
  config$ncores <- 1L
  one <- comparable(phase9a_run("group", config, FALSE))
  expect_identical(comparable(phase9a_run("group", config, FALSE)), one)
  config$ncores <- 2L
  two <- comparable(phase9a_run("group", config, FALSE))
  expect_identical(two, one)
  expect_identical(comparable(phase9a_run("group", config, FALSE)), two)
  config$ncores <- 1L
  expect_identical(comparable(phase9a_run("group", config, FALSE)), one)
  invisible(phase9a_run("annotation", phase9a_configs$annotation$annot_fixed, FALSE))
  expect_identical(comparable(phase9a_run("group", config, FALSE)), one)
})

test_that("Phase 9D1 fixtures retain supported group policy coverage", {
  expect_identical(phase9a_configs$group$group_one$normalize, TRUE)
  expect_identical(phase9a_configs$group$group_chains$normalize, FALSE)
  expect_identical(phase9a_configs$group$group_explicit$seeds, c(41L, 42L))
  expect_identical(phase9a_configs$group$group_chains$keep, TRUE)
  expect_identical(phase9a_configs$group$group_explicit$keep, FALSE)
  fixture <- phase9a_inputs(1L)
  expect_identical(unname(fixture$group), c("coding", "background", "coding", "background"))
})

test_that("Phase 9D1 protects adjacent backends, wrappers, signatures, and schema", {
  protected <- c(
    "src/st_cpg_omp_csr.cpp" = "92dafc0266d5a0e72aea000224154cef",
    "src/blr_csr_bayesc_types.h" = "e5975c311c69fe536db57dd21f01334f",
    "src/blr_csr_bayesc_core_impl.h" = "f7c617cbfc172639c1f8aea1bd8b1876",
    "src/st_cpg_omp_csr_prior.cpp" = "cce51072da6ddc3c18d58ab3b1f3c6df",
    "src/blr_csr_prior_bayesc_types.h" = "f109314a249e7b084b57aa339d27fdf6",
    "src/blr_csr_prior_bayesc_core_impl.h" = "bba60cef54adbf30b05f29be20bff286",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/blr_csr_bayesr_types.h" = "bf1d4b73065207ca361c7abdab3cb253",
    "src/blr_csr_bayesr_core_impl.h" = "4dac6bef2df917613df8e1a827640303",
    "src/st_sbayesrc_omp_csr.cpp" = "8c1b03d8f5b93e6831ccbed856c77ead",
    "src/blr_csr_sbayesrc_types.h" = "103b2a1282b99069be963d4ff3da15c8",
    "src/blr_csr_sbayesrc_core_impl.h" = "d06ec2a530e8c914201ee22b6be65739",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/st_block_eigen.h" = "bec3bc1e41841ab77747e34dc9818574",
    "src/st_cpg_omp_individual.cpp" = "667a0445503ef9f6b23dbab1e0114b4d",
    "src/st_cpg_omp_individual_scheduled.cpp" = "0d726fe3faf5deec887381c1458ab6b6",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp" = "85a5e45e03c59ce62654496a2f076fe9",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp" = "5904c60b32165a7ae73bfc9d6c0f920c",
    "src/mtblr.cpp" = "419472a9d17afbf39edfcafb98bba459",
    "src/mt_cpg.cpp" = "49a2c308b127de69cfe3bdf9df2be227",
    "src/mt_cpg_arma.cpp" = "f911293210e4a29017f64a92769ec814",
    "src/mt_cpg_omp.cpp" = "4c2e24988bd3151674be3c8982a36118",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e",
    "src/RcppExports.cpp" = "b4859db0f6308fa7e38051ddcf32d245",
    "R/RcppExports.R" = "9d13ea00b326c7e0cd606194d13a8bca",
    "R/sparse_ld_bed_helper.R" = "142f674bbe43063280f2dda25fa30a64",
    "docs/dev/stblr_raw_schema.md" = "82ac9ba4b7d8edc6f3e16ee3a26d8466"
  )
  actual <- unname(tools::md5sum(vapply(names(protected), phase9d1_path, character(1))))
  expect_identical(actual, unname(protected))
})
