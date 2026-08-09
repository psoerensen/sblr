test_that("BayesC reusable-engine no-op policy preserves frozen trajectories", {
  root <- blr_repo_path()
  oldwd <- setwd(root)
  on.exit(setwd(oldwd), add = TRUE)
  source("tests/testthat/fixtures/st-bayesc-csr-reference.R", local = TRUE)
  for (name in names(st_bayesc_csr_reference_configurations)) {
    config <- st_bayesc_csr_reference_configurations[[name]]
    objects <- st_bayesc_csr_reference_objects(config)
    # The complete native raw object freezes every sampler trajectory/state and
    # schema field. Formatted fits contain data$operator$prefix, a fresh
    # tempfile value not normalized by the historical fixture helper, so their
    # serialized hashes are intentionally not used for byte identity.
    observed <- st_bayesc_csr_reference_md5(objects$raw)
    expect_identical(
      observed,
      unname(st_bayesc_csr_pre_engine_extraction_hashes[[name]][["raw"]]),
      info = name
    )
  }
})
