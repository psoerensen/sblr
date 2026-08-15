bayesr_engine_reference_fixture <- function() {
  prefix <- tempfile("bayesr-engine-reference-")
  m <- 4L
  writeLines(c(paste0("n_variants=", m), "nnz=0"),
             paste0(prefix, ".meta.txt"))
  writeBin(rep(as.raw(0), 8L * (m + 1L)),
           paste0(prefix, ".row_ptr.u64.bin"))
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))
  ids <- paste0("m", seq_len(m))
  list(
    prefix = prefix,
    stats = list(
      yy = stats::setNames(79, "T1"),
      ww = list(T1 = stats::setNames(rep(79, m), ids)),
      wy = list(T1 = stats::setNames(c(4, -3, 2, 1), ids)),
      n = 80L,
      m = m
    )
  )
}

test_that("BayesR reusable-engine no-op policy preserves frozen trajectory", {
  fixture <- bayesr_engine_reference_fixture()
  on.exit(unlink(paste0(fixture$prefix, c(
    ".meta.txt", ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin"))), add = TRUE)

  fit <- blr_with_legacy_execution(function() {
    stblr_csr(
      stats = fixture$stats,
      Glist = list(sparseLD = list(prefix = fixture$prefix)),
      method = "sbayesr", mixture_var = c(0, 0.1, 1),
      nit = 12L, nburn = 4L, nchains = 2L, ncores = 1L, seed = 902L,
      keep_chains = TRUE, convergence = "none"
    )
  })
  trajectory <- fit[c(
    "bm", "dm", "b", "d", "b_final", "d_final", "vbs", "vgs",
    "ves", "vle", "vld", "pis", "pi", "pim",
    "component_probabilities", "dm_component_mean", "chains", "ld_swap",
    "ld_swap_chains"
  )]
  serialized <- tempfile(fileext = ".rds")
  on.exit(unlink(serialized), add = TRUE)
  saveRDS(trajectory, serialized, version = 3)
  # pkgload's development namespace changes serialization metadata relative
  # to the installed-library fresh-process hash. The installed hash was frozen
  # before extraction and reproduced exactly afterward; accept only these two
  # known encodings of the same complete trajectory object.
  expect_true(unname(tools::md5sum(serialized)) %in% c(
    "a210f6ab3144f08ba5796ac01da30329",
    "5cb387613e7f2285c09a030f150ae8c8"
  ))
})
