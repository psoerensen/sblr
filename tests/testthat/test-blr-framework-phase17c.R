phase17c_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
  winslash = "/", mustWork = TRUE)
phase17c_text <- function(path) paste(readLines(file.path(phase17c_root, path),
  warn = FALSE), collapse = "\n")
source(file.path(phase17c_root,
  "tests/testthat/fixtures/blr_phase17b_mt_default",
  "blr-phase17b-mt-default-reference.R"))
source(file.path(phase17c_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  "blr-phase17c-mt-default-corrected-reference.R"))

phase17c_reference <- function(id) readRDS(file.path(phase17c_root,
  "tests/testthat/fixtures/blr_phase17c_mt_default_corrected",
  sprintf("config-%d.rds", id)))

phase17c_public_source <- function() {
  core <- phase17c_text("src/mtblr.cpp")
  start <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr(",
    core, fixed = TRUE)[1]
  finish <- regexpr("std::vector<std::vector<std::vector<double>>>  mtblr_hybrid(",
    core, fixed = TRUE)[1]
  substr(core, start, finish - 1L)
}

test_that("corrected raw and formatted references retain exact structure", {
  for (id in 1:3) {
    ref <- phase17c_reference(id)
    expect_equal(phase17c_mt_capture(id, FALSE), ref$raw, tolerance = 1e-12,
      info = paste("raw config", id))
    expect_equal(phase17c_mt_capture(id, TRUE), ref$fit, tolerance = 1e-12,
      info = paste("formatted config", id))
    expect_identical(ref$metadata$reference_mode,
      "structure_exact_numeric_tolerance")
    expect_identical(ref$metadata$numeric_tolerance, 1e-12)
    expect_true(isTRUE(ref$metadata$structure_exact))
    expect_identical(length(ref$raw), 20L)
  }
})

test_that("fixed update controls preserve supplied B, E, and pi", {
  config <- phase17c_mt_config(2L)
  fit <- phase17c_mt_capture(2L, TRUE)
  expect_identical(unname(fit$vb), unname(config$vb))
  expect_identical(unname(fit$ve), unname(config$ve))
  expect_identical(unname(fit$pi), c(.8, rep(.2 / 3, 3)))
  expect_true(all(fit$covb == 0))
  expect_true(all(fit$cove == 0))
  expect_true(all(fit$pim == 0))
  expect_true(all(is.finite(unlist(fit))))
})

test_that("updated-B final trajectory is preserved from Phase 17B", {
  stable_raw <- c(3:10, 14:17)
  stable_fit <- c("wy", "r", "b", "d", "o", "vbs", "vgs", "ves",
    "vb", "vg", "ve", "pi")
  for (id in c(1L, 3L)) {
    legacy <- readRDS(file.path(phase17c_root,
      "tests/testthat/fixtures/blr_phase17b_mt_default",
      sprintf("config-%d.rds", id)))
    current_raw <- phase17c_mt_capture(id, FALSE)
    current_fit <- phase17c_mt_capture(id, TRUE)
    for (field in stable_raw)
      expect_equal(current_raw[[field]], legacy$raw[[field]], tolerance = 1e-12)
    for (field in stable_fit)
      expect_equal(current_fit[[field]], legacy$fit[[field]], tolerance = 1e-12)
  }
})

test_that("retained-count metadata follows the corrected iteration policy", {
  for (id in 1:3) {
    config <- phase17c_mt_config(id)
    expected <- phase17c_mt_retained_counts(config)
    metadata <- phase17c_reference(id)$metadata
    for (field in names(expected))
      expect_identical(metadata[[field]], expected[[field]])
    expect_identical(expected$marker_retained_count,
      as.integer(ceiling(config$nit / config$nthin)))
  }
})

test_that("updated probability and covariance summaries use their own counts", {
  for (id in c(1L, 3L)) {
    fit <- phase17c_mt_capture(id, TRUE)
    config <- phase17c_mt_config(id)
    retained_rows <- (config$nburn + 1L):(config$nburn + config$nit)
    expect_equal(sum(unname(fit$pim)), 1, tolerance = 1e-12)
    expect_equal(diag(fit$covb), colMeans(fit$vbs[retained_rows, , drop = FALSE]),
      tolerance = 1e-12)
    expect_equal(diag(fit$covg), colMeans(fit$vgs[retained_rows, , drop = FALSE]),
      tolerance = 1e-12)
    expect_equal(diag(fit$cove), colMeans(fit$ves[retained_rows, , drop = FALSE]),
      tolerance = 1e-12)
  }
})

test_that("first post-burn iteration is eligible for marker summaries", {
  config <- phase17c_mt_config(2L)
  config$nit <- 1L
  config$nburn <- 0L
  config$nthin <- 2L
  raw <- phase17c_mt_native_raw(config)
  expect_equal(raw[[1]], Map(function(effect, state) effect * (state > 0),
    raw[[5]], raw[[6]]), tolerance = 1e-12)
  expect_equal(raw[[2]], lapply(raw[[6]], function(state) as.numeric(state > 0)),
    tolerance = 1e-12)
})

test_that("same-process, thread environment, and intervening fits are stable", {
  a <- phase17c_mt_capture(1L, TRUE)
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
  invisible(phase17c_mt_capture(2L, TRUE))
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else
    Sys.setenv(OMP_NUM_THREADS = old), add = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "1")
  one <- phase17c_mt_capture(1L, TRUE)
  Sys.setenv(OMP_NUM_THREADS = "2")
  expect_equal(phase17c_mt_capture(1L, TRUE), one, tolerance = 1e-12)
  source(file.path(phase17c_root,
    "tests/testthat/fixtures/blr-phase5a-bayesr-reference.R"),
    local = environment())
  invisible(phase5a_bayesr_run(phase5a_bayesr_configs$one_chain, raw = FALSE))
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
  source(file.path(phase17c_root,
    "tests/testthat/fixtures/blr-phase13a-bed-bayesr-reference.R"),
    local = environment())
  invisible(phase13a_capture(ncores = 1L, nchains = 1L, seed = 71L))
  expect_equal(phase17c_mt_capture(1L, TRUE), a, tolerance = 1e-12)
})

test_that("fresh-process corrected reference is stable", {
  skip_if_not(identical(Sys.getenv("SBLR_RUN_PHASE17C_FRESH"), "true"))
  skip_if_not_installed("callr")
  observed <- callr::r(function(root) {
    setwd(root); pkgload::load_all(".", compile = FALSE, quiet = TRUE)
    source(paste0("tests/testthat/fixtures/blr_phase17b_mt_default/",
      "blr-phase17b-mt-default-reference.R"))
    source(paste0("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/",
      "blr-phase17c-mt-default-corrected-reference.R"))
    phase17c_mt_capture(1L, TRUE)
  }, list(root = phase17c_root))
  expect_equal(observed, phase17c_reference(1L)$fit, tolerance = 1e-12)
})

test_that("corrected scientific identities and legacy schema are stable", {
  for (id in 1:3) {
    fit <- phase17c_mt_capture(id, TRUE)
    expect_identical(dim(fit$bm), c(4L, 2L))
    expect_identical(rownames(fit$bm), paste0("M", 1:4))
    expect_identical(colnames(fit$bm), c("TraitA", "TraitB"))
    expect_true(all(is.finite(unlist(fit))))
    expect_true(all(fit$dm >= 0 & fit$dm <= 1))
    expect_true(all(fit$d %in% 0:1))
    for (field in c("covb", "covg", "cove", "vb", "vg", "ve")) {
      x <- fit[[field]]
      expect_equal(x, t(x), tolerance = 1e-12)
      expect_true(min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
    }
  }
  expect_source_count('fit <- .Call("_sblr_mtblr"',
    phase17c_text("R/interface_mtblr.R"), 1L)
})

test_that("public source freezes guards, conditions, and distinct counts", {
  public <- phase17c_public_source()
  expect_source_count("if (updateB) {", public, 1L)
  expect_source_count("sampleBset(nt, m, nub, B", public, 1L)
  expect_source_count("sampleB_latent(nt, m, nub, B", public, 1L)
  expect_source_count("if(updateB && method==4)", public, 1L)
  expect_source_count("it >= nburn", public, 10L)
  expect_source_count("(it - nburn) % nthin", public, 2L)
  for (name in c("marker_retained_count", "covb_retained_count",
      "covg_retained_count", "cove_retained_count", "pi_retained_count"))
    expect_true(source_match_count(name, public) >= 3L)
  expect_source_forbidden(public, c("it>nburn", "it > nburn",
    "bm[t][i]/nsamples", "pis[i]/nit"))
  expect_source_count("std::mt19937 gen(seed);", public, 1L)
  expect_source_forbidden(public, c("omp_get_thread_num", "static std::mt19937",
    "thread_local"))
})

test_that("alternative and canonical backend files remain protected", {
  protected <- c(
    "src/mt_cpg.cpp" = "49a2c308b127de69cfe3bdf9df2be227",
    "src/mt_cpg_arma.cpp" = "f911293210e4a29017f64a92769ec814",
    "src/mt_cpg_omp.cpp" = "4c2e24988bd3151674be3c8982a36118",
    "src/mt_cpg_omp_csr.cpp" = "aec85896b5c30db3014efaeb5e3c3a96",
    "src/st_block_eigen.cpp" = "49f0a62c9fe235967a264b0f8de144a7",
    "src/st_block_eigen.h" = "bec3bc1e41841ab77747e34dc9818574",
    "src/st_cpg_omp_csr.cpp" = "92dafc0266d5a0e72aea000224154cef",
    "src/st_cpg_omp_csr_bayesr.cpp" = "0a005f9d5a19037285fd4869fdc4dcf0",
    "src/st_sbayesrc_omp_csr.cpp" = "8c1b03d8f5b93e6831ccbed856c77ead",
    "R/RcppExports.R" = "9d13ea00b326c7e0cd606194d13a8bca",
    "src/RcppExports.cpp" = "b4859db0f6308fa7e38051ddcf32d245",
    "NAMESPACE" = "f5b6ee37a3972aa436357bdc8f602f4e")
  expect_identical(unname(tools::md5sum(file.path(phase17c_root,
    names(protected)))), unname(protected))
  omp <- phase17c_text("src/mt_cpg_omp.cpp")
  expect_match(omp, "seed + 100000 * it + omp_get_thread_num()", fixed = TRUE)

  core <- phase17c_text("src/mtblr.cpp")
  region_md5 <- function(start, finish = NULL) {
    first <- regexpr(start, core, fixed = TRUE)[1]
    last <- if (is.null(finish)) nchar(core) + 1L else
      regexpr(finish, core, fixed = TRUE)[1]
    path <- tempfile(fileext = ".cpp")
    on.exit(unlink(path), add = TRUE)
    writeChar(substr(core, first, last - 1L), path, eos = NULL, useBytes = TRUE)
    unname(tools::md5sum(path))
  }
  hybrid <- "std::vector<std::vector<std::vector<double>>>  mtblr_hybrid("
  eigen <- "std::vector<std::vector<std::vector<double>>>  mtblr_eigen("
  expect_identical(region_md5(hybrid, eigen),
    "03ed109628b874223db109d2ec654827")
  expect_identical(region_md5(eigen), "4e5c38ede3345de10a684ab38470bf7b")
})
