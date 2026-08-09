logvar_interface_fixture <- function() {
  prefix <- tempfile("logvar-interface-")
  m <- 6L
  writeLines(c(paste0("n_variants=", m), "nnz=0"),
             paste0(prefix, ".meta.txt"))
  writeBin(rep(as.raw(0), 8L * (m + 1L)),
           paste0(prefix, ".row_ptr.u64.bin"))
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))
  ids <- paste0("m", seq_len(m))
  list(
    prefix = prefix, ids = ids,
    stats = list(
      yy = stats::setNames(99, "T1"),
      ww = list(T1 = stats::setNames(rep(99, m), ids)),
      wy = list(T1 = stats::setNames(c(5, -4, 3, -2, 1, 0.5), ids)),
      n = 100L, m = m))
}

logvar_interface_cleanup <- function(fixture) {
  unlink(paste0(fixture$prefix, c(
    ".meta.txt", ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin")))
}

test_that("log-variance preprocessing follows the frozen mixed annotation contract", {
  ids <- paste0("m", 1:6)
  A <- data.frame(binary = c(0, 0, 0, 1, 1, 1),
                  continuous = c(2, 4, 5, 8, 10, 13), row.names = ids)
  out <- sblr:::.stblr_preprocess_logvar_annotations(A, ids)
  expect_equal(unname(out$X[, "binary"]), A$binary - mean(A$binary), tolerance = 0)
  expect_equal(mean(out$X[, "binary"]), 0, tolerance = 1e-15)
  expect_equal(mean(out$X[, "continuous"]), 0, tolerance = 1e-15)
  expect_equal(stats::sd(out$X[, "continuous"]), 1, tolerance = 1e-15)
  expect_identical(out$transform$type, c("binary", "continuous"))
  expect_identical(out$transform$transform,
                   c("center_only", "center_and_sd"))
  expect_equal(mean(log(exp(drop(out$X %*% c(0.4, -0.2))))), 0,
               tolerance = 1e-14)
})

test_that("log-variance preprocessing rejects unidentified designs", {
  ids <- paste0("m", 1:6)
  expect_error(sblr:::.stblr_preprocess_logvar_annotations(
    matrix(1, 6, 1, dimnames = list(ids, "intercept")), ids), "intercept")
  expect_error(sblr:::.stblr_preprocess_logvar_annotations(
    cbind(a = 1:6, b = 1:6), ids), "rank deficient")
  expect_error(sblr:::.stblr_preprocess_logvar_annotations(
    cbind(a = 1:6, b = rep(2, 6)), ids), "constant")
  bad <- matrix(1:6, 6, 1, dimnames = list(ids, "a")); bad[2] <- NA_real_
  expect_error(sblr:::.stblr_preprocess_logvar_annotations(bad, ids),
               "missing or non-finite")
  mismatched <- matrix(1:6, 6, 1,
                       dimnames = list(rev(ids), "a"))
  expect_error(sblr:::.stblr_preprocess_logvar_annotations(mismatched, ids),
               "exactly match")
})

test_that("public annotation interface dispatches BayesC-LV and BayesR-LV", {
  fixture <- logvar_interface_fixture()
  on.exit(logvar_interface_cleanup(fixture), add = TRUE)
  A <- data.frame(binary = c(0, 0, 0, 1, 1, 1),
                  continuous = c(2, 4, 5, 8, 10, 13),
                  row.names = fixture$ids)
  for (method in c("sbayesc", "sbayesr")) {
    fit <- stblr_csr_annot(
      fixture$stats, Glist = list(sparseLD = list(prefix = fixture$prefix)),
      annotations = A, annotation_model = "log_variance", method = method,
      nit = 20L, nburn = 5L, nchains = 2L, ncores = 1L, seed = 650L,
      convergence = "none")
    expect_identical(fit$model, paste0(method, "_logvar"))
    expect_identical(fit$annotation_model, "log_variance")
    expect_identical(fit$input$theta_prior_sd, 0.7)
    expect_identical(dim(fit$theta), c(2L, 1L))
    expect_identical(dim(fit$marker_prior_scale), c(6L, 1L))
    expect_true(all(is.finite(fit$marker_prior_scale)))
    expect_true(all(fit$marker_prior_scale > 0))
    expect_true(all(c("mean", "sd", "median", "lower", "upper", "Rhat",
                      "bulk_ESS", "tail_ESS", "MCSE") %in%
                    names(fit$theta_summary)))
    expect_true(all(c("theta_updates", "min_log_q", "max_log_q") %in%
                    names(fit$diagnostics$logvar)))
  }
})

test_that("public log-variance controls fail early and clearly", {
  fixture <- logvar_interface_fixture()
  on.exit(logvar_interface_cleanup(fixture), add = TRUE)
  A <- matrix(1:6, 6, 1, dimnames = list(fixture$ids, "a"))
  common <- list(
    stats = fixture$stats,
    Glist = list(sparseLD = list(prefix = fixture$prefix)),
    annotations = A, annotation_model = "log_variance", nit = 2L, nburn = 1L)
  expect_error(do.call(stblr_csr_annot, c(common, list(method = "sbayesrc"))),
               "sbayesc.*sbayesr")
  expect_error(do.call(stblr_csr_annot, c(common, list(
    method = "sbayesc", theta_prior_sd = 0))), "positive finite")
  expect_error(do.call(stblr_csr_annot, c(common, list(
    method = "sbayesc", maf_effect_s = 0))), "does not support")
})
