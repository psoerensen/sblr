test_that("individual BayesRC helpers are internal", {
  expect_true(exists(".stblr_bed_bayesrc_native", envir = asNamespace("sblr")))
  expect_false(".stblr_bed_bayesrc_native" %in% getNamespaceExports("sblr"))
  expect_true(exists(".bayesr_pi_to_probit_stick_intercepts", envir = asNamespace("sblr")))
  expect_false(".bayesr_pi_to_probit_stick_intercepts" %in% getNamespaceExports("sblr"))
})

test_that("individual BayesRC uses the shared BED utility header", {
  source_path <- blr_repo_path(
    "src", "stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc.cpp")
  source <- paste(readLines(source_path, warn = FALSE), collapse = "\n")
  expect_match(source, '#include "st_bed_bayesr_common.h"', fixed = TRUE)
  expect_false(grepl(
    'include "stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp"',
    source,
    fixed = TRUE
  ))
})

test_that("intercept conversion reproduces BayesR component probabilities", {
  target <- c(0.95, 0.03, 0.015, 0.005)
  intercept <- sblr:::.bayesr_pi_to_probit_stick_intercepts(target)
  stick <- stats::pnorm(intercept[1L, ])
  got <- numeric(length(target))
  remaining <- 1
  for (k in seq_along(stick)) {
    got[k] <- remaining * (1 - stick[k])
    remaining <- remaining * stick[k]
  }
  got[length(got)] <- remaining
  expect_equal(got, target, tolerance = 1e-12)
})

test_that("probit stick conversion handles multiple priors and rejects invalid input", {
  priors <- list(
    c(0.95, 0.03, 0.015, 0.005),
    c(0.70, 0.20, 0.08, 0.02),
    rep(0.25, 4)
  )
  for (target in priors) {
    stick <- stats::pnorm(sblr:::.bayesr_pi_to_probit_stick_intercepts(target))
    got <- numeric(length(target))
    remaining <- 1
    for (k in seq_along(stick)) {
      got[k] <- remaining * (1 - stick[k])
      remaining <- remaining * stick[k]
    }
    got[length(got)] <- remaining
    expect_true(all(is.finite(got) & got >= 0 & got <= 1))
    expect_equal(sum(got), 1, tolerance = 1e-12)
    expect_equal(got, target / sum(target), tolerance = 1e-12)
  }
  for (bad in list(1, c(-1, 2), c(0, 1), c(NA, 1), c(Inf, 1))) {
    expect_error(sblr:::.bayesr_pi_to_probit_stick_intercepts(bad), "pi")
  }
  for (bad_floor in c(0, -1, 0.5, 1, NA, Inf)) {
    expect_error(
      sblr:::.bayesr_pi_to_probit_stick_intercepts(c(0.5, 0.5), bad_floor),
      "pi_floor"
    )
  }
})

test_that("individual BayesRC native symbol is registered when compiled", {
  ok <- tryCatch({
    getNativeSymbolInfo(
      "_sblr_stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc",
      PACKAGE = "sblr"
    )
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native individual BayesRC symbol is not loaded")
  expect_true(is.function(sblr:::.stblr_bed_bayesrc_native))
})

make_individual_bayesrc_fixture <- function(
    dosage = rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1)),
    y = matrix(c(-1, 0, 1, -0.5, 0.5, 1.5), ncol = 1L)) {
  bed <- tempfile(fileext = ".bed")
  code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(j) {
    z <- unname(code[as.character(dosage[j, ])])
    z <- c(z, rep(0L, (-length(z)) %% 4L))
    vapply(seq(1L, length(z), 4L), function(i) {
      sum(z[i:(i + 3L)] * c(1L, 4L, 16L, 64L))
    }, integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), bed)
  list(
    bed = bed,
    dosage = dosage,
    n = ncol(dosage),
    m = nrow(dosage),
    af = rowMeans(dosage) / 2,
    y = y,
    gamma = c(0, 0.01, 0.1, 1),
    pi = c(0.95, 0.03, 0.015, 0.005)
  )
}

skip_if_no_individual_bayesrc <- function() {
  ok <- tryCatch({
    getNativeSymbolInfo(
      "_sblr_stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc",
      PACKAGE = "sblr"
    )
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "native individual BayesRC symbol is not loaded")
}

make_individual_bayesrc_args <- function(
    fixture = make_individual_bayesrc_fixture(), y = fixture$y,
    updateAlpha = FALSE, nchains = 1L, keep_chains = FALSE,
    ncores = 1L, seed = 17L, nit = 6L, nburn = 2L, nthin = 1L,
    A = matrix(1, fixture$m, 1L)) {
  nt <- ncol(y)
  list(
    bed_files = fixture$bed, n = fixture$n, cls = list(seq_len(fixture$m)), y = y,
    b_init = replicate(nt, numeric(fixture$m), simplify = FALSE),
    sets = rep(1L, fixture$m), rows = NULL, af = list(fixture$af), scale = TRUE,
    B = diag(0.1, nt), E = diag(1, nt),
    ssb_prior = replicate(nt, rep(0.05, nt), simplify = FALSE),
    sse_prior = replicate(nt, rep(0.5, nt), simplify = FALSE),
    A = A, gamma = fixture$gamma,
    annot_alpha_init = rbind(
      sblr:::.bayesr_pi_to_probit_stick_intercepts(fixture$pi),
      matrix(0, max(0, ncol(A) - 1L), length(fixture$gamma) - 1L)
    ),
    annot_sigma_sq_alpha_init = rep(1, length(fixture$gamma) - 1L),
    intercept_prior_resolved = sblr:::.sbayesrc_resolve_intercept_prior(
      fixture$pi)$native,
    updateAlpha = updateAlpha, annot_alpha_update_every = 1L,
    updateB = TRUE, updateE = TRUE, adjE = 0.9,
    nit = nit, nburn = nburn, nthin = nthin, rebuild_every = 1L,
    return_wy = TRUE, return_r = TRUE, read_block_size = 2L,
    nchains = nchains, keep_chains = keep_chains,
    ncores = ncores, seed = seed
  )
}

run_individual_bayesrc <- function(...) {
  do.call(sblr:::.stblr_bed_bayesrc_native, make_individual_bayesrc_args(...))
}

expect_individual_bayesrc_raw <- function(raw, m, K, nt, ntrace) {
  expect_s3_class(raw, "stblr_raw_v1")
  expect_identical(raw$schema$class, "stblr_raw")
  expect_identical(as.integer(raw$schema$version), 1L)
  expect_identical(raw$meta$model, "bayesrc")
  expect_identical(raw$meta$backend, "bed_bayesrc")
  expect_identical(raw$meta$data_level, "individual")
  expect_identical(raw$meta$prior_type, "annotation_component")
  expect_true(isTRUE(raw$meta$annotations))
  expect_false(isTRUE(raw$meta$scheduled))
  expect_equal(
    setdiff(
      c("schema", "meta", "marker", "trace", "variance", "pi",
        "diagnostics", "chains", "prior", "group", "annotation",
        "component", "selection"),
      names(raw)
    ),
    character()
  )
  for (field in c("bm", "dm", "b", "state")) expect_equal(dim(raw$marker[[field]]), c(m, nt))
  for (field in c("vbs", "vgs", "ves", "vle", "vld", "pis")) {
    expect_equal(dim(raw$trace[[field]]), c(ntrace, nt))
    expect_true(all(is.finite(raw$trace[[field]])))
  }
  expect_equal(raw$trace$vld, raw$trace$vgs - raw$trace$vle, tolerance = 1e-12)
  expect_true(all(raw$trace$vbs > 0 & raw$trace$ves > 0))
  expect_identical(raw$component$names[1L], "gamma_0.00")
  expect_equal(length(raw$component$names), K)
  expect_equal(as.numeric(raw$component$mixture_var), c(0, 0.01, 0.1, 1))
  for (t in seq_len(nt)) {
    cp <- raw$component$prob[[t]]
    expect_equal(dim(cp), c(m, K))
    expect_true(all(is.finite(cp) & cp >= 0 & cp <= 1))
    expect_equal(rowSums(cp), rep(1, m), tolerance = 1e-12)
    expect_equal(raw$marker$dm[, t], 1 - cp[, 1L], tolerance = 1e-12)
    expect_equal(
      raw$component$dm_component_mean[, t],
      drop(cp %*% seq.int(0, K - 1L)),
      tolerance = 1e-12
    )
    expect_equal(raw$component$ncomp[t, ], colSums(cp), tolerance = 1e-12)
    expect_equal(sum(raw$component$ncomp[t, ]), m, tolerance = 1e-12)
    prior <- raw$annotation$marker_prior_final[[t]]
    expect_equal(dim(prior), c(m, K))
    expect_true(all(is.finite(prior) & prior >= 0 & prior <= 1))
    expect_equal(rowSums(prior), rep(1, m), tolerance = 1e-12)
    expect_equal(raw$pi$final[t, ], colMeans(prior), tolerance = 1e-12)
  }
  expect_true(isTRUE(raw$diagnostics$full_sweeps))
  expect_false(isTRUE(raw$diagnostics$adaptive_skipping))
}

test_that("individual BayesRC native backend runs and reduces to fixed-pi BayesR", {
  bayesrc_ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  bayesr_ok <- tryCatch({
    getNativeSymbolInfo("_sblr_stblr_cpg_omp_bed_marker_scheduled_chains_bayesr", PACKAGE = "sblr")
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(bayesrc_ok && bayesr_ok, "native BED BayesRC/BayesR symbols are not loaded")

  x <- make_individual_bayesrc_fixture()
  intercept <- sblr:::.bayesr_pi_to_probit_stick_intercepts(x$pi)
  common <- list(
    bed_files = x$bed, n = 6L, cls = list(1:2), y = x$y,
    b_init = list(c(0, 0)), sets = c(1L, 1L), rows = NULL,
    af = list(c(0.5, 0.5)), scale = TRUE,
    B = matrix(0.1, 1L, 1L), E = matrix(1, 1L, 1L),
    ssb_prior = list(0.05), sse_prior = list(0.5),
    nub = 4, nue = 4, updateB = TRUE, updateE = TRUE,
    adjE = 0.9, nit = 5L, nburn = 2L, nthin = 1L,
    rebuild_every = 1L, return_wy = TRUE, return_r = TRUE,
    read_block_size = 2L, nchains = 2L, ncores = 1L, seed = 17L
  )
  raw_rc <- do.call(sblr:::.stblr_bed_bayesrc_native, c(common, list(
    A = matrix(1, 2L, 1L), gamma = x$gamma,
    annot_alpha_init = intercept,
    annot_sigma_sq_alpha_init = rep(1, length(x$gamma) - 1L),
    intercept_prior_resolved = sblr:::.sbayesrc_resolve_intercept_prior(
      x$pi)$native,
    updateAlpha = FALSE, annot_alpha_update_every = 10L
  )))
  expect_s3_class(raw_rc, "stblr_raw_v1")
  expect_identical(raw_rc$meta$model, "bayesrc")
  expect_identical(raw_rc$meta$backend, "bed_bayesrc")
  expect_true(isTRUE(raw_rc$diagnostics$full_sweeps))
  expect_false(isTRUE(raw_rc$diagnostics$adaptive_skipping))
  expect_null(raw_rc$chains)
  expect_identical(raw_rc$meta$nchains, 2L)
  cp <- raw_rc$component$prob[[1L]]
  expect_identical(colnames(cp), NULL)
  expect_true(all(is.finite(cp) & cp >= 0 & cp <= 1))
  expect_equal(rowSums(cp), rep(1, nrow(cp)), tolerance = 1e-12)
  expect_equal(raw_rc$marker$dm[, 1L], 1 - cp[, 1L], tolerance = 1e-12)
  expect_identical(raw_rc$component$names[1L], "gamma_0.00")
  expect_equal(raw_rc$pi$final[1L, ], x$pi, tolerance = 1e-12)
  expect_equal(unname(raw_rc$annotation$alpha_mean[[1L]]), unname(intercept), tolerance = 1e-12)
  expect_equal(dim(raw_rc$annotation$sigmaSqAlpha_mean), c(3L, 1L))

  raw_r <- do.call(sblr:::stblr_cpg_omp_bed_marker_scheduled_chains_bayesr, c(common, list(
    pi = x$pi, c = x$gamma, alpha = rep(1, length(x$pi)),
    updatePi = FALSE, full_sweep_every = 1L, null_skip_base = 1L,
    null_skip_max = 1L, candidate_threshold = 0,
    candidate_lifetime = 0L, skip_nulls_burnin_only = FALSE,
    progress_every = 0L
  )))
  expect_equal(raw_rc$marker$bm, raw_r$marker$bm, tolerance = 1e-12)
  expect_equal(raw_rc$marker$dm, raw_r$marker$dm, tolerance = 1e-12)
  expect_equal(raw_rc$component$prob[[1L]], raw_r$component$prob[[1L]], tolerance = 1e-12)
  for (field in c("b", "state")) {
    expect_equal(raw_rc$marker[[field]], raw_r$marker[[field]], tolerance = 1e-12)
  }
  for (field in c("vbs", "vgs", "ves", "vle", "vld")) {
    expect_equal(raw_rc$trace[[field]], raw_r$trace[[field]], tolerance = 1e-12)
  }
  expect_equal(raw_rc$diagnostics$log_cpo, raw_r$diagnostics$log_cpo, tolerance = 1e-12)
  expect_equal(
    raw_rc$diagnostics$mean_log_cpo,
    raw_r$diagnostics$mean_log_cpo,
    tolerance = 1e-12
  )
})

test_that("BED BayesRC component, annotation, variance, CPO, and schema outputs are coherent", {
  skip_if_no_individual_bayesrc()
  raw <- run_individual_bayesrc(updateAlpha = TRUE, nit = 6L, nburn = 2L)
  expect_individual_bayesrc_raw(raw, m = 2L, K = 4L, nt = 1L, ntrace = 8L)
  expect_length(raw$annotation$alpha_mean, 1L)
  expect_length(raw$annotation$alpha_final, 1L)
  expect_equal(dim(raw$annotation$alpha_mean[[1L]]), c(1L, 3L))
  expect_equal(dim(raw$annotation$alpha_final[[1L]]), c(1L, 3L))
  expect_true(all(is.finite(raw$annotation$alpha_mean[[1L]])))
  expect_true(all(is.finite(raw$annotation$alpha_final[[1L]])))
  expect_equal(dim(raw$annotation$sigmaSqAlpha_mean), c(3L, 1L))
  expect_equal(dim(raw$annotation$sigmaSqAlpha_final), c(3L, 1L))
  expect_true(all(is.finite(raw$annotation$sigmaSqAlpha_mean) &
                  raw$annotation$sigmaSqAlpha_mean > 0))
  expect_true(all(is.finite(raw$annotation$sigmaSqAlpha_final) &
                  raw$annotation$sigmaSqAlpha_final > 0))
  expect_identical(raw$diagnostics$annotation_updates_per_chain, 8L)
  expect_true(all(is.finite(raw$trace$pis) & raw$trace$pis >= 0 & raw$trace$pis <= 1))
  expect_length(raw$diagnostics$log_cpo, 1L)
  expect_length(raw$diagnostics$mean_log_cpo, 1L)
  expect_length(raw$diagnostics$nsamples, 1L)
  expect_true(all(is.finite(raw$diagnostics$log_cpo)))
  expect_equal(
    raw$diagnostics$mean_log_cpo,
    raw$diagnostics$log_cpo / raw$diagnostics$n_used,
    tolerance = 1e-12
  )
  expect_equal(raw$diagnostics$nsamples, 6)

  fixed <- run_individual_bayesrc(updateAlpha = FALSE, nit = 6L, nburn = 2L, nthin = 2L)
  prior_nonnull <- mean(1 - fixed$annotation$marker_prior_final[[1L]][, 1L])
  expect_equal(fixed$trace$pis[, 1L], rep(prior_nonnull, 8L), tolerance = 1e-12)
  expect_equal(fixed$diagnostics$nsamples, 3)
  expect_false(isTRUE(all.equal(
    mean(fixed$trace$pis[, 1L]),
    sum(fixed$component$ncomp[1L, -1L]) / 2
  )))
})

test_that("BED BayesRC learns directional annotation enrichment", {
  skip_if_no_individual_bayesrc()
  set.seed(7041)
  n <- 80L
  m <- 16L
  dosage <- do.call(rbind, lapply(seq_len(m), function(j) {
    stats::rbinom(n, 2L, stats::runif(1L, 0.2, 0.45))
  }))
  expect_true(all(apply(dosage, 1L, stats::var) > 0))
  expect_identical(anyDuplicated(lapply(seq_len(m), function(j) dosage[j, ])), 0L)
  annotated <- seq_len(6L)
  causal <- seq_len(4L)
  x <- t(scale(t(dosage), center = 2 * rowMeans(dosage) / 2,
               scale = sqrt(2 * (rowMeans(dosage) / 2) * (1 - rowMeans(dosage) / 2))))
  set.seed(7042)
  y <- drop(crossprod(x[causal, , drop = FALSE], c(1.8, -1.5, 1.3, -1.1))) +
    stats::rnorm(n, sd = 0.45)
  fixture <- make_individual_bayesrc_fixture(dosage, matrix(y, ncol = 1L))
  annotation <- cbind(intercept = 1, enriched = as.integer(seq_len(m) %in% annotated))
  fit_enrichment <- function(ncores) run_individual_bayesrc(
    fixture = fixture, y = fixture$y, A = annotation,
    updateAlpha = TRUE, nchains = 2L, ncores = ncores,
    nit = 160L, nburn = 40L, seed = 331L
  )
  raw <- fit_enrichment(1L)
  repeat_raw <- fit_enrichment(1L)
  two_core_raw <- fit_enrichment(2L)
  enrichment_output <- function(x) list(
    marker = x$marker[c("bm", "dm")],
    component_prob = x$component$prob,
    alpha_mean = x$annotation$alpha_mean,
    marker_prior_final = x$annotation$marker_prior_final
  )
  expect_identical(enrichment_output(raw), enrichment_output(repeat_raw))
  expect_identical(enrichment_output(raw), enrichment_output(two_core_raw))
  prior <- raw$annotation$marker_prior_final[[1L]]
  nonnull <- 1 - prior[, 1L]
  expect_gt(mean(nonnull[annotated]), mean(nonnull[-annotated]))
  expect_gt(raw$annotation$alpha_mean[[1L]][2L, 1L], 0)
  expect_true(all(is.finite(raw$annotation$alpha_mean[[1L]])))
})

test_that("BED BayesRC supports separate traits and retained chains", {
  skip_if_no_individual_bayesrc()
  fixture <- make_individual_bayesrc_fixture()
  y2 <- cbind(fixture$y[, 1L], -0.5 * fixture$y[, 1L] + c(0, 0.1, 0, -0.1, 0, 0.2))
  raw <- run_individual_bayesrc(
    fixture = fixture, y = y2, updateAlpha = TRUE,
    nchains = 2L, keep_chains = TRUE, nit = 4L, nburn = 2L
  )
  expect_individual_bayesrc_raw(raw, m = 2L, K = 4L, nt = 2L, ntrace = 6L)
  expect_identical(raw$meta$nt, 2L)
  expect_identical(raw$meta$nchains, 2L)
  expect_length(raw$component$prob, 2L)
  expect_length(raw$annotation$alpha_mean, 2L)
  expect_length(raw$diagnostics$log_cpo, 2L)
  expect_length(raw$chains, 2L)
  for (t in 1:2) {
    expect_length(raw$chains[[t]], 2L)
    for (chain in raw$chains[[t]]) {
      expect_equal(
        setdiff(c("bm", "dm", "b", "state", "comp_prob", "alpha",
                  "sigmaSqAlpha", "pis"), names(chain)),
        character()
      )
      expect_equal(dim(chain$bm), c(1L, 2L))
      expect_equal(dim(chain$dm), c(1L, 2L))
      expect_equal(dim(chain$comp_prob), c(2L, 4L))
      expect_equal(dim(chain$alpha), c(1L, 3L))
    }
    expect_equal(
      raw$marker$bm[, t],
      drop(Reduce(`+`, lapply(raw$chains[[t]], `[[`, "bm")) / 2),
      tolerance = 1e-12
    )
    expect_equal(
      raw$marker$dm[, t],
      drop(Reduce(`+`, lapply(raw$chains[[t]], `[[`, "dm")) / 2),
      tolerance = 1e-12
    )
    expect_equal(
      raw$component$prob[[t]],
      Reduce(`+`, lapply(raw$chains[[t]], `[[`, "comp_prob")) / 2,
      tolerance = 1e-12
    )
    expect_equal(
      raw$annotation$alpha_mean[[t]],
      Reduce(`+`, lapply(raw$chains[[t]], `[[`, "alpha")) / 2,
      tolerance = 1e-12
    )
    expect_equal(
      raw$annotation$sigmaSqAlpha_mean[, t],
      drop(Reduce(`+`, lapply(raw$chains[[t]], `[[`, "sigmaSqAlpha")) / 2),
      tolerance = 1e-12
    )
  }
  expect_false(isTRUE(all.equal(raw$marker$bm[, 1L], raw$marker$bm[, 2L])))
  expect_null(run_individual_bayesrc(nchains = 2L, keep_chains = FALSE)$chains)
})

test_that("BED BayesRC is reproducible by seed and OpenMP thread count", {
  skip_if_no_individual_bayesrc()
  one <- run_individual_bayesrc(updateAlpha = TRUE, nchains = 2L, ncores = 1L, seed = 91L)
  repeat_one <- run_individual_bayesrc(updateAlpha = TRUE, nchains = 2L, ncores = 1L, seed = 91L)
  two <- run_individual_bayesrc(updateAlpha = TRUE, nchains = 2L, ncores = 2L, seed = 91L)
  different <- run_individual_bayesrc(updateAlpha = TRUE, nchains = 2L, ncores = 1L, seed = 92L)
  stochastic_output <- function(x) list(
    marker = x$marker[c("bm", "dm", "b", "state")],
    trace = x$trace[c("vbs", "vgs", "ves", "vle", "vld", "pis")],
    component_prob = x$component$prob,
    alpha_mean = x$annotation$alpha_mean,
    alpha_final = x$annotation$alpha_final,
    sigma_mean = x$annotation$sigmaSqAlpha_mean,
    sigma_final = x$annotation$sigmaSqAlpha_final,
    log_cpo = x$diagnostics$log_cpo,
    mean_log_cpo = x$diagnostics$mean_log_cpo
  )
  expect_identical(stochastic_output(one), stochastic_output(repeat_one))
  expect_identical(stochastic_output(one), stochastic_output(two))

  two_first <- run_individual_bayesrc(updateAlpha = TRUE, nchains = 2L, ncores = 2L, seed = 91L)
  one_after_two <- run_individual_bayesrc(updateAlpha = TRUE, nchains = 2L, ncores = 1L, seed = 91L)
  expect_identical(stochastic_output(one), stochastic_output(two_first))
  expect_identical(stochastic_output(one), stochastic_output(one_after_two))

  invisible(run_individual_bayesrc(updateAlpha = FALSE, nchains = 1L, ncores = 1L, seed = 812L))
  after_unrelated <- run_individual_bayesrc(updateAlpha = TRUE, nchains = 2L, ncores = 1L, seed = 91L)
  expect_identical(stochastic_output(one), stochastic_output(after_unrelated))
  expect_false(isTRUE(all.equal(one$marker$bm, different$marker$bm)))
})

test_that("BED BayesRC follows optional wy and residual conventions", {
  skip_if_no_individual_bayesrc()
  raw <- run_individual_bayesrc()
  expect_equal(dim(raw$marker$wy), c(2L, 1L))
  expect_equal(dim(raw$marker$r), c(2L, 1L))
  expect_true(all(is.finite(raw$marker$wy)))
  expect_true(all(is.finite(raw$marker$r)))
  args <- make_individual_bayesrc_args()
  args$return_wy <- args$return_r <- FALSE
  absent <- do.call(sblr:::.stblr_bed_bayesrc_native, args)
  expect_null(absent$marker$wy)
  expect_null(absent$marker$r)
})

test_that("BED BayesRC validates native inputs clearly", {
  skip_if_no_individual_bayesrc()
  base <- make_individual_bayesrc_args()
  check_bad <- function(name, value, pattern) {
    args <- base
    args[[name]] <- value
    expect_error(do.call(sblr:::.stblr_bed_bayesrc_native, args), pattern)
  }
  check_bad("A", matrix(1, 1L, 1L), "A")
  check_bad("A", matrix(numeric(), 2L, 0L), "A")
  check_bad("A", matrix(c(1, NA), 2L), "A")
  check_bad("A", matrix(c(1, Inf), 2L), "A")
  check_bad("annot_alpha_init", matrix(0, 2L, 3L), "annot_alpha_init")
  check_bad("annot_alpha_init", matrix(0, 1L, 2L), "annot_alpha_init")
  check_bad("annot_alpha_init", matrix(c(NA, 0, 0), 1L), "annotation initial")
  check_bad("annot_sigma_sq_alpha_init", c(1, 1), "annot_sigma")
  check_bad("annot_sigma_sq_alpha_init", c(1, 0, 1), "variances positive")
  check_bad("gamma", 0, "gamma")
  check_bad("gamma", c(0.1, 1), "exact zero")
  check_bad("gamma", c(0, 0), "positive")
  check_bad("gamma", c(0, -1), "positive")
  check_bad("gamma", c(0, Inf), "positive")
  for (bad in c(0, -1, 1, NA, Inf)) check_bad("pi_floor", bad, "pi_floor")
  check_bad("annot_alpha_update_every", 0L, "annot_alpha_update_every")
  check_bad("nit", 0L, "MCMC")
  check_bad("nburn", -1L, "MCMC")
  check_bad("nthin", 0L, "MCMC")
  check_bad("nchains", 0L, "MCMC")
  check_bad("ncores", 0L, "MCMC")
  check_bad("y", matrix(1, 5L, 1L), "phenotype")
  check_bad("b_init", list(c(0, 0), c(0, 0)), "sets or b_init")
  check_bad("b_init", list(0), "each b_init")
  check_bad("sets", 1L, "sets or b_init")
  check_bad("B", diag(2), "B and E")
  check_bad("E", diag(2), "B and E")
  check_bad("bed_files", tempfile(fileext = ".bed"), "Could not open BED")
  check_bad("cls", list(integer()), "no markers selected")
})
