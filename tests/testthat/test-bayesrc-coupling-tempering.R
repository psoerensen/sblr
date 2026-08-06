coupling_tempering_control <- function(prior, enabled = TRUE, swap_every = 5L,
                                       update_variance = TRUE) {
  if (nrow(prior) == 3L) {
    prior <- rbind(
      prior,
      update_sigmaSqAlpha = rep(as.integer(update_variance), ncol(prior))
    )
  }
  if (nrow(prior) == 4L) {
    prior <- rbind(
      prior,
      allocation_updates_per_cycle = rep(1, ncol(prior)),
      annotation_updates_per_cycle = rep(1, ncol(prior))
    )
  }
  rbind(
    prior,
    coupling_tempering = rep(as.integer(enabled), ncol(prior)),
    coupling_swap_every = rep(if (enabled) swap_every else 0L, ncol(prior))
  )
}

make_coupling_bed_fixture <- function(
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
    bed = bed, dosage = dosage, n = ncol(dosage), m = nrow(dosage),
    af = rowMeans(dosage) / 2, y = y,
    gamma = c(0, 0.01, 0.1, 1), pi = c(0.95, 0.03, 0.015, 0.005)
  )
}

make_coupling_bed_args <- function(
    fixture = make_coupling_bed_fixture(), y = fixture$y,
    updateAlpha = FALSE, nchains = 1L, keep_chains = FALSE,
    ncores = 1L, seed = 17L, nit = 6L, nburn = 2L, nthin = 1L,
    A = matrix(1, fixture$m, 1L)) {
  nt <- ncol(y)
  list(
    bed_files = fixture$bed, n = fixture$n, cls = list(seq_len(fixture$m)),
    y = y, b_init = replicate(nt, numeric(fixture$m), simplify = FALSE),
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
      fixture$pi
    )$native,
    updateAlpha = updateAlpha, annot_alpha_update_every = 1L,
    updateB = TRUE, updateE = TRUE, adjE = 0.9,
    nit = nit, nburn = nburn, nthin = nthin, rebuild_every = 1L,
    return_wy = TRUE, return_r = TRUE, read_block_size = 2L,
    nchains = nchains, keep_chains = keep_chains,
    ncores = ncores, seed = seed
  )
}

test_that("coupling endpoints reproduce baseline and learned BayesRC priors", {
  annotation <- cbind(intercept = 1, binary = c(0, 1, 0, 1))
  mixture <- c(0.80, 0.15, 0.05)
  gamma <- c(0, 0.1, 1)
  baseline <- drop(sblr:::.bayesr_pi_to_probit_stick_intercepts(mixture))
  alpha <- rbind(baseline + c(0.2, -0.1), c(0.7, -0.4))

  at_one <- sblr:::.st_bayesrc_tempered_probabilities(
    annotation, alpha, baseline, 1, 1e-12
  )
  ordinary <- sblr:::sbayesrc_marker_pi(annotation, alpha, gamma)
  expect_equal(unname(at_one), unname(ordinary), tolerance = 1e-14)

  at_zero <- sblr:::.st_bayesrc_tempered_probabilities(
    annotation, alpha, baseline, 0, 1e-12
  )
  expect_equal(at_zero, matrix(mixture, nrow(annotation), 3L, byrow = TRUE),
               tolerance = 1e-12)
  for (coupling in c(0, 0.5, 1)) {
    probability <- sblr:::.st_bayesrc_tempered_probabilities(
      annotation, alpha, baseline, coupling, 1e-12
    )
    expect_true(all(is.finite(probability) & probability > 0))
    expect_equal(rowSums(probability), rep(1, nrow(annotation)),
                 tolerance = 1e-14)
  }
})

test_that("coupling swap ratio satisfies the exact allocation-prior identity", {
  annotation <- cbind(intercept = 1, binary = c(0, 1, 0, 1, 1))
  baseline <- qnorm(0.2)
  alpha_lower <- matrix(c(baseline, 0.3), 2L, 1L)
  alpha_upper <- matrix(c(baseline + 0.5, -0.7), 2L, 1L)
  component_lower <- c(0, 0, 1, 0, 1)
  component_upper <- c(1, 0, 0, 1, 1)
  log_ratio <- sblr:::.st_bayesrc_coupling_swap_log_ratio(
    annotation, alpha_lower, component_lower, alpha_upper, component_upper,
    baseline, 0.5, 1, 1e-12
  )
  probability <- function(alpha, coupling, component) {
    p <- sblr:::.st_bayesrc_tempered_probabilities(
      annotation, alpha, baseline, coupling, 1e-12
    )
    sum(log(p[cbind(seq_len(nrow(p)), component + 1L)]))
  }
  reference <- probability(alpha_upper, 0.5, component_upper) +
    probability(alpha_lower, 1, component_lower) -
    probability(alpha_lower, 0.5, component_lower) -
    probability(alpha_upper, 1, component_upper)
  expect_equal(log_ratio, reference, tolerance = 1e-13)
  reverse <- sblr:::.st_bayesrc_coupling_swap_log_ratio(
    annotation, alpha_upper, component_upper, alpha_lower, component_lower,
    baseline, 0.5, 1, 1e-12
  )
  expect_equal(reverse, -log_ratio, tolerance = 1e-13)
  expect_equal(
    min(1, exp(log_ratio)), exp(log_ratio) * min(1, exp(reverse)),
    tolerance = 1e-13
  )
})

test_that("disabled coupling is exactly RNG and output neutral", {
  base <- make_coupling_bed_args(
    updateAlpha = TRUE, nit = 5L, nburn = 1L, keep_chains = TRUE, seed = 9761L
  )
  disabled <- base
  disabled$intercept_prior_resolved <- coupling_tempering_control(
    base$intercept_prior_resolved, enabled = FALSE
  )
  ordinary <- do.call(sblr:::.stblr_bed_bayesrc_native, base)
  explicit_disabled <- do.call(sblr:::.stblr_bed_bayesrc_native, disabled)
  expect_identical(explicit_disabled, ordinary)
})

test_that("tiny BED coupling ensemble is finite, reproducible, and target-only", {
  fixture <- make_coupling_bed_fixture(
    dosage = rbind(
      c(0, 1, 2, 0, 1, 2, 0, 1),
      c(2, 1, 0, 1, 2, 0, 1, 2),
      c(0, 0, 1, 1, 2, 2, 1, 0),
      c(2, 1, 2, 0, 1, 0, 2, 1)
    ),
    y = matrix(c(-1, -0.3, 0.8, -0.5, 1.1, 0.4, -0.8, 0.3), ncol = 1L)
  )
  annotation <- cbind(intercept = 1, binary = c(0, 0, 1, 1))
  args <- make_coupling_bed_args(
    fixture = fixture, A = annotation, updateAlpha = TRUE,
    nit = 12L, nburn = 2L, keep_chains = TRUE, seed = 9762L
  )
  args$convergence_markers <- 0:(fixture$m - 1L)
  args$convergence_annotations <- TRUE
  args$convergence_b <- TRUE
  args$convergence_d <- TRUE
  args$convergence_component <- TRUE
  args$intercept_prior_resolved <- coupling_tempering_control(
    args$intercept_prior_resolved, enabled = TRUE, swap_every = 2L
  )
  first <- do.call(sblr:::.stblr_bed_bayesrc_native, args)
  second <- do.call(sblr:::.stblr_bed_bayesrc_native, args)
  first$chains[[1L]][[1L]]$coupling_tempering$transition_seconds <- NULL
  first$chains[[1L]][[1L]]$coupling_tempering$swap_seconds <- NULL
  second$chains[[1L]][[1L]]$coupling_tempering$transition_seconds <- NULL
  second$chains[[1L]][[1L]]$coupling_tempering$swap_seconds <- NULL
  expect_identical(second, first)
  chain <- first$chains[[1L]][[1L]]
  tempering <- chain$coupling_tempering
  expect_equal(dim(tempering$replica_identity), c(12L, 3L))
  expect_equal(dim(tempering$active_count), c(12L, 3L))
  expect_equal(dim(tempering$expected_active_count), c(12L, 3L))
  expect_equal(ncol(tempering$swap), 9L)
  expect_true(all(is.finite(tempering$swap)))
  expect_true(all(tempering$swap[, 4L] >= 0 & tempering$swap[, 4L] <= 1))
  expect_true(all(is.finite(chain$convergence_trace$alpha)))
  expect_true(all(is.finite(chain$convergence_trace$sigmaSqAlpha)))
  expect_true(all(chain$convergence_trace$component >= 0L &
                    chain$convergence_trace$component < length(fixture$gamma)))
  expect_equal(dim(chain$convergence_trace$component), c(12L, fixture$m))
})
