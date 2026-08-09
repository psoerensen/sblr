logvar_bayesc_fixture_args <- function(
  update_ld_swap = FALSE,
  keep_chains = TRUE,
  nchains = 2L
) {
  root <- blr_repo_path()
  oldwd <- setwd(root)
  on.exit(setwd(oldwd), add = TRUE)
  source("tests/testthat/fixtures/st-bayesc-csr-reference.R", local = TRUE)
  prefix <- st_bayesc_csr_reference_write_csr()
  inputs <- st_bayesc_csr_reference_inputs(1L)
  stats <- inputs$stats
  vy <- as.numeric(stats$yy) / (stats$n - 1)
  list(
    wy = stats$wy,
    ww = stats$ww,
    yy = stats$yy,
    b_init = list(rep(0, stats$m)),
    d_init = list(rep(0, stats$m)),
    use_d_init = FALSE,
    r_init = stats$wy,
    use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE,
    ld_prefix = prefix,
    B = diag((vy * 0.3) / (stats$m * 0.5), 1),
    E = diag(vy * 0.7, 1),
    ssb_prior = list(((4 - 2) / 4) * (vy * 0.3) / (stats$m * 0.5)),
    sse_prior = list(((4 - 2) / 4) * (vy * 0.7)),
    pi = c(0.5, 0.5),
    nub = 4,
    nue = 4,
    updateB = TRUE,
    updateE = TRUE,
    updatePi = TRUE,
    adjE = 0,
    n = stats$n,
    nit = 30L,
    nburn = 10L,
    nthin = 1L,
    pi_prior_a = 1,
    pi_prior_b = 1,
    ncores = 1L,
    seed = 907L,
    nchains = nchains,
    keep_chains = keep_chains,
    chain_seeds = if (nchains == 2L) c(1907L, 2907L) else integer(),
    updateLDswap = update_ld_swap,
    ld_swap_prob = if (update_ld_swap) 1 else 0,
    ld_swap_r2 = 0.1,
    ld_swap_max_friends = 10L,
    ld_swap_moves = if (update_ld_swap) 2L else 0L,
    convergence_markers = 0:2,
    convergence_b = TRUE,
    convergence_d = TRUE
  )
}

logvar_bayesc_core_trajectory <- function(raw) {
  diagnostics <- raw$diagnostics
  diagnostics$seconds_mean[] <- 0
  diagnostics$seconds_max[] <- 0
  diagnostics$logvar <- NULL
  list(
    marker = raw$marker,
    trace = raw$trace,
    variance = raw$variance,
    pi = raw$pi,
    diagnostics = diagnostics,
    chains = raw$chains
  )
}

test_that("BayesC-LV theta zero is trajectory-identical to ordinary BayesC", {
  args <- logvar_bayesc_fixture_args(update_ld_swap = FALSE)
  ordinary <- do.call(sblr:::stblr_cpg_omp_csr, args)
  lv <- do.call(sblr:::stblr_cpg_omp_csr_logvar_bayesc, c(args, list(
    annotation = matrix(c(-1, 0, 1), ncol = 1),
    theta_init = matrix(0, 1, 1),
    theta_prior_sd = 0.7,
    updateTheta = FALSE
  )))

  expect_identical(
    logvar_bayesc_core_trajectory(lv),
    logvar_bayesc_core_trajectory(ordinary)
  )
  expect_identical(drop(lv$annotation$theta), 0)
  expect_identical(drop(lv$annotation$marker_prior_scale), rep(1, 3))
  expect_identical(lv$diagnostics$logvar$theta_updates, 0)
})

test_that("fixed-theta BayesC-LV equals the fixed marker-scale route", {
  args <- logvar_bayesc_fixture_args(update_ld_swap = TRUE, nchains = 1L)
  X <- matrix(c(-1, 0, 1), ncol = 1)
  theta <- 0.35
  q <- exp(drop(X) * theta)
  ordinary <- do.call(sblr:::stblr_cpg_omp_csr, c(args, list(
    maf_effect_s_prior_scale = q
  )))
  lv <- do.call(sblr:::stblr_cpg_omp_csr_logvar_bayesc, c(args, list(
    annotation = X,
    theta_init = matrix(theta, 1, 1),
    theta_prior_sd = 0.7,
    updateTheta = FALSE
  )))

  expect_identical(
    logvar_bayesc_core_trajectory(lv),
    logvar_bayesc_core_trajectory(ordinary)
  )
  expect_equal(drop(lv$annotation$theta), theta, tolerance = 1e-15)
  expect_equal(drop(lv$annotation$marker_prior_scale), q, tolerance = 1e-15)
  expect_gt(lv$diagnostics$ld_swap[1, 1], 0)
})

test_that("learned BayesC-LV returns finite theta, q, and ESS diagnostics", {
  args <- logvar_bayesc_fixture_args(update_ld_swap = FALSE)
  lv <- do.call(sblr:::stblr_cpg_omp_csr_logvar_bayesc, c(args, list(
    annotation = matrix(c(-1, 0, 1), ncol = 1),
    theta_init = matrix(0, 1, 1),
    theta_prior_sd = 0.7,
    updateTheta = TRUE
  )))

  expect_true(all(is.finite(lv$annotation$theta)))
  expect_true(all(is.finite(lv$annotation$theta_trace)))
  expect_true(all(is.finite(lv$annotation$marker_prior_scale)))
  expect_true(all(lv$annotation$marker_prior_scale > 0))
  expect_equal(lv$diagnostics$logvar$theta_updates, 80)
  expect_gte(
    lv$diagnostics$logvar$mean_likelihood_evaluations_per_update, 0
  )
  expect_true(is.finite(lv$diagnostics$logvar$min_log_q))
  expect_true(is.finite(lv$diagnostics$logvar$max_log_q))
})

logvar_bayesc_oracle_ess <- function(theta, loglik, prior_sd) {
  current <- loglik(theta)
  nu <- stats::rnorm(length(theta), 0, prior_sd)
  threshold <- current + log(stats::runif(1))
  angle <- stats::runif(1, 0, 2 * pi)
  lower <- angle - 2 * pi
  upper <- angle
  evaluations <- contractions <- 0L
  for (step in seq_len(10000L)) {
    proposal <- theta * cos(angle) + nu * sin(angle)
    value <- loglik(proposal)
    evaluations <- evaluations + 1L
    if (is.finite(value) && value > threshold) {
      return(list(theta = proposal, evaluations = evaluations,
                  contractions = contractions))
    }
    contractions <- contractions + 1L
    if (angle < 0) lower <- angle else upper <- angle
    angle <- stats::runif(1, lower, upper)
  }
  stop("test oracle ESS exceeded its bracket guard")
}

logvar_bayesc_oracle_chain <- function(
  xy, yy, wi, X, iterations, burnin, seed, vb_init, ve_init,
  vb_prior, ve_prior
) {
  set.seed(seed)
  m <- length(xy)
  theta <- 0
  q <- rep(1, m)
  b <- numeric(m)
  state <- integer(m)
  pi_active <- 0.1
  theta_draw <- numeric(iterations - burnin)
  pi_draw <- vb_draw <- ve_draw <- theta_draw
  b_sum <- pip_sum <- q_sum <- numeric(m)
  vb <- vb_init
  ve <- ve_init
  saved <- 0L
  for (iteration in seq_len(iterations)) {
    for (marker in sample.int(m)) {
      vbi <- vb * q[marker]
      denominator <- ve + wi[marker] * vbi
      log_weight <- c(
        log(1 - pi_active),
        log(pi_active) + 0.5 * log(ve / denominator) +
          0.5 * xy[marker]^2 * vbi / (ve * denominator)
      )
      probability <- exp(log_weight - max(log_weight))
      probability <- probability / sum(probability)
      state[marker] <- sample.int(2L, 1L, prob = probability) - 1L
      if (state[marker] > 0L) {
        lhs <- wi[marker] + ve / vbi
        b[marker] <- stats::rnorm(1, xy[marker] / lhs, sqrt(ve / lhs))
      } else b[marker] <- 0
    }
    active <- which(state > 0L)
    pi_active <- stats::rbeta(1, 1 + length(active), 4 + m - length(active))
    square <- if (length(active)) sum(b[active]^2 / q[active]) else 0
    vb <- (square + 4 * vb_prior) /
      stats::rchisq(1, length(active) + 4)
    likelihood <- function(value) {
      if (!length(active)) return(0)
      eta <- drop(X[active, , drop = FALSE] %*% value)
      if (any(!is.finite(eta)) || any(abs(eta) > 700)) return(-Inf)
      -0.5 * sum(eta + b[active]^2 / (vb * exp(eta)))
    }
    if (length(active)) {
      theta <- logvar_bayesc_oracle_ess(theta, likelihood, 0.7)$theta
    } else theta <- stats::rnorm(1, 0, 0.7)
    q <- exp(drop(X %*% theta))
    sse <- max(yy - 2 * sum(b * xy) + wi[1] * sum(b^2), 1e-12)
    ve <- (sse + 4 * ve_prior) / stats::rchisq(1, wi[1] + 4)
    if (iteration > burnin) {
      saved <- saved + 1L
      theta_draw[saved] <- theta
      pi_draw[saved] <- pi_active
      vb_draw[saved] <- vb
      ve_draw[saved] <- ve
      b_sum <- b_sum + b
      pip_sum <- pip_sum + state
      q_sum <- q_sum + q
    }
  }
  list(theta = theta_draw, pi = pi_draw, vb = vb_draw, ve = ve_draw,
       b = b_sum / saved, pip = pip_sum / saved, q = q_sum / saved)
}

logvar_write_diagonal_csr <- function(m) {
  prefix <- tempfile("logvar_bayesc_oracle_csr_")
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), rep(0, m + 1L))
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), integer())
  writeBin(numeric(), paste0(prefix, ".values.f32.bin"), size = 4)
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    paste0("n_variants=", m), "nnz=0", "triangle=upper",
    "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"
  ), paste0(prefix, ".meta.txt"))
  prefix
}

test_that("learned BayesC-LV agrees with the frozen independent R oracle", {
  skip_if_not_installed("coda")
  set.seed(6207)
  m <- 120L
  n <- 900L
  X <- matrix(drop(scale(stats::rnorm(m))), ncol = 1)
  theta_true <- 0.6
  vb_true <- 0.05
  q_true <- exp(drop(X) * theta_true)
  state_true <- stats::rbinom(m, 1, 0.35)
  b_true <- numeric(m)
  b_true[state_true > 0] <- stats::rnorm(
    sum(state_true), 0, sqrt(vb_true * q_true[state_true > 0]))
  ve_true <- sum(b_true^2)
  noise_score <- stats::rnorm(m, 0, sqrt(n * ve_true))
  xy <- n * b_true + noise_score
  yy <- n * sum(b_true^2) + 2 * sum(b_true * noise_score) +
    sum(noise_score^2) / n + ve_true * stats::rchisq(1, n - m)
  iterations <- 1400L
  burnin <- 400L
  seeds_r <- 7101:7104
  oracle <- lapply(seeds_r, function(current_seed) {
    logvar_bayesc_oracle_chain(
      xy, yy, rep(n, m), X, iterations, burnin, current_seed,
      vb_init = vb_true * 0.7, ve_init = ve_true * 1.3,
      vb_prior = vb_true, ve_prior = ve_true)
  })

  prefix <- logvar_write_diagonal_csr(m)
  production <- sblr:::stblr_cpg_omp_csr_logvar_bayesc(
    wy = list(xy), ww = list(rep(n, m)), yy = yy,
    b_init = list(numeric(m)), d_init = list(numeric(m)),
    use_d_init = FALSE, r_init = list(xy), use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE, ld_prefix = prefix,
    B = matrix(vb_true * 0.7, 1, 1), E = matrix(ve_true * 1.3, 1, 1),
    ssb_prior = list(vb_true), sse_prior = list(ve_true),
    pi = c(0.9, 0.1), nub = 4, nue = 4,
    updateB = TRUE, updateE = TRUE, updatePi = TRUE, adjE = 0,
    n = n, nit = iterations - burnin, nburn = burnin, nthin = 1L,
    pi_prior_a = 1, pi_prior_b = 4, ncores = 2L, seed = 8100L,
    nchains = 4L, keep_chains = TRUE, chain_seeds = 8101:8104,
    updateLDswap = FALSE, ld_swap_prob = 0, ld_swap_r2 = 0.8,
    ld_swap_max_friends = 50L, ld_swap_moves = 0L,
    annotation = X, theta_init = matrix(0, 1, 1),
    theta_prior_sd = 0.7, updateTheta = TRUE,
    convergence_markers = integer(), convergence_b = FALSE,
    convergence_d = FALSE
  )

  prod_theta <- production$annotation$theta_trace[
    (burnin + 1L):iterations, 1, , drop = FALSE]
  prod_theta <- matrix(prod_theta, nrow = iterations - burnin, ncol = 4)
  oracle_theta <- do.call(cbind, lapply(oracle, `[[`, "theta"))
  positive_sequence_ess <- function(draws) {
    vapply(seq_len(ncol(draws)), function(chain) {
      value <- draws[, chain]
      correlation <- drop(stats::acf(
        value, lag.max = min(200L, length(value) %/% 2L),
        plot = FALSE, demean = TRUE)$acf)[-1]
      first_negative <- which(correlation < 0)[1]
      if (!is.na(first_negative)) correlation <- correlation[seq_len(first_negative - 1L)]
      length(value) / max(1, 1 + 2 * sum(correlation))
    }, numeric(1))
  }
  prod_ess <- sum(positive_sequence_ess(prod_theta))
  oracle_ess <- sum(positive_sequence_ess(oracle_theta))
  rhat <- function(draws) {
    count <- nrow(draws)
    within <- mean(apply(draws, 2, stats::var))
    between <- count * stats::var(colMeans(draws))
    sqrt((((count - 1) / count) * within + between / count) / within)
  }
  prod_mean <- mean(prod_theta)
  oracle_mean <- mean(oracle_theta)
  mc_tolerance <- 3 * sqrt(
    stats::var(as.numeric(prod_theta)) / prod_ess +
      stats::var(as.numeric(oracle_theta)) / oracle_ess) + 0.08

  expect_lte(abs(prod_mean - oracle_mean), mc_tolerance)
  expect_equal(sd(as.numeric(prod_theta)), sd(as.numeric(oracle_theta)), tolerance = 0.12)
  expect_lt(rhat(prod_theta), 1.15)
  expect_lt(rhat(oracle_theta), 1.15)
  expect_gt(prod_ess, 50)
  expect_gt(oracle_ess, 50)

  oracle_q <- Reduce(`+`, lapply(oracle, `[[`, "q")) / 4
  oracle_b <- Reduce(`+`, lapply(oracle, `[[`, "b")) / 4
  oracle_pip <- Reduce(`+`, lapply(oracle, `[[`, "pip")) / 4
  expect_gt(stats::cor(drop(production$annotation$marker_prior_scale), oracle_q), 0.98)
  expect_gt(stats::cor(drop(production$marker$bm), oracle_b), 0.9)
  expect_gt(stats::cor(drop(production$marker$dm), oracle_pip), 0.75)
  expect_equal(production$pi$mean[1, 2],
               mean(unlist(lapply(oracle, `[[`, "pi"))), tolerance = 0.1)
  expect_equal(mean(production$trace$vbs[(burnin + 1):iterations, 1]),
               mean(unlist(lapply(oracle, `[[`, "vb"))), tolerance = 0.03)
  expect_equal(mean(production$trace$ves[(burnin + 1):iterations, 1]),
               mean(unlist(lapply(oracle, `[[`, "ve"))), tolerance = 0.15 * ve_true)
  q_chain <- production$annotation$marker_prior_scale_chain
  expect_gt(min(stats::cor(q_chain)), 0.98)
  expect_gt(production$diagnostics$logvar$theta_updates, 0)
  expect_gte(production$diagnostics$logvar$mean_likelihood_evaluations_per_update, 1)
})
