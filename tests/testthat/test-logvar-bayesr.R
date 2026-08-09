logvar_bayesr_fixture <- function() {
  prefix <- tempfile("logvar-bayesr-")
  m <- 4L
  writeLines(c(paste0("n_variants=", m), "nnz=0"),
             paste0(prefix, ".meta.txt"))
  writeBin(rep(as.raw(0), 8L * (m + 1L)),
           paste0(prefix, ".row_ptr.u64.bin"))
  file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
  file.create(paste0(prefix, ".values.f32.bin"))
  list(prefix = prefix, m = m)
}

logvar_bayesr_diagonal_csr <- function(m) {
  fixture <- logvar_bayesr_fixture()
  if (fixture$m != m) {
    unlink(paste0(fixture$prefix, c(
      ".meta.txt", ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
      ".values.f32.bin")))
    prefix <- tempfile("logvar-bayesr-oracle-")
    writeLines(c(paste0("n_variants=", m), "nnz=0"),
               paste0(prefix, ".meta.txt"))
    writeBin(rep(as.raw(0), 8L * (m + 1L)),
             paste0(prefix, ".row_ptr.u64.bin"))
    file.create(paste0(prefix, ".col_idx.u32.0based.bin"))
    file.create(paste0(prefix, ".values.f32.bin"))
    return(prefix)
  }
  fixture$prefix
}

logvar_bayesr_common <- function(fixture) {
  list(
    wy = list(c(4, -3, 2, 1)), ww = list(rep(79, fixture$m)), yy = 79,
    b_init = list(rep(0, fixture$m)), comp_init = list(rep(0, fixture$m)),
    use_comp_init = TRUE, r_init = list(c(4, -3, 2, 1)),
    use_r_init = TRUE, rebuild_r_before_updateE = FALSE,
    ld_prefix = fixture$prefix, B = matrix(0.04), E = matrix(0.7),
    ssb_prior = list(0.04), sse_prior = list(0.7),
    pi = c(0.7, 0.2, 0.1), mixture_var = c(0, 0.1, 1),
    alpha = c(1, 1, 1), nub = 4, nue = 4,
    updateB = TRUE, updateE = TRUE, updatePi = TRUE, adjE = 0,
    n = 80L, nit = 20L, nburn = 5L, nthin = 1L, ncores = 1L,
    seed = 920L, nchains = 2L, keep_chains = TRUE,
    chain_seeds = 921:922, updateE_start = 0L, updateE_every = 1L,
    updateLDswap = FALSE, ld_swap_prob = 0, ld_swap_r2 = 0.8,
    ld_swap_max_friends = 50L, ld_swap_moves = 0L,
    convergence_markers = integer(), convergence_probability = FALSE,
    convergence_b = FALSE, convergence_d = FALSE,
    convergence_component = FALSE
  )
}

logvar_bayesr_trajectory <- function(raw) {
  diagnostics <- raw$diagnostics
  diagnostics$logvar <- NULL
  list(
    marker = raw$marker, trace = raw$trace, variance = raw$variance,
    pi = raw$pi, component = raw$component, chains = raw$chains,
    diagnostics = diagnostics
  )
}

test_that("BayesR-LV theta zero is an exact ordinary BayesR trajectory", {
  fixture <- logvar_bayesr_fixture()
  on.exit(unlink(paste0(fixture$prefix, c(
    ".meta.txt", ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin"))), add = TRUE)
  common <- logvar_bayesr_common(fixture)
  ordinary <- do.call(stblr_cpg_omp_csr_bayesr, c(common, list(
    maf_effect_s_prior_scale = NULL, estimate_maf_effect_s = FALSE,
    maf_effect_s_init = 0, maf_effect_s_prior = c(-3, 2),
    maf_effect_s_proposal_sd = 0.35, maf_effect_s_log_h = NULL)))
  lv <- do.call(stblr_cpg_omp_csr_logvar_bayesr, c(common, list(
    annotation = matrix(c(-1.5, -0.5, 0.5, 1.5), ncol = 1),
    theta_init = matrix(0), theta_prior_sd = 0.7, updateTheta = FALSE)))
  expect_identical(logvar_bayesr_trajectory(lv),
                   logvar_bayesr_trajectory(ordinary))
  expect_identical(drop(lv$annotation$theta), 0)
  expect_identical(drop(lv$annotation$marker_prior_scale), rep(1, fixture$m))
  expect_identical(lv$diagnostics$logvar$theta_updates, 0)
})

test_that("BayesR-LV fixed theta is exact fixed marker prior scale BayesR", {
  fixture <- logvar_bayesr_fixture()
  on.exit(unlink(paste0(fixture$prefix, c(
    ".meta.txt", ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin"))), add = TRUE)
  common <- logvar_bayesr_common(fixture)
  X <- matrix(c(-1.5, -0.5, 0.5, 1.5), ncol = 1)
  theta <- 0.3
  q <- drop(exp(X * theta))
  fixed <- do.call(stblr_cpg_omp_csr_bayesr, c(common, list(
    maf_effect_s_prior_scale = q, estimate_maf_effect_s = FALSE,
    maf_effect_s_init = 0, maf_effect_s_prior = c(-3, 2),
    maf_effect_s_proposal_sd = 0.35, maf_effect_s_log_h = NULL)))
  lv <- do.call(stblr_cpg_omp_csr_logvar_bayesr, c(common, list(
    annotation = X, theta_init = matrix(theta), theta_prior_sd = 0.7,
    updateTheta = FALSE)))
  expect_identical(logvar_bayesr_trajectory(lv),
                   logvar_bayesr_trajectory(fixed))
  expect_equal(as.numeric(lv$annotation$marker_prior_scale), q,
               tolerance = 1e-15)
})

test_that("BayesR-LV learned theta returns finite shared-kernel diagnostics", {
  fixture <- logvar_bayesr_fixture()
  on.exit(unlink(paste0(fixture$prefix, c(
    ".meta.txt", ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin"))), add = TRUE)
  common <- logvar_bayesr_common(fixture)
  lv <- do.call(stblr_cpg_omp_csr_logvar_bayesr, c(common, list(
    annotation = matrix(c(-1.5, -0.5, 0.5, 1.5), ncol = 1),
    theta_init = matrix(0), theta_prior_sd = 0.7, updateTheta = TRUE)))
  expect_true(all(is.finite(lv$annotation$theta_trace)))
  expect_true(all(is.finite(lv$annotation$marker_prior_scale)))
  expect_gt(lv$diagnostics$logvar$theta_updates, 0)
  expect_gte(lv$diagnostics$logvar$max_likelihood_evaluations, 1)
  expect_lte(lv$diagnostics$logvar$min_log_q,
             lv$diagnostics$logvar$max_log_q)
})

logvar_bayesr_oracle_ess <- function(theta, loglik, prior_sd) {
  current <- loglik(theta)
  nu <- stats::rnorm(length(theta), 0, prior_sd)
  threshold <- current + log(stats::runif(1))
  angle <- stats::runif(1, 0, 2 * pi)
  lower <- angle - 2 * pi
  upper <- angle
  for (step in seq_len(10000L)) {
    proposal <- theta * cos(angle) + nu * sin(angle)
    value <- loglik(proposal)
    if (is.finite(value) && value > threshold) return(proposal)
    if (angle < 0) lower <- angle else upper <- angle
    angle <- stats::runif(1, lower, upper)
  }
  stop("BayesR test oracle ESS exceeded its bracket guard")
}

logvar_bayesr_oracle_chain <- function(
    xy, yy, wi, X, gamma, iterations, burnin, seed, vb_init, ve_init,
    vb_prior, ve_prior) {
  set.seed(seed)
  m <- length(xy)
  K <- length(gamma)
  theta <- 0
  q <- rep(1, m)
  b <- numeric(m)
  component <- integer(m)
  pi_mix <- rep(1 / K, K)
  theta_draw <- numeric(iterations - burnin)
  pi_draw <- matrix(0, iterations - burnin, K)
  vb_draw <- ve_draw <- theta_draw
  b_sum <- pip_sum <- q_sum <- numeric(m)
  component_sum <- matrix(0, m, K)
  vb <- vb_init
  ve <- ve_init
  saved <- 0L
  for (iteration in seq_len(iterations)) {
    for (marker in sample.int(m)) {
      log_weight <- rep(-Inf, K)
      log_weight[1] <- log(max(pi_mix[1], 1e-300))
      for (k in 2:K) {
        vk <- vb * gamma[k] * q[marker]
        denominator <- ve + wi[marker] * vk
        log_weight[k] <- log(max(pi_mix[k], 1e-300)) +
          0.5 * log(ve / denominator) +
          0.5 * xy[marker]^2 * vk / (ve * denominator)
      }
      probability <- exp(log_weight - max(log_weight))
      probability <- probability / sum(probability)
      component[marker] <- sample.int(K, 1L, prob = probability) - 1L
      if (component[marker] > 0L) {
        vk <- vb * gamma[component[marker] + 1L] * q[marker]
        lhs <- wi[marker] + ve / vk
        b[marker] <- stats::rnorm(1, xy[marker] / lhs, sqrt(ve / lhs))
      } else b[marker] <- 0
    }
    counts <- tabulate(component + 1L, nbins = K)
    pi_gamma <- stats::rgamma(K, shape = 1 + counts, rate = 1)
    pi_mix <- pi_gamma / sum(pi_gamma)
    active <- which(component > 0L)
    square <- if (length(active)) sum(
      b[active]^2 /
        (gamma[component[active] + 1L] * q[active])) else 0
    vb <- (square + 4 * vb_prior) /
      stats::rchisq(1, length(active) + 4)
    likelihood <- function(value) {
      if (!length(active)) return(0)
      eta <- drop(X[active, , drop = FALSE] %*% value)
      if (any(!is.finite(eta)) || any(abs(eta) > 700)) return(-Inf)
      -0.5 * sum(eta + b[active]^2 /
        (vb * gamma[component[active] + 1L] * exp(eta)))
    }
    if (length(active)) {
      theta <- logvar_bayesr_oracle_ess(theta, likelihood, 0.7)
    } else theta <- stats::rnorm(1, 0, 0.7)
    q <- exp(drop(X %*% theta))
    sse <- max(yy - 2 * sum(b * xy) + wi[1] * sum(b^2), 1e-12)
    ve <- (sse + 4 * ve_prior) / stats::rchisq(1, wi[1] + 4)
    if (iteration > burnin) {
      saved <- saved + 1L
      theta_draw[saved] <- theta
      pi_draw[saved, ] <- pi_mix
      vb_draw[saved] <- vb
      ve_draw[saved] <- ve
      b_sum <- b_sum + b
      pip_sum <- pip_sum + (component > 0L)
      q_sum <- q_sum + q
      for (k in seq_len(K))
        component_sum[, k] <- component_sum[, k] +
          (component == (k - 1L))
    }
  }
  list(theta = theta_draw, pi = pi_draw, vb = vb_draw, ve = ve_draw,
       b = b_sum / saved, pip = pip_sum / saved, q = q_sum / saved,
       component = component_sum / saved)
}

test_that("learned BayesR-LV agrees with the frozen independent R oracle", {
  set.seed(6307)
  m <- 120L
  n <- 900L
  gamma <- c(0, 0.1, 1)
  X <- matrix(drop(scale(stats::rnorm(m))), ncol = 1)
  theta_true <- 0.6
  vb_true <- 0.05
  q_true <- exp(drop(X) * theta_true)
  component_true <- sample.int(3L, m, replace = TRUE,
                               prob = c(0.55, 0.3, 0.15)) - 1L
  b_true <- numeric(m)
  active <- which(component_true > 0L)
  b_true[active] <- stats::rnorm(
    length(active), 0,
    sqrt(vb_true * gamma[component_true[active] + 1L] * q_true[active]))
  ve_true <- sum(b_true^2)
  noise_score <- stats::rnorm(m, 0, sqrt(n * ve_true))
  xy <- n * b_true + noise_score
  yy <- n * sum(b_true^2) + 2 * sum(b_true * noise_score) +
    sum(noise_score^2) / n + ve_true * stats::rchisq(1, n - m)
  iterations <- 1400L
  burnin <- 400L
  oracle <- lapply(7201:7204, function(current_seed)
    logvar_bayesr_oracle_chain(
      xy, yy, rep(n, m), X, gamma, iterations, burnin, current_seed,
      vb_init = vb_true * 0.7, ve_init = ve_true * 1.3,
      vb_prior = vb_true, ve_prior = ve_true))

  prefix <- logvar_bayesr_diagonal_csr(m)
  on.exit(unlink(paste0(prefix, c(
    ".meta.txt", ".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin"))), add = TRUE)
  production <- stblr_cpg_omp_csr_logvar_bayesr(
    wy = list(xy), ww = list(rep(n, m)), yy = yy,
    b_init = list(numeric(m)), comp_init = list(numeric(m)),
    use_comp_init = TRUE, r_init = list(xy), use_r_init = FALSE,
    rebuild_r_before_updateE = FALSE, ld_prefix = prefix,
    B = matrix(vb_true * 0.7), E = matrix(ve_true * 1.3),
    ssb_prior = list(vb_true), sse_prior = list(ve_true),
    pi = rep(1 / 3, 3), mixture_var = gamma, alpha = rep(1, 3),
    nub = 4, nue = 4, updateB = TRUE, updateE = TRUE,
    updatePi = TRUE, adjE = 0, n = n,
    nit = iterations - burnin, nburn = burnin, nthin = 1L,
    ncores = 2L, seed = 8200L, nchains = 4L, keep_chains = TRUE,
    chain_seeds = 8201:8204, updateE_start = 0L, updateE_every = 1L,
    updateLDswap = FALSE, ld_swap_prob = 0, ld_swap_r2 = 0.8,
    ld_swap_max_friends = 50L, ld_swap_moves = 0L,
    annotation = X, theta_init = matrix(0), theta_prior_sd = 0.7,
    updateTheta = TRUE, convergence_markers = integer(),
    convergence_probability = FALSE, convergence_b = FALSE,
    convergence_d = FALSE, convergence_component = FALSE)

  prod_theta <- matrix(production$annotation$theta_trace[
    (burnin + 1L):iterations, 1, , drop = FALSE],
    nrow = iterations - burnin, ncol = 4)
  oracle_theta <- do.call(cbind, lapply(oracle, `[[`, "theta"))
  positive_sequence_ess <- function(draws) {
    vapply(seq_len(ncol(draws)), function(chain) {
      correlation <- drop(stats::acf(
        draws[, chain], lag.max = min(200L, nrow(draws) %/% 2L),
        plot = FALSE, demean = TRUE)$acf)[-1]
      first_negative <- which(correlation < 0)[1]
      if (!is.na(first_negative))
        correlation <- correlation[seq_len(first_negative - 1L)]
      nrow(draws) / max(1, 1 + 2 * sum(correlation))
    }, numeric(1))
  }
  rhat <- function(draws) {
    count <- nrow(draws)
    within <- mean(apply(draws, 2, stats::var))
    between <- count * stats::var(colMeans(draws))
    sqrt((((count - 1) / count) * within + between / count) / within)
  }
  prod_ess <- sum(positive_sequence_ess(prod_theta))
  oracle_ess <- sum(positive_sequence_ess(oracle_theta))
  tolerance <- 3 * sqrt(
    stats::var(as.numeric(prod_theta)) / prod_ess +
      stats::var(as.numeric(oracle_theta)) / oracle_ess) + 0.1
  expect_lte(abs(mean(prod_theta) - mean(oracle_theta)), tolerance)
  expect_equal(sd(as.numeric(prod_theta)), sd(as.numeric(oracle_theta)),
               tolerance = 0.15)
  expect_lt(rhat(prod_theta), 1.15)
  expect_lt(rhat(oracle_theta), 1.15)
  expect_gt(prod_ess, 50)
  expect_gt(oracle_ess, 50)

  oracle_q <- Reduce(`+`, lapply(oracle, `[[`, "q")) / 4
  oracle_b <- Reduce(`+`, lapply(oracle, `[[`, "b")) / 4
  oracle_pip <- Reduce(`+`, lapply(oracle, `[[`, "pip")) / 4
  oracle_component <- Reduce(`+`, lapply(oracle, `[[`, "component")) / 4
  expect_gt(stats::cor(drop(production$annotation$marker_prior_scale),
                       oracle_q), 0.98)
  expect_gt(stats::cor(drop(production$marker$bm), oracle_b), 0.85)
  expect_gt(stats::cor(drop(production$marker$dm), oracle_pip), 0.7)
  expect_gt(stats::cor(as.numeric(production$component$prob[[1]]),
                       as.numeric(oracle_component)), 0.7)
  expect_equal(as.numeric(production$pi$mean),
               colMeans(do.call(rbind, lapply(oracle, `[[`, "pi"))),
               tolerance = 0.1)
  expect_equal(mean(production$trace$vbs[(burnin + 1):iterations, 1]),
               mean(unlist(lapply(oracle, `[[`, "vb"))), tolerance = 0.04)
  expect_equal(mean(production$trace$ves[(burnin + 1):iterations, 1]),
               mean(unlist(lapply(oracle, `[[`, "ve"))),
               tolerance = 0.15 * ve_true)
  expect_gt(min(stats::cor(
    production$annotation$marker_prior_scale_chain)), 0.98)
  expect_gt(production$diagnostics$logvar$theta_updates, 0)
  expect_gte(
    production$diagnostics$logvar$mean_likelihood_evaluations_per_update, 1)
})
