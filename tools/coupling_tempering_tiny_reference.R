#!/usr/bin/env Rscript

# Independent near-exact reference for the development-only BayesRC coupling
# tempering kernel. The tolerances below are declared before sampling.
tiny_tolerance <- list(
  quadrature_pip = 2e-4,
  quadrature_active = 5e-4,
  monte_carlo_pip = 0.035,
  monte_carlo_active = 0.15,
  monte_carlo_alpha = 0.18
)

write_tiny_bed <- function(dosage) {
  path <- tempfile(fileext = ".bed")
  code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    z <- unname(code[as.character(dosage[marker, ])])
    z <- c(z, rep(0L, (-length(z)) %% 4L))
    vapply(seq(1L, length(z), 4L), function(first) {
      sum(z[first:(first + 3L)] * c(1L, 4L, 16L, 64L))
    }, integer(1L))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
  path
}

log_mvn_zero <- function(y, covariance) {
  factor <- chol(covariance)
  solved <- backsolve(factor, y, transpose = TRUE)
  -0.5 * (length(y) * log(2 * pi) +
    2 * sum(log(diag(factor))) + sum(solved^2))
}

allocation_grid <- function(markers) {
  t(vapply(0:(2^markers - 1L), function(value) {
    as.integer(intToBits(value))[seq_len(markers)]
  }, integer(markers)))
}

tiny_exact_reference <- function(dosage, y, annotation, mixture, vb, ve,
                                 quadrature_order) {
  stopifnot(requireNamespace("statmod", quietly = TRUE))
  markers <- nrow(dosage)
  af <- rowMeans(dosage) / 2
  scale <- sqrt(2 * af * (1 - af))
  X <- t((dosage - 2 * af) / scale)
  allocation <- allocation_grid(markers)
  mu <- qnorm(mixture[2L])
  quadrature <- statmod::gauss.quad.prob(quadrature_order, dist = "normal")
  nodes <- expand.grid(intercept = quadrature$nodes, slope = quadrature$nodes)
  weights <- as.vector(outer(quadrature$weights, quadrature$weights))
  nodes$intercept <- nodes$intercept + mu

  log_likelihood <- apply(allocation, 1L, function(state) {
    active <- which(state == 1L)
    covariance <- diag(ve, length(y))
    if (length(active)) {
      covariance <- covariance + vb * tcrossprod(X[, active, drop = FALSE])
    }
    log_mvn_zero(y, covariance)
  })

  integrated_prior <- numeric(nrow(allocation))
  alpha_numerator <- matrix(0, nrow(allocation), 2L)
  eta <- outer(nodes$intercept, rep(1, markers)) +
    outer(nodes$slope, annotation[, 2L])
  log_active_probability <- pnorm(eta, log.p = TRUE)
  log_null_probability <- pnorm(eta, lower.tail = FALSE, log.p = TRUE)
  for (state_index in seq_len(nrow(allocation))) {
    state <- allocation[state_index, ]
    log_terms <- log_null_probability
    active <- which(state == 1L)
    if (length(active)) {
      log_terms[, active] <- log_active_probability[, active, drop = FALSE]
    }
    conditional <- exp(rowSums(log_terms))
    weighted <- weights * conditional
    integrated_prior[state_index] <- sum(weighted)
    alpha_numerator[state_index, ] <- colSums(weighted * as.matrix(nodes))
  }
  log_weight <- log_likelihood + log(integrated_prior)
  posterior <- exp(log_weight - max(log_weight))
  posterior <- posterior / sum(posterior)
  alpha_conditional <- alpha_numerator / integrated_prior
  list(
    allocation = allocation,
    probability = posterior,
    pip = colSums(posterior * allocation),
    active = sum(posterior * rowSums(allocation)),
    alpha = colSums(posterior * alpha_conditional),
    active_probability = rowsum(posterior, rowSums(allocation), reorder = FALSE),
    X = X
  )
}

tiny_native_args <- function(dosage, y, annotation, mixture, vb, ve, seed,
                             iterations) {
  bed <- write_tiny_bed(dosage)
  markers <- nrow(dosage)
  baseline <- sblr:::.bayesr_pi_to_probit_stick_intercepts(mixture)
  prior <- sblr:::.sbayesrc_resolve_intercept_prior(mixture)$native
  prior <- rbind(
    prior,
    update_sigmaSqAlpha = rep(0, ncol(prior)),
    allocation_updates_per_cycle = rep(1, ncol(prior)),
    annotation_updates_per_cycle = rep(1, ncol(prior)),
    coupling_tempering = rep(1, ncol(prior)),
    coupling_swap_every = rep(5, ncol(prior))
  )
  list(
    bed_files = bed, n = ncol(dosage), cls = list(seq_len(markers)), y = matrix(y),
    b_init = list(numeric(markers)), sets = rep(1L, markers), rows = NULL,
    af = list(rowMeans(dosage) / 2), scale = TRUE, B = matrix(vb), E = matrix(ve),
    ssb_prior = list(0.05), sse_prior = list(0.5), A = annotation,
    gamma = c(0, 1), annot_alpha_init = rbind(baseline, 0),
    annot_sigma_sq_alpha_init = 1, intercept_prior_resolved = prior,
    updateAlpha = TRUE, annot_alpha_update_every = 1L,
    updateB = FALSE, updateE = FALSE, adjE = 0,
    nit = iterations, nburn = 0L, nthin = 1L, rebuild_every = 20L,
    return_wy = FALSE, return_r = FALSE, read_block_size = 8L,
    nchains = 1L, keep_chains = TRUE, ncores = 1L, seed = seed,
    convergence_markers = 0:(markers - 1L), convergence_annotations = TRUE,
    convergence_b = TRUE, convergence_d = TRUE, convergence_component = TRUE
  )
}

summarise_tiny_native <- function(raw, burnin) {
  chain <- raw$chains[[1L]][[1L]]
  trace <- chain$convergence_trace
  keep <- seq.int(burnin + 1L, nrow(trace$component))
  component <- trace$component[keep, , drop = FALSE]
  list(
    pip = colMeans(component > 0L),
    active = mean(rowSums(component > 0L)),
    alpha = colMeans(trace$alpha[keep, , drop = FALSE]),
    swaps = chain$coupling_tempering$swap,
    identities = chain$coupling_tempering$replica_identity[keep, , drop = FALSE]
  )
}

run_tiny_coupling_reference <- function(output_root = file.path(
    "results", "local", "study06_bed_coupling_tempering_screen", "tiny")) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 0, 1, 2, 1, 0, 2),
    c(2, 1, 0, 2, 1, 0, 2, 1, 0, 1, 2, 0),
    c(0, 0, 1, 1, 2, 2, 0, 1, 2, 2, 1, 0),
    c(2, 2, 1, 1, 0, 0, 2, 1, 0, 0, 1, 2),
    c(0, 1, 0, 1, 2, 1, 2, 0, 2, 1, 0, 2),
    c(2, 1, 2, 1, 0, 1, 0, 2, 0, 1, 2, 0),
    c(0, 2, 1, 0, 2, 1, 1, 0, 2, 1, 2, 0),
    c(2, 0, 1, 2, 0, 1, 1, 2, 0, 1, 0, 2)
  )
  annotation <- cbind(intercept = 1, binary = rep(c(0, 1), each = 4L))
  mixture <- c(0.8, 0.2)
  vb <- 0.20
  ve <- 1
  y <- c(-0.7, -0.3, 1.0, -0.4, 0.9, 0.2, -0.8, 0.4, 1.1, 0.1, -0.6, 0.5)

  exact31 <- tiny_exact_reference(dosage, y, annotation, mixture, vb, ve, 31L)
  exact41 <- tiny_exact_reference(dosage, y, annotation, mixture, vb, ve, 41L)
  quadrature_error <- c(
    pip = max(abs(exact31$pip - exact41$pip)),
    active = abs(exact31$active - exact41$active),
    alpha = max(abs(exact31$alpha - exact41$alpha))
  )
  quadrature_pass <- quadrature_error[["pip"]] <= tiny_tolerance$quadrature_pip &&
    quadrature_error[["active"]] <= tiny_tolerance$quadrature_active

  seeds <- c(81101L, 81202L, 81303L, 81404L)
  iterations <- 50000L
  burnin <- 10000L
  draws <- lapply(seeds, function(seed) {
    args <- tiny_native_args(dosage, y, annotation, mixture, vb, ve, seed, iterations)
    summarise_tiny_native(do.call(sblr:::.stblr_bed_bayesrc_native, args), burnin)
  })
  native <- list(
    pip = Reduce(`+`, lapply(draws, `[[`, "pip")) / length(draws),
    active = mean(vapply(draws, `[[`, numeric(1L), "active")),
    alpha = Reduce(`+`, lapply(draws, `[[`, "alpha")) / length(draws)
  )
  monte_carlo_error <- c(
    pip = max(abs(native$pip - exact41$pip)),
    active = abs(native$active - exact41$active),
    alpha = max(abs(native$alpha - exact41$alpha))
  )
  monte_carlo_pass <- monte_carlo_error[["pip"]] <= tiny_tolerance$monte_carlo_pip &&
    monte_carlo_error[["active"]] <= tiny_tolerance$monte_carlo_active &&
    monte_carlo_error[["alpha"]] <= tiny_tolerance$monte_carlo_alpha
  acceptance <- do.call(rbind, lapply(seq_along(draws), function(index) {
    swap <- draws[[index]]$swaps
    data.frame(seed = seeds[index], pair = swap[, 2L], accepted = swap[, 3L],
               probability = swap[, 4L])
  }))
  result <- list(
    specification = list(markers = nrow(dosage), samples = ncol(dosage),
      components = 2L, sticks = 1L, mixture = mixture, vb = vb, ve = ve),
    tolerance = tiny_tolerance, exact = exact41[c("pip", "active", "alpha")],
    native = native, quadrature_error = quadrature_error,
    monte_carlo_error = monte_carlo_error,
    quadrature_pass = quadrature_pass, monte_carlo_pass = monte_carlo_pass,
    seeds = seeds, iterations = iterations, burnin = burnin,
    swap_summary = aggregate(cbind(accepted, probability) ~ pair, acceptance, mean)
  )
  saveRDS(result, file.path(output_root, "tiny_reference_result.rds"))
  write.csv(data.frame(marker = seq_len(nrow(dosage)), exact_pip = exact41$pip,
    native_pip = native$pip), file.path(output_root, "tiny_pip_comparison.csv"),
    row.names = FALSE)
  print(result)
  if (!quadrature_pass || !monte_carlo_pass) {
    stop("Tiny coupling-tempering posterior validation failed.")
  }
  invisible(result)
}

if (sys.nframe() == 0L) {
  devtools::load_all(quiet = TRUE)
  run_tiny_coupling_reference()
}
