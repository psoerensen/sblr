repo_root <- normalizePath(Sys.getenv("SBLR_ROOT", "."), winslash = "/",
                           mustWork = TRUE)
sys.source(file.path(repo_root, "R", "sbayesrc-block-px-reference.R"),
           envir = globalenv())
output_root <- file.path(repo_root, "results", "local", "sbayesrc_block_px")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

fixture <- local({
  annotation <- cbind(intercept = 1, enriched = c(0, 0, 0, 0, 1, 1, 1, 1))
  list(marker_count = 8L, annotation = annotation,
       diagonal = c(1.6, 1.8, 2.0, 2.2, 1.7, 1.9, 2.1, 2.3),
       score = c(-0.15, 0.2, 0.05, -0.3, 1.25, 0.95, 1.1, 0.8),
       residual_variance = 1, marker_variance = 0.7,
       intercept_mean = qnorm(0.15), intercept_precision = 1,
       sigma_sq_alpha = 1)
})

reference <- jsonlite::read_json(file.path(
  repo_root, "results", "local", "bayesrc_coordinated_transition",
  "oracle_result.json"), simplifyVector = TRUE)$exact

bayes_factor <- sqrt(fixture$residual_variance /
  (fixture$residual_variance + fixture$diagonal * fixture$marker_variance)) *
  exp(0.5 * fixture$score^2 * fixture$marker_variance /
    (fixture$residual_variance *
      (fixture$residual_variance + fixture$diagonal * fixture$marker_variance)))
active_mean <- fixture$score /
  (fixture$diagonal + fixture$residual_variance / fixture$marker_variance)
active_variance <- fixture$residual_variance /
  (fixture$diagonal + fixture$residual_variance / fixture$marker_variance)

draw_truncated <- function(mean, positive) {
  if (positive) {
    p0 <- pnorm(-mean)
    mean + qnorm(runif(1L, p0, 1))
  } else {
    p0 <- pnorm(-mean)
    mean + qnorm(runif(1L, 0, p0))
  }
}

draw_marker_state <- function(alpha) {
  continuation <- pnorm(drop(fixture$annotation %*% alpha))
  active_probability <- continuation * bayes_factor /
    (1 - continuation + continuation * bayes_factor)
  component <- as.integer(runif(fixture$marker_count) < active_probability)
  beta <- numeric(fixture$marker_count)
  active <- which(component == 1L)
  beta[active] <- rnorm(length(active), active_mean[active],
                        sqrt(active_variance[active]))
  list(component = component, beta = beta)
}

draw_latent <- function(component, alpha) {
  mean <- drop(fixture$annotation %*% alpha)
  vapply(seq_len(fixture$marker_count), function(marker) {
    draw_truncated(mean[[marker]], component[[marker]] == 1L)
  }, numeric(1L))
}

draw_alpha_scalar <- function(latent, alpha) {
  residual <- latent - drop(fixture$annotation %*% alpha)
  for (coefficient in seq_along(alpha)) {
    x <- fixture$annotation[, coefficient]
    prior_precision <- if (coefficient == 1L)
      fixture$intercept_precision else 1 / fixture$sigma_sq_alpha
    prior_mean <- if (coefficient == 1L) fixture$intercept_mean else 0
    rhs <- sum(x * residual) + sum(x^2) * alpha[[coefficient]] +
      prior_precision * prior_mean
    variance <- 1 / (sum(x^2) + prior_precision)
    updated <- rnorm(1L, variance * rhs, sqrt(variance))
    residual <- residual + x * (alpha[[coefficient]] - updated)
    alpha[[coefficient]] <- updated
  }
  alpha
}

draw_alpha_blocked <- function(latent) {
  conditional <- .sbayesrc_px_alpha_conditional(
    latent, fixture$annotation,
    c(fixture$intercept_mean, 0),
    c(fixture$intercept_precision, 1 / fixture$sigma_sq_alpha))
  drop(conditional$mean + t(chol(conditional$covariance)) %*% rnorm(2L))
}

run_chain <- function(seed, kernel, iterations = 30000L, burnin = 5000L,
                      thinning = 5L, proposal_sd = 0.45) {
  set.seed(seed)
  alpha <- c(fixture$intercept_mean, 0)
  retained <- (iterations - burnin) %/% thinning
  alpha_trace <- matrix(NA_real_, retained, 2L)
  component_trace <- matrix(NA_integer_, retained, fixture$marker_count)
  beta_trace <- matrix(NA_real_, retained, fixture$marker_count)
  attempted <- accepted <- 0L
  scale_sum <- alpha_jump_sum <- 0
  draw <- 0L
  for (iteration in seq_len(iterations)) {
    marker_state <- draw_marker_state(alpha)
    latent <- draw_latent(marker_state$component, alpha)
    old_alpha <- alpha
    if (kernel == "ordinary") {
      alpha <- draw_alpha_scalar(latent, alpha)
    } else {
      attempted <- attempted + 1L
      log_scale <- rnorm(1L, 0, proposal_sd)
      log_ratio <- .sbayesrc_px_log_scale_ratio(
        log_scale, latent, fixture$annotation,
        c(fixture$intercept_mean, 0),
        c(fixture$intercept_precision, 1 / fixture$sigma_sq_alpha))
      if (log(runif(1L)) < min(0, log_ratio)) {
        latent <- exp(log_scale) * latent
        accepted <- accepted + 1L
        scale_sum <- scale_sum + abs(log_scale)
      }
      alpha <- draw_alpha_blocked(latent)
    }
    alpha_jump_sum <- alpha_jump_sum + sqrt(sum((alpha - old_alpha)^2))
    if (iteration > burnin && (iteration - burnin) %% thinning == 0L) {
      draw <- draw + 1L
      alpha_trace[draw, ] <- alpha
      component_trace[draw, ] <- marker_state$component
      beta_trace[draw, ] <- marker_state$beta
    }
  }
  list(alpha = alpha_trace, component = component_trace, beta = beta_trace,
       attempts = attempted, acceptances = accepted,
       acceptance_rate = if (attempted) accepted / attempted else NA_real_,
       mean_accepted_abs_log_scale = if (accepted) scale_sum / accepted else NA_real_,
       mean_alpha_jump = alpha_jump_sum / iterations)
}

summarise <- function(chains) {
  alpha <- do.call(rbind, lapply(chains, `[[`, "alpha"))
  component <- do.call(rbind, lapply(chains, `[[`, "component"))
  beta <- do.call(rbind, lapply(chains, `[[`, "beta"))
  pip <- colMeans(component)
  active <- rowSums(component)
  active_distribution <- tabulate(active + 1L,
    nbins = fixture$marker_count + 1L) / length(active)
  list(
    alpha_mean = colMeans(alpha),
    alpha_variance = apply(alpha, 2L, var),
    alpha_quantile = apply(alpha, 2L, quantile, c(.025, .5, .975)),
    pip = pip,
    beta_mean = colMeans(beta),
    active_distribution = active_distribution,
    alpha_chain_mcse = apply(vapply(chains, function(x) colMeans(x$alpha),
      numeric(2L)), 1L, sd) / sqrt(length(chains)),
    pip_chain_mcse = apply(vapply(chains, function(x) colMeans(x$component),
      numeric(fixture$marker_count)), 1L, sd) / sqrt(length(chains)),
    attempts = sum(vapply(chains, `[[`, integer(1L), "attempts")),
    acceptances = sum(vapply(chains, `[[`, integer(1L), "acceptances")),
    acceptance_rate = mean(vapply(chains, `[[`, numeric(1L), "acceptance_rate"),
      na.rm = TRUE),
    mean_accepted_abs_log_scale = mean(vapply(chains, `[[`, numeric(1L),
      "mean_accepted_abs_log_scale"), na.rm = TRUE),
    mean_alpha_jump = mean(vapply(chains, `[[`, numeric(1L), "mean_alpha_jump"))
  )
}

compare_reference <- function(summary) {
  alpha_mcse <- pmax(summary$alpha_chain_mcse, 0.01)
  pip_mcse <- pmax(summary$pip_chain_mcse, 0.0025)
  list(
    alpha_error = summary$alpha_mean - reference$alpha_mean,
    alpha_tolerance = 4 * alpha_mcse,
    pip_max_error = max(abs(summary$pip - reference$pip)),
    pip_max_tolerance = 4 * max(pip_mcse),
    beta_max_error = max(abs(summary$beta_mean - reference$beta_mean)),
    active_distribution_l1 = sum(abs(
      summary$active_distribution - reference$active_distribution))
  )
}

seeds <- c(22011L, 22022L, 22033L, 22044L)
started <- proc.time()[[3L]]
ordinary <- lapply(seeds, run_chain, kernel = "ordinary")
px <- lapply(seeds, run_chain, kernel = "px")
ordinary_summary <- summarise(ordinary)
px_summary <- summarise(px)
result <- list(
  schema = "sblr-block-px-oracle-v1",
  package_sha = system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
                        stdout = TRUE),
  seeds = seeds,
  settings = list(iterations = 30000L, burnin = 5000L, thinning = 5L,
                  proposal_sd = 0.45),
  exact = reference,
  ordinary = ordinary_summary,
  px = px_summary,
  comparison = list(ordinary = compare_reference(ordinary_summary),
                    px = compare_reference(px_summary)),
  runtime_seconds = proc.time()[[3L]] - started)
jsonlite::write_json(result, file.path(output_root, "px_oracle_result.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 16)
saveRDS(list(ordinary = ordinary, px = px),
        file.path(output_root, "px_oracle_chains.rds"), compress = "xz")
print(result[c("ordinary", "px", "comparison", "runtime_seconds")])
