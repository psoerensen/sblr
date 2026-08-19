repo_root <- normalizePath(Sys.getenv("SBLR_ROOT", "."), winslash = "/",
                           mustWork = TRUE)
sys.source(file.path(repo_root, "R", "bayesrc-coordinated-transition-reference.R"),
           envir = globalenv())
output_root <- file.path(repo_root, "results", "local",
                         "bayesrc_coordinated_transition")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

weighted_quantile <- function(x, weight, probability) {
  order <- order(x)
  x <- x[order]
  cumulative <- cumsum(weight[order]) / sum(weight)
  vapply(probability, function(p) x[which(cumulative >= p)[1L]], numeric(1L))
}

log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

fixture <- local({
  annotation <- cbind(intercept = 1, enriched = c(0, 0, 0, 0, 1, 1, 1, 1))
  list(marker_count = 8L, annotation = annotation,
       diagonal = c(1.6, 1.8, 2.0, 2.2, 1.7, 1.9, 2.1, 2.3),
       score = c(-0.15, 0.2, 0.05, -0.3, 1.25, 0.95, 1.1, 0.8),
       residual_variance = 1, marker_variance = 0.7,
       gamma = c(0, 1), intercept_mean = qnorm(0.15),
       intercept_sd = 1, sigma_sq_alpha = 1)
})

active_bayes_factor <- function(x) {
  variance <- x$marker_variance * x$gamma[[2L]]
  sqrt(x$residual_variance /
         (x$residual_variance + x$diagonal * variance)) *
    exp(0.5 * x$score^2 * variance /
          (x$residual_variance *
             (x$residual_variance + x$diagonal * variance)))
}
bayes_factor <- active_bayes_factor(fixture)
active_mean <- fixture$score /
  (fixture$diagonal + fixture$residual_variance /
     (fixture$marker_variance * fixture$gamma[[2L]]))
active_variance <- fixture$residual_variance /
  (fixture$diagonal + fixture$residual_variance /
     (fixture$marker_variance * fixture$gamma[[2L]]))

# Exact two-dimensional quadrature on a preregistered regular grid. Tail mass
# outside this range is checked against the proper Gaussian prior below.
alpha0 <- seq(-4.5, 2.5, length.out = 181L)
alpha1 <- seq(-3.5, 3.5, length.out = 181L)
grid <- expand.grid(intercept = alpha0, enriched = alpha1)
eta <- outer(grid$intercept, rep(1, fixture$marker_count)) +
  outer(grid$enriched, fixture$annotation[, 2L])
continuation <- pnorm(eta)
collapsed <- sweep(continuation, 2L, bayes_factor, `*`) + 1 - continuation
log_weight <- rowSums(log(collapsed)) +
  dnorm(grid$intercept, fixture$intercept_mean,
        fixture$intercept_sd, log = TRUE) +
  dnorm(grid$enriched, 0, sqrt(fixture$sigma_sq_alpha), log = TRUE)
log_weight <- log_weight - log_sum_exp(log_weight)
weight <- exp(log_weight)
conditional_pip <- sweep(continuation, 2L, bayes_factor, `*`) / collapsed
exact_pip <- colSums(conditional_pip * weight)
exact_beta_mean <- exact_pip * active_mean
exact_alpha_mean <- c(sum(grid$intercept * weight),
                      sum(grid$enriched * weight))
exact_alpha_variance <- c(sum((grid$intercept - exact_alpha_mean[[1L]])^2 * weight),
                          sum((grid$enriched - exact_alpha_mean[[2L]])^2 * weight))
exact_alpha_quantile <- rbind(
  intercept = weighted_quantile(grid$intercept, weight, c(.025, .5, .975)),
  enriched = weighted_quantile(grid$enriched, weight, c(.025, .5, .975)))

active_distribution <- numeric(fixture$marker_count + 1L)
for (row in seq_len(nrow(grid))) {
  polynomial <- c(1, rep(0, fixture$marker_count))
  for (marker in seq_len(fixture$marker_count)) {
    p <- conditional_pip[row, marker]
    next_polynomial <- polynomial * (1 - p)
    next_polynomial[2:length(next_polynomial)] <-
      next_polynomial[2:length(next_polynomial)] +
      polynomial[1:(length(polynomial) - 1L)] * p
    polynomial <- next_polynomial
  }
  active_distribution <- active_distribution + weight[[row]] * polynomial
}

draw_truncated <- function(mu, active) {
  boundary <- pnorm(-mu)
  if (active) {
    u <- runif(1L, boundary, 1)
  } else {
    u <- runif(1L, 0, boundary)
  }
  mu + qnorm(pmin(pmax(u, .Machine$double.xmin),
                  1 - .Machine$double.eps))
}

draw_alpha_given_component <- function(component, alpha) {
  current_mean <- drop(fixture$annotation %*% alpha)
  latent <- vapply(seq_len(fixture$marker_count), function(marker) {
    draw_truncated(current_mean[[marker]], component[[marker]] > 0L)
  }, numeric(1L))
  prior_precision <- diag(c(1 / fixture$intercept_sd^2,
                            1 / fixture$sigma_sq_alpha), 2L)
  precision <- crossprod(fixture$annotation) + prior_precision
  rhs <- crossprod(fixture$annotation, latent) +
    prior_precision %*% c(fixture$intercept_mean, 0)
  covariance <- solve(precision)
  mean <- drop(covariance %*% rhs)
  mean + drop(t(chol(covariance)) %*% rnorm(2L))
}

draw_component_beta <- function(alpha, marker) {
  p <- pnorm(drop(fixture$annotation[marker, , drop = FALSE] %*% alpha))
  active_probability <- p * bayes_factor[[marker]] /
    (1 - p + p * bayes_factor[[marker]])
  active <- runif(1L) < active_probability
  list(component = as.integer(active),
       beta = if (active) rnorm(1L, active_mean[[marker]],
                               sqrt(active_variance[[marker]])) else 0)
}

ordinary_step <- function(state) {
  for (marker in seq_len(fixture$marker_count)) {
    draw <- draw_component_beta(state$alpha, marker)
    state$component[[marker]] <- draw$component
    state$beta[[marker]] <- draw$beta
  }
  state$alpha <- draw_alpha_given_component(state$component, state$alpha)
  state
}

log_alpha_prior <- function(alpha) {
  dnorm(alpha[[1L]], fixture$intercept_mean, fixture$intercept_sd, log = TRUE) +
    dnorm(alpha[[2L]], 0, sqrt(fixture$sigma_sq_alpha), log = TRUE)
}

coordinated_step <- function(state, proposal_scale = 0.25, subset_size = 2L) {
  subset <- sort(sample.int(fixture$marker_count, subset_size))
  outside <- setdiff(seq_len(fixture$marker_count), subset)
  alpha_new <- state$alpha + rnorm(2L, 0, proposal_scale)
  p_old <- pnorm(drop(fixture$annotation %*% state$alpha))
  p_new <- pnorm(drop(fixture$annotation %*% alpha_new))
  z_old <- 1 - p_old + p_old * bayes_factor
  z_new <- 1 - p_new + p_new * bayes_factor
  outside_log <- sum(ifelse(state$component[outside] > 0L,
    log(p_new[outside]) - log(p_old[outside]),
    log1p(-p_new[outside]) - log1p(-p_old[outside])))
  log_mh <- log_alpha_prior(alpha_new) - log_alpha_prior(state$alpha) +
    outside_log + sum(log(z_new[subset]) - log(z_old[subset]))
  proposed_component <- state$component
  proposed_beta <- state$beta
  for (marker in subset) {
    draw <- draw_component_beta(alpha_new, marker)
    proposed_component[[marker]] <- draw$component
    proposed_beta[[marker]] <- draw$beta
  }
  alpha_jump <- sqrt(sum((alpha_new - state$alpha)^2))
  accepted <- log(runif(1L)) < min(0, log_mh)
  active_before <- sum(state$component > 0L)
  active_proposed <- sum(proposed_component > 0L)
  if (accepted) {
    state$alpha <- alpha_new
    state$component <- proposed_component
    state$beta <- proposed_beta
  }
  state$diagnostic <- c(attempt = 1, accept = accepted,
    proposed_jump = active_proposed - active_before,
    accepted_jump = if (accepted) active_proposed - active_before else 0,
    log_mh = log_mh, alpha_jump = alpha_jump)
  state
}

run_chain <- function(seed, kernel, iterations = 30000L, burnin = 5000L,
                      thinning = 5L) {
  set.seed(seed)
  state <- list(alpha = c(fixture$intercept_mean, 0),
                component = integer(fixture$marker_count),
                beta = numeric(fixture$marker_count), diagnostic = numeric())
  retained <- (iterations - burnin) %/% thinning
  alpha <- matrix(NA_real_, retained, 2L)
  component <- matrix(NA_integer_, retained, fixture$marker_count)
  beta <- matrix(NA_real_, retained, fixture$marker_count)
  diagnostic <- matrix(0, iterations, 6L)
  colnames(diagnostic) <- c("attempt", "accept", "proposed_jump",
    "accepted_jump", "log_mh", "alpha_jump")
  draw <- 0L
  for (iteration in seq_len(iterations)) {
    if (kernel == "ordinary") {
      state <- ordinary_step(state)
    } else if (kernel == "coordinated") {
      state <- coordinated_step(state)
      diagnostic[iteration, ] <- state$diagnostic
    } else {
      state <- ordinary_step(state)
      state <- coordinated_step(state)
      diagnostic[iteration, ] <- state$diagnostic
    }
    if (iteration > burnin && (iteration - burnin) %% thinning == 0L) {
      draw <- draw + 1L
      alpha[draw, ] <- state$alpha
      component[draw, ] <- state$component
      beta[draw, ] <- state$beta
    }
  }
  list(alpha = alpha, component = component, beta = beta,
       diagnostic = diagnostic)
}

seeds <- c(12011L, 12022L, 12033L, 12044L)
kernels <- c("ordinary", "coordinated", "combined")
started <- proc.time()[[3L]]
chains <- lapply(kernels, function(kernel) {
  lapply(seeds, run_chain, kernel = kernel)
})
names(chains) <- kernels
runtime <- proc.time()[[3L]] - started

summarise_kernel <- function(chain_list) {
  alpha <- do.call(rbind, lapply(chain_list, `[[`, "alpha"))
  component <- do.call(rbind, lapply(chain_list, `[[`, "component"))
  beta <- do.call(rbind, lapply(chain_list, `[[`, "beta"))
  chain_alpha_mean <- do.call(rbind, lapply(chain_list,
    function(x) colMeans(x$alpha)))
  chain_pip <- do.call(rbind, lapply(chain_list,
    function(x) colMeans(x$component > 0L)))
  diagnostic <- do.call(rbind, lapply(chain_list, `[[`, "diagnostic"))
  attempted <- sum(diagnostic[, "attempt"])
  accepted <- sum(diagnostic[, "accept"])
  list(
    alpha_mean = colMeans(alpha), alpha_variance = apply(alpha, 2L, var),
    alpha_quantile = apply(alpha, 2L, quantile, c(.025, .5, .975)),
    pip = colMeans(component > 0L), beta_mean = colMeans(beta),
    active_distribution = tabulate(rowSums(component > 0L) + 1L,
      nbins = fixture$marker_count + 1L) / nrow(component),
    alpha_chain_mcse = apply(chain_alpha_mean, 2L, sd) / sqrt(length(chain_list)),
    pip_chain_mcse = apply(chain_pip, 2L, sd) / sqrt(length(chain_list)),
    attempts = attempted, acceptances = accepted,
    acceptance_rate = if (attempted) accepted / attempted else NA_real_,
    nonzero_proposed_jumps = sum(diagnostic[, "proposed_jump"] != 0),
    nonzero_accepted_jumps = sum(diagnostic[, "accepted_jump"] != 0),
    largest_accepted_jump = if (attempted)
      max(abs(diagnostic[, "accepted_jump"])) else NA_real_)
}
summary <- lapply(chains, summarise_kernel)

# Extensive-overlap audit for a fixed two-marker exact subset. A second scale
# shrinks as 1/sqrt(m), showing the cost in alpha movement required to retain
# local acceptance as marker count grows.
scaling_counts <- c(8L, 50L, 100L, 250L, 500L, 1500L, 37991L)
scaling <- do.call(rbind, lapply(scaling_counts, function(marker_count) {
  set.seed(41000L + marker_count)
  enriched <- rep(c(0, 1), length.out = marker_count)
  annotation <- cbind(1, enriched)
  alpha <- c(fixture$intercept_mean, 0.8)
  p_old <- pnorm(drop(annotation %*% alpha))
  component <- rbinom(marker_count, 1L, p_old)
  bf <- rep(bayes_factor, length.out = marker_count)
  evaluate_scale <- function(scale) {
    log_mh <- numeric(1000L)
    alpha_jump <- numeric(1000L)
    for (attempt in seq_along(log_mh)) {
      subset <- sort(sample.int(marker_count, min(2L, marker_count)))
      outside <- setdiff(seq_len(marker_count), subset)
      proposed <- alpha + rnorm(2L, 0, scale)
      p_new <- pnorm(drop(annotation %*% proposed))
      z_old <- 1 - p_old + p_old * bf
      z_new <- 1 - p_new + p_new * bf
      outside_log <- sum(ifelse(component[outside] > 0L,
        log(p_new[outside]) - log(p_old[outside]),
        log1p(-p_new[outside]) - log1p(-p_old[outside])))
      log_mh[[attempt]] <- log_alpha_prior(proposed) - log_alpha_prior(alpha) +
        outside_log + sum(log(z_new[subset]) - log(z_old[subset]))
      alpha_jump[[attempt]] <- sqrt(sum((proposed - alpha)^2))
    }
    c(mean_log_mh = mean(log_mh), median_log_mh = median(log_mh),
      mean_acceptance = mean(pmin(1, exp(pmin(log_mh, 0)))),
      median_acceptance = median(pmin(1, exp(pmin(log_mh, 0)))),
      mean_alpha_jump = mean(alpha_jump))
  }
  fixed <- evaluate_scale(0.25)
  local <- evaluate_scale(0.25 * sqrt(8 / marker_count))
  data.frame(marker_count = marker_count, scale_policy = c("fixed", "local"),
             rbind(fixed, local), row.names = NULL)
}))

exact <- list(alpha_mean = exact_alpha_mean,
              alpha_variance = exact_alpha_variance,
              alpha_quantile = exact_alpha_quantile,
              pip = exact_pip, beta_mean = exact_beta_mean,
              active_distribution = active_distribution,
              prior_grid_tail_bound = c(
                intercept = pnorm(min(alpha0), fixture$intercept_mean,
                                  fixture$intercept_sd) +
                  pnorm(max(alpha0), fixture$intercept_mean,
                        fixture$intercept_sd, lower.tail = FALSE),
                enriched = 2 * pnorm(max(abs(alpha1)), 0,
                                     sqrt(fixture$sigma_sq_alpha),
                                     lower.tail = FALSE)))

comparison <- lapply(summary, function(x) list(
  alpha_error = x$alpha_mean - exact$alpha_mean,
  alpha_tolerance = 4 * x$alpha_chain_mcse + 0.02,
  pip_max_error = max(abs(x$pip - exact$pip)),
  pip_max_tolerance = max(4 * x$pip_chain_mcse + 0.01),
  beta_max_error = max(abs(x$beta_mean - exact$beta_mean)),
  active_distribution_l1 = sum(abs(x$active_distribution -
                                      exact$active_distribution))))

result <- list(
  schema = "sblr-bayesrc-coordinated-transition-oracle-v1",
  package_sha = system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
                        stdout = TRUE),
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  fixture = fixture, exact = exact, kernels = summary,
  comparison = comparison, scaling = scaling,
  runtime_seconds = runtime,
  decision = "AA-R5")
jsonlite::write_json(result, file.path(output_root, "oracle_result.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = 16,
                     null = "null")
saveRDS(list(result = result, chains = chains),
        file.path(output_root, "oracle_chains.rds"), compress = "xz")
print(result[c("exact", "comparison", "scaling", "runtime_seconds",
               "decision")])
