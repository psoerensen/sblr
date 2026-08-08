.sbs_em_ref_gh <- function(order = 15L) {
 index <- seq_len(order - 1L)
 jacobi <- diag(0, order)
 jacobi[cbind(index, index + 1L)] <- sqrt(index / 2)
 jacobi[cbind(index + 1L, index)] <- sqrt(index / 2)
 decomposition <- eigen(jacobi, symmetric = TRUE)
 ordering <- order(decomposition$values)
 list(
  node = decomposition$values[ordering],
  weight = sqrt(pi) * decomposition$vectors[1L, ordering]^2
 )
}

.sbs_em_ref_stick_marginal <- function(
  A, eligible, success, selected, tau2, intercept_mean, intercept_sd,
  order = 15L, floor = 1e-12
) {
 columns <- c(1L, which(selected) + 1L)
 design <- A[, columns, drop = FALSE]
 prior_mean <- c(intercept_mean, rep(0, sum(selected)))
 prior_sd <- c(intercept_sd, rep(sqrt(tau2), sum(selected)))
 rule <- .sbs_em_ref_gh(order)
 grid <- expand.grid(rep(list(seq_len(order)), length(columns)))
 log_term <- numeric(nrow(grid))
 coefficient_sum <- numeric(length(columns))
 for (row in seq_len(nrow(grid))) {
  location <- as.integer(grid[row, ])
  parameter <- prior_mean + sqrt(2) * prior_sd * rule$node[location]
  probability <- pmin(pmax(
   stats::pnorm(drop(design %*% parameter)), floor
  ), 1 - floor)
  log_likelihood <- sum(
   success * log(probability) +
    (eligible - success) * log1p(-probability)
  )
  log_weight <- sum(log(rule$weight[location] / sqrt(pi)))
  log_term[row] <- log_weight + log_likelihood
 }
 log_marginal <- .sbayesrc_s_em_log_sum_exp(log_term)
 normalized <- exp(log_term - log_marginal)
 for (row in seq_len(nrow(grid))) {
  location <- as.integer(grid[row, ])
  parameter <- prior_mean + sqrt(2) * prior_sd * rule$node[location]
  coefficient_sum <- coefficient_sum + normalized[row] * parameter
 }
 full_mean <- numeric(ncol(A))
 full_mean[columns] <- coefficient_sum
 list(log_marginal = log_marginal, mean = full_mean)
}

.sbs_em_ref_exact <- function(
  A, responsibility, pi_a, tau2, intercept_prior_resolved,
  order = 15L
) {
 soft <- .sbayesrc_mcem_soft_stick_information(responsibility)
 stick_count <- ncol(responsibility) - 1L
 annotation_count <- ncol(A) - 1L
 intercept <- .sbayesrc_mcem_intercept_prior(
  intercept_prior_resolved, stick_count
 )
 state <- as.matrix(expand.grid(rep(list(0:1), annotation_count)))
 colnames(state) <- colnames(A)[-1L]
 log_weight <- numeric(nrow(state))
 alpha_mean <- array(0, c(nrow(A), 0L, 0L))
 model_alpha <- array(0, c(ncol(A), stick_count, nrow(state)))
 for (model in seq_len(nrow(state))) {
  delta <- state[model, ]
  selected <- delta == 1L
  log_weight[model] <- sum(delta) * log(pi_a) +
   (annotation_count - sum(delta)) * log1p(-pi_a)
  for (stick in seq_len(stick_count)) {
   marginal <- .sbs_em_ref_stick_marginal(
    A, soft$eligible[, stick], soft$success[, stick], selected,
    tau2[stick], intercept$mean[stick],
    1 / sqrt(intercept$precision[stick]), order
   )
   log_weight[model] <- log_weight[model] + marginal$log_marginal
   model_alpha[, stick, model] <- marginal$mean
  }
 }
 probability <- exp(log_weight - .sbayesrc_s_em_log_sum_exp(log_weight))
 list(
  state = state,
  model_probability = probability,
  annotation_pip_eb = drop(crossprod(probability, state)),
  alpha_model_average = apply(
   model_alpha * rep(probability, each = ncol(A) * stick_count),
   c(1L, 2L), sum
  ),
  log_weight = log_weight
 )
}

.sbs_em_ref_fixture <- function(correlated = FALSE, empty_last_stick = FALSE) {
 set.seed(55103L)
 marker_count <- 72L
 signal <- rep(c(-1, 1), length.out = marker_count)
 null <- rep(c(-1, -1, 1, 1), length.out = marker_count)
 if (isTRUE(correlated)) null <- 0.92 * signal + 0.08 * null
 A <- cbind(intercept = 1, signal = signal, candidate = null)
 alpha <- matrix(c(
  -0.35, 0.80, 0.00,
   0.10, 0.55, 0.00
 ), nrow = 3L, ncol = 2L)
 responsibility <- .sbayesrc_mcem_component_prior(A, alpha)
 if (isTRUE(empty_last_stick)) {
  responsibility[,] <- 0
  responsibility[, 1L] <- 1
 }
 intercept_prior <- rbind(
  type = c(0, 0), mean = c(-0.25, 0.15), precision = c(1, 1)
 )
 list(
  A = A,
  responsibility = responsibility,
  alpha = alpha,
  alpha_start = matrix(0, ncol(A), 2L),
  delta = c(1L, 0L),
  pi_a = 0.30,
  tau2 = c(0.8, 0.8),
  intercept_prior = intercept_prior
 )
}

.sbs_em_run_csr <- function(fixture, delta_start = 0L,
                            seed = 8810L, sweeps = 350L,
                            burn = 100L, outer = 8L) {
 sblr:::.stblr_mcem_sbayesrc_s_csr(
  fixture$wy, fixture$ww, fixture$yy, fixture$b, fixture$component,
  fixture$r, fixture$prefix, fixture$B, fixture$E, fixture$ssb,
  fixture$sse, fixture$A, fixture$gamma, fixture$alpha_start,
  delta_init = delta_start, pi_a = 0.30,
  tau2 = fixture$sigmaSqAlpha,
  intercept_prior_resolved = fixture$intercept_prior$native,
  n = fixture$sample_size, inner_sweeps = sweeps, inner_burn = burn,
  selection_sweeps = 800L, selection_burn = 200L,
  final_sweeps = sweeps, final_burn = burn,
  max_outer = outer, seed = seed, ncores = 1L
 )
}

.sbs_em_run_block <- function(fixture, delta_start = 0L,
                              alpha_start = matrix(0, 2L, 2L),
                              seed = 8820L, updateB = FALSE,
                              updateE = FALSE, sweeps = 350L,
                              burn = 120L, outer = 8L) {
 sblr:::.stblr_mcem_sbayesrc_s_block_eigen(
  stats = fixture$stats, Glist = fixture$Glist,
  annotation = fixture$A, block_start = 1L,
  B = fixture$B, E = fixture$E,
  ssb_prior = fixture$ssb_prior, sse_prior = fixture$sse_prior,
  gamma = fixture$gamma, alpha_init = alpha_start,
  delta_init = delta_start, pi_a = 0.30, tau2 = c(1, 1),
  intercept_prior_resolved = fixture$intercept_prior,
  representation = "low_rank", eigen_prop = 0.999999,
  residual_policy = "gctb_block",
  block_ve_mode = if (isTRUE(updateE)) "allMixVe" else "fixVe",
  updateB = updateB, updateE = updateE,
  inner_sweeps = sweeps, inner_burn = burn,
  selection_sweeps = 800L, selection_burn = 200L,
  final_sweeps = sweeps, final_burn = burn,
  max_outer = outer, seed = seed, ncores = 1L
 )
}

.sbs_em_run_block_as_csr <- function(
  fixture, delta_start = 0L, alpha_start = matrix(0, 2L, 2L),
  seed = 8830L, sweeps = 350L, burn = 120L, outer = 8L
) {
 marker_count <- nrow(fixture$A)
 sblr:::.stblr_mcem_sbayesrc_s_csr(
  fixture$stats$wy, fixture$stats$ww, as.numeric(fixture$stats$yy[[1L]]),
  list(rep(0, marker_count)), list(rep(0L, marker_count)),
  list(as.numeric(fixture$stats$wy[[1L]])), fixture$prefix,
  fixture$B, fixture$E, fixture$ssb_prior, fixture$sse_prior,
  fixture$A, fixture$gamma, alpha_start, delta_init = delta_start,
  pi_a = 0.30, tau2 = c(1, 1),
  intercept_prior_resolved = fixture$intercept_prior,
  n = fixture$stats$n, inner_sweeps = sweeps, inner_burn = burn,
  selection_sweeps = 800L, selection_burn = 200L,
  final_sweeps = sweeps, final_burn = burn,
  max_outer = outer, seed = seed, ncores = 1L
 )
}
