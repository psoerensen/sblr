phase17t_postburn <- function(chains, nburn) {
  stopifnot(is.matrix(chains), nburn >= 0, nburn < ncol(chains))
  chains[, seq.int(nburn + 1L, ncol(chains)), drop = FALSE]
}

phase17t_split_chains <- function(chains) {
  stopifnot(is.matrix(chains))
  half <- ncol(chains) %/% 2L
  if (half < 1L) return(matrix(numeric(), 2L * nrow(chains), 0L))
  rbind(chains[, seq_len(half), drop = FALSE],
        chains[, ncol(chains) - half + seq_len(half), drop = FALSE])
}

phase17t_rank_normalize <- function(x) {
  ranks <- rank(as.vector(x), ties.method = "average")
  stats::qnorm((ranks - 3 / 8) / (length(ranks) + 1 / 4))
}

phase17t_rhat_basic <- function(split) {
  m <- nrow(split); n <- ncol(split)
  within <- apply(split, 1L, stats::var)
  w <- mean(within)
  if (!is.finite(w) || w <= 0) return(NA_real_)
  b <- n * stats::var(rowMeans(split))
  sqrt((((n - 1) / n) * w + b / n) / w)
}

phase17t_rank_rhat <- function(chains) {
  split <- phase17t_split_chains(chains)
  z <- matrix(phase17t_rank_normalize(split), nrow(split), ncol(split))
  phase17t_rhat_basic(z)
}

phase17t_folded_rhat <- function(chains) {
  folded <- abs(chains - stats::median(as.vector(chains)))
  phase17t_rank_rhat(folded)
}

phase17t_rhat <- function(chains) {
  max(phase17t_rank_rhat(chains), phase17t_folded_rhat(chains), na.rm = TRUE)
}

phase17t_autocovariance <- function(x) {
  n <- length(x)
  variance <- stats::var(x)
  if (variance == 0) return(rep(0, n))
  size <- 2L * stats::nextn(n)
  centered <- c(x - mean(x), rep(0, size - n))
  transformed <- fft(centered)
  autocov <- Re(fft(Mod(transformed)^2, inverse = TRUE))[seq_len(n)]
  autocov / autocov[1L] * variance * (n - 1) / n
}

phase17t_ess_split <- function(split) {
  n <- ncol(split); m <- nrow(split)
  if (m < 2L || n < 3L) return(NA_real_)
  x <- t(split)
  acov <- vapply(seq_len(m), function(i) phase17t_autocovariance(x[, i]),
                  numeric(n))
  acov_means <- rowMeans(acov)
  w <- acov_means[1L] * n / (n - 1)
  var_plus <- w * (n - 1) / n + stats::var(colMeans(x))
  if (!is.finite(var_plus) || var_plus <= 0) return(NA_real_)
  rho <- numeric(n); at <- 0L; even <- 1
  rho[1L] <- even
  odd <- 1 - (w - acov_means[2L]) / var_plus
  rho[2L] <- odd
  while (at < n - 5L && is.finite(even + odd) && even + odd > 0) {
    at <- at + 2L
    even <- 1 - (w - acov_means[at + 1L]) / var_plus
    odd <- 1 - (w - acov_means[at + 2L]) / var_plus
    if (is.finite(even + odd) && even + odd >= 0) {
      rho[at + 1L] <- even
      rho[at + 2L] <- odd
    }
  }
  max_at <- at
  if (even > 0) rho[max_at + 1L] <- even
  at <- 0L
  while (at <= max_at - 4L) {
    at <- at + 2L
    if (rho[at + 1L] + rho[at + 2L] > rho[at - 1L] + rho[at]) {
      rho[at + 1L] <- (rho[at - 1L] + rho[at]) / 2
      rho[at + 2L] <- rho[at + 1L]
    }
  }
  retained <- if (max_at == 0L) rho[1L] else sum(rho[seq_len(max_at)])
  tau <- -1 + 2 * retained + rho[max_at + 1L]
  bound <- 1 / log10(m * n)
  if (tau < bound) tau <- bound
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  m * n / tau
}

phase17t_ess_bulk <- function(chains) {
  split <- phase17t_split_chains(chains)
  z <- matrix(phase17t_rank_normalize(split), nrow(split), ncol(split))
  phase17t_ess_split(z)
}

phase17t_ess_tail <- function(chains) {
  pooled <- as.vector(chains)
  limits <- stats::quantile(pooled, c(0.05, 0.95), names = FALSE, type = 7)
  values <- vapply(limits, function(q) {
    phase17t_ess_split(phase17t_split_chains(1L * (chains <= q)))
  }, numeric(1))
  if (anyNA(values)) NA_real_ else min(values)
}

phase17t_ess_mean <- function(chains) {
  phase17t_ess_split(phase17t_split_chains(chains))
}

phase17t_mcse_mean <- function(chains) {
  posterior_sd <- stats::sd(as.vector(chains))
  ess <- phase17t_ess_mean(chains)
  if (!is.finite(posterior_sd) || posterior_sd <= 0 || !is.finite(ess)) {
    return(c(ess_mean = NA_real_, posterior_sd = posterior_sd,
             mcse_mean = NA_real_, mcse_mean_over_sd = NA_real_))
  }
  mcse <- posterior_sd / sqrt(ess)
  c(ess_mean = ess, posterior_sd = posterior_sd, mcse_mean = mcse,
    mcse_mean_over_sd = mcse / posterior_sd)
}

phase17t_status <- function(chains, updated = TRUE) {
  if (!updated) return("not_updated")
  if (any(!is.finite(chains))) return("nonfinite")
  if (nrow(chains) < 2L) return("unavailable_single_chain")
  if (ncol(chains) < 4L) return("insufficient_draws")
  if (length(unique(as.vector(chains))) == 1L) return("constant")
  chain_constant <- apply(chains, 1L, function(x) length(unique(x)) == 1L)
  if (any(chain_constant)) return("constant_chain_mismatch")
  if (nrow(chains) < 4L) "computed_fewer_than_four_chains" else "computed"
}

phase17t_flags <- function(rhat, ess_bulk, ess_tail, mcse_relative, nchains,
                           rhat_threshold = 1.01,
                           ess_per_chain_threshold = 100,
                           mcse_threshold = 0.05) {
  c(rhat_flag = is.finite(rhat) && rhat > rhat_threshold,
    ess_bulk_flag = is.finite(ess_bulk) &&
      ess_bulk < ess_per_chain_threshold * nchains,
    ess_tail_flag = is.finite(ess_tail) &&
      ess_tail < ess_per_chain_threshold * nchains,
    mcse_flag = is.finite(mcse_relative) && mcse_relative > mcse_threshold)
}

phase17t_overview <- function(summary) {
  computed <- summary$status %in% c("computed", "computed_fewer_than_four_chains")
  flagged <- rowSums(summary[c("rhat_flag", "ess_bulk_flag",
                               "ess_tail_flag", "mcse_flag")]) > 0L
  list(overall_status = if (!any(computed)) "unavailable" else if (any(flagged))
         "warning" else if (all(computed)) "ok" else "partial",
       n_computed = sum(computed), n_unavailable = sum(!computed),
       n_flagged = sum(flagged), max_rhat = max(summary$rhat, na.rm = TRUE),
       min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
       min_ess_tail = min(summary$ess_tail, na.rm = TRUE),
       max_mcse_mean_over_sd = max(summary$mcse_mean_over_sd, na.rm = TRUE),
       fewer_than_four_chains = summary$nchains[1L] < 4L)
}

phase17t_memory <- function(nchains, nit, nt, nmodels, selected_markers) {
  q <- nt * (nt + 1) / 2
  c(tier1_trace_bytes = 8 * nchains * nit * 3 * nt,
    covariance_trace_bytes = 8 * nchains * nit * 3 * q,
    probability_trace_bytes = 8 * nchains * nit * 2,
    full_pi_trace_bytes = 8 * nchains * nit * nmodels,
    selected_b_trace_bytes = 8 * nchains * nit * selected_markers * nt,
    selected_d_trace_bytes = 4 * nchains * nit * selected_markers * nt,
    per_quantity_workspace_bytes = 8 * nchains * nit)
}

phase17t_scope_contract <- function() {
  list(tier1 = c("B_diag", "G_diag", "E_diag"),
       tier2_requires_new_traces = c("B_lower", "G_lower", "E_lower",
                                    "pi_null", "pi_active"),
       tier3_opt_in = c("selected_marker_b", "selected_marker_d"),
       all_markers_default = FALSE, full_pi_default = FALSE,
       keep_chains_required = FALSE, diagnostic_thinning = FALSE,
       trace_retention = "independent_bundle")
}

phase17t_fixtures <- function() {
  half <- seq(-2, 2, length.out = 20)
  base <- t(vapply(0:3, function(i) {
    perm <- c(half[(seq_along(half) + 5L * i - 1L) %% 20L + 1L])
    c(perm, rev(perm))
  }, numeric(40)))
  list(well_mixed = base,
       shifted = base + c(0, 0, 1.5, 1.5),
       scales = base * c(1, 1, 2.5, 2.5),
       drift = base + rep(seq(0, 2, length.out = 40), each = 4),
       positive_ar = t(vapply(1:4, function(i) cumsum(sin((1:40 + i) / 7)),
                              numeric(40))),
       negative_ar = t(vapply(1:4, function(i) (-1)^(1:40) + i / 100,
                              numeric(40))),
       poor_tail = rbind(base[1:3, ], c(base[4, 1:35], rep(8, 5))),
       tied = matrix(rep(c(0, 0, 1, 1), 40), 4),
       binary = matrix(rep(c(0, 1), 80), 4),
       constant = matrix(2, 4, 40),
       one_constant = rbind(rep(0, 40), base[2:4, ]),
       nonfinite = { x <- base; x[1, 1] <- NA_real_; x },
       one_chain = base[1, , drop = FALSE],
       two_chains = base[1:2, , drop = FALSE],
       four_chains = base,
       odd = base[, 1:39], short = base[, 1:3])
}
