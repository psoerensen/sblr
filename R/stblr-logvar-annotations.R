.stblr_preprocess_logvar_annotations <- function(annotations, marker_ids) {
  if (!is.matrix(annotations) && !is.data.frame(annotations)) {
    stop("log_variance annotations must be a matrix or data frame.", call. = FALSE)
  }
  original <- annotations
  A <- as.matrix(annotations)
  if (!is.numeric(A)) {
    stop("log_variance annotation columns must be numeric or logical.", call. = FALSE)
  }
  storage.mode(A) <- "double"
  if (nrow(A) != length(marker_ids) || ncol(A) < 1L) {
    stop("log_variance annotations must have exactly m rows and at least one column.",
         call. = FALSE)
  }
  annotation_rows <- rownames(original)
  default_rows <- !is.null(annotation_rows) &&
    identical(annotation_rows, as.character(seq_len(nrow(A))))
  if (!is.null(annotation_rows) && !default_rows &&
      !identical(as.character(annotation_rows), as.character(marker_ids))) {
    stop("log_variance annotation row names must exactly match final marker order.",
         call. = FALSE)
  }
  if (anyNA(A) || any(!is.finite(A))) {
    stop("log_variance annotations must not contain missing or non-finite values.",
         call. = FALSE)
  }
  annotation_names <- colnames(A)
  if (is.null(annotation_names)) annotation_names <- paste0("annotation_", seq_len(ncol(A)))
  if (any(!nzchar(annotation_names)) || anyDuplicated(annotation_names)) {
    stop("log_variance annotation column names must be non-empty and unique.",
         call. = FALSE)
  }
  X <- matrix(0, nrow(A), ncol(A), dimnames = list(marker_ids, annotation_names))
  means <- sds <- numeric(ncol(A))
  types <- transforms <- character(ncol(A))
  for (column in seq_len(ncol(A))) {
    value <- A[, column]
    if (all(value == 1)) {
      stop("log_variance annotations must not contain an explicit all-ones intercept.",
           call. = FALSE)
    }
    unique_value <- unique(value)
    if (length(unique_value) < 2L) {
      stop("log_variance annotation columns must not be constant.", call. = FALSE)
    }
    binary <- all(unique_value %in% c(0, 1))
    means[column] <- mean(value)
    if (binary) {
      X[, column] <- value - means[column]
      sds[column] <- 1
      types[column] <- "binary"
      transforms[column] <- "center_only"
    } else {
      scale_sd <- stats::sd(value)
      if (!is.finite(scale_sd) || scale_sd <= 0) {
        stop("continuous log_variance annotation columns must have positive finite SD.",
             call. = FALSE)
      }
      X[, column] <- (value - means[column]) / scale_sd
      sds[column] <- scale_sd
      types[column] <- "continuous"
      transforms[column] <- "center_and_sd"
    }
  }
  if (qr(X, tol = sqrt(.Machine$double.eps))$rank < ncol(X)) {
    stop("log_variance annotation design is rank deficient or contains duplicate columns.",
         call. = FALSE)
  }
  list(
    X = X,
    transform = data.frame(
      annotation = annotation_names, type = types, mean = means, sd = sds,
      transform = transforms, stringsAsFactors = FALSE),
    marker_alignment = if (is.null(annotation_rows) || default_rows)
      "positional_exact_length" else "exact_marker_ids"
  )
}

.stblr_logvar_theta_init <- function(theta_init, p, nt, annotation_names,
                                     trait_names) {
  if (is.null(theta_init)) theta_init <- matrix(0, p, nt)
  if (is.vector(theta_init) && nt == 1L) theta_init <- matrix(theta_init, p, 1L)
  theta_init <- as.matrix(theta_init)
  if (!identical(dim(theta_init), c(p, nt)) || anyNA(theta_init) ||
      any(!is.finite(theta_init))) {
    stop("theta_init must be NULL or a finite annotation-by-trait matrix.",
         call. = FALSE)
  }
  dimnames(theta_init) <- list(annotation_names, trait_names)
  theta_init
}

.stblr_attach_logvar_output <- function(fit, raw, transform, marker_ids,
                                        trait_names, method, theta_prior_sd,
                                        updateTheta) {
  theta <- as.matrix(raw$annotation$theta)
  rownames(theta) <- transform$annotation
  colnames(theta) <- trait_names
  ratio <- as.matrix(raw$annotation$variance_ratio)
  rownames(ratio) <- transform$annotation
  colnames(ratio) <- trait_names
  q <- as.matrix(raw$annotation$marker_prior_scale)
  rownames(q) <- marker_ids
  colnames(q) <- trait_names
  fit$theta <- theta
  fit$annotation_variance_ratio <- ratio
  fit$marker_prior_scale <- q
  fit$annotation_transform <- transform
  trace <- raw$annotation$theta_trace
  trace_dim <- dim(trace)
  burn <- as.integer(fit$input$nburn %||% 0L)
  draws <- seq.int(burn + 1L, trace_dim[1L])
  nchains <- trace_dim[3L] %/% length(trait_names)
  summary_rows <- vector("list", nrow(theta) * ncol(theta))
  chain_mean <- array(
    NA_real_, c(nrow(theta), length(trait_names), nchains),
    dimnames = list(transform$annotation, trait_names, paste0("chain", seq_len(nchains))))
  row <- 0L
  for (trait in seq_along(trait_names)) for (annotation in seq_len(nrow(theta))) {
    tasks <- (trait - 1L) * nchains + seq_len(nchains)
    value <- matrix(trace[draws, annotation, tasks, drop = FALSE],
                    nrow = length(draws), ncol = nchains)
    chain_mean[annotation, trait, ] <- colMeans(value)
    diagnostic <- .blr_convergence_scalar(value, updated = updateTheta)
    row <- row + 1L
    summary_rows[[row]] <- data.frame(
      annotation = transform$annotation[annotation], trait = trait_names[trait],
      mean = mean(value), sd = stats::sd(as.numeric(value)),
      median = stats::median(value),
      lower = stats::quantile(value, 0.025, names = FALSE),
      upper = stats::quantile(value, 0.975, names = FALSE),
      Rhat = diagnostic$rhat, bulk_ESS = diagnostic$ess_bulk,
      tail_ESS = diagnostic$ess_tail, MCSE = diagnostic$mcse_mean,
      stringsAsFactors = FALSE)
  }
  fit$theta_summary <- do.call(rbind, summary_rows)
  fit$theta_chain_mean <- chain_mean
  if (isTRUE(fit$input$keep_chains)) fit$theta_trace <- trace
  fit$diagnostics$logvar <- raw$diagnostics$logvar
  fit$input$method <- method
  fit$input$model <- paste0(method, "_logvar")
  fit$input$backend <- paste0("csr_logvar_", sub("sbayes", "bayes", method))
  fit$input$data_level <- "summary"
  fit$input$annotations <- TRUE
  fit$input$annotation_model <- "log_variance"
  fit$input$annotation_policy <- "log_variance"
  fit$input$probability_policy <- "global"
  fit$input$prior_kernel <- sub("sbayes", "bayes", method)
  fit$input$effect_scale <- "annotation_log_variance"
  fit$input$effect_scale_policy <- "annotation_log_variance"
  fit$input$theta_prior_sd <- theta_prior_sd
  fit$input$updateTheta <- updateTheta
  fit
}

.stblr_csr_logvar_bayesc <- function(
    stats, ld_prefix, annotation_info, theta_prior_sd, theta_init, updateTheta,
    h2 = 0.3, adjE = 0.9, nit = 1000, nburn = 100, nthin = 1,
    ncores = 1L, seed = 1L, nchains = 1L, keep_chains = FALSE,
    chain_seeds = NULL, updateB = TRUE, updateE = TRUE, updatePi = TRUE,
    updateLDswap = FALSE, ld_swap_prob = 0.05, ld_swap_r2 = 0.8,
    ld_swap_max_friends = 50L, ld_swap_moves = 1L, nub = 4, nue = 4,
    pi_init = NULL, pi_vb_init = NULL, pi_prior_mean = NULL,
    pi_prior_strength = NULL, pi_prior_a = NULL, pi_prior_b = NULL,
    b_init = NULL, d_init = NULL, use_d_init = FALSE, r_init = NULL,
    use_r_init = FALSE, rebuild_r_before_updateE = FALSE,
    .convergence_spec = NULL) {
  chain <- .stblr_validate_annotation_chain_args(nchains, keep_chains, chain_seeds)
  dims <- .stblr_get_nt_m_names(stats)
  .stblr_validate_stats(stats, dims$nt, dims$m)
  arch <- .stblr_resolve_architecture(
    pi_marker = pi_init %||% 0.001, pi_init = pi_init,
    pi_vb_init = pi_vb_init, pi_prior_mean = pi_prior_mean,
    pi_prior_strength = pi_prior_strength, pi_prior_a = pi_prior_a,
    pi_prior_b = pi_prior_b)
  pri <- .stblr_make_csr_variance_priors(
    stats, dims$n, dims$m, dims$nt, h2, nub, nue, arch$pi_vb_init,
    arch$pi_prior_mean, dims$trait_names)
  state <- .stblr_init_marker_state(dims$nt, dims$m, b_init, d_init)
  residual <- .stblr_init_r_state(
    stats, dims$nt, dims$m, use_r_init, r_init)
  theta_init <- .stblr_logvar_theta_init(
    theta_init, ncol(annotation_info$X), dims$nt,
    colnames(annotation_info$X), dims$trait_names)
  raw <- stblr_cpg_omp_csr_logvar_bayesc(
    stats$wy, stats$ww, stats$yy, state$b_init, state$d_init, use_d_init,
    residual, use_r_init, rebuild_r_before_updateE, ld_prefix, pri$B, pri$E,
    pri$ssb_prior_list, pri$sse_prior_list, arch$pi, nub, nue, updateB,
    updateE, updatePi, adjE, rep(as.integer(dims$n), dims$nt),
    as.integer(nit), as.integer(nburn), as.integer(nthin), arch$pi_prior_a,
    arch$pi_prior_b, as.integer(ncores), as.integer(seed), chain$nchains,
    chain$keep_chains, chain$chain_seeds, updateLDswap, ld_swap_prob,
    ld_swap_r2, as.integer(ld_swap_max_friends), as.integer(ld_swap_moves),
    annotation_info$X, theta_init, theta_prior_sd, updateTheta,
    .convergence_spec$markers %||% integer(),
    isTRUE(.convergence_spec$b), isTRUE(.convergence_spec$d))
  fit <- .as_stblr_fit(raw, dims$trait_names, dims$variable_names)
  fit$input <- c(list(
    n = dims$n, m = dims$m, nt = dims$nt, h2 = h2, nub = nub, nue = nue,
    nit = nit, nburn = nburn, nthin = nthin, ncores = ncores, seed = seed,
    nchains = chain$nchains, keep_chains = chain$keep_chains,
    chain_seeds = if (length(chain$chain_seeds)) chain$chain_seeds else NULL,
    updateB = updateB, updateE = updateE, updatePi = updatePi,
    updateLDswap = updateLDswap, ld_prefix = ld_prefix), arch)
  .stblr_attach_logvar_output(
    fit, raw, annotation_info$transform, dims$variable_names,
    dims$trait_names, "sbayesc", theta_prior_sd, updateTheta)
}

.stblr_csr_logvar_bayesr <- function(
    stats, ld_prefix, annotation_info, theta_prior_sd, theta_init, updateTheta,
    h2 = 0.3, adjE = 0.9, nit = 1000, nburn = 100, nthin = 1,
    ncores = 1L, seed = 1L, nchains = 1L, keep_chains = FALSE,
    chain_seeds = NULL, updateB = TRUE, updateE = TRUE, updatePi = TRUE,
    updateLDswap = FALSE, ld_swap_prob = 0.05, ld_swap_r2 = 0.8,
    ld_swap_max_friends = 50L, ld_swap_moves = 1L, nub = 4, nue = 4,
    pi = NULL, mixture_var = c(0, 0.01, 0.1, 1), alpha = NULL,
    comp_init = NULL, use_comp_init = FALSE, r_init = NULL,
    use_r_init = FALSE, rebuild_r_before_updateE = FALSE,
    updateE_start = 0L, updateE_every = 1L, .convergence_spec = NULL) {
  chain <- .stblr_validate_annotation_chain_args(nchains, keep_chains, chain_seeds)
  dims <- .stblr_get_nt_m_names(stats)
  .stblr_validate_stats(stats, dims$nt, dims$m)
  if (!is.numeric(mixture_var) || length(mixture_var) < 2L ||
      mixture_var[1L] != 0 || any(!is.finite(mixture_var)) ||
      any(mixture_var[-1L] <= 0))
    stop("mixture_var must start with 0 and have positive non-null components.")
  if (is.null(pi)) {
    active <- 0.001
    pi <- c(1 - active, rep(active / (length(mixture_var) - 1L),
                            length(mixture_var) - 1L))
  }
  if (length(pi) != length(mixture_var) || any(!is.finite(pi)) ||
      any(pi < 0) || sum(pi) <= 0) stop("pi must match mixture_var.")
  pi <- pi / sum(pi)
  if (is.null(alpha)) alpha <- pmax(pi * 5e5, .Machine$double.eps)
  if (length(alpha) != length(pi) || any(!is.finite(alpha)) || any(alpha <= 0))
    stop("alpha must be positive and match mixture_var.")
  vy <- as.numeric(stats$yy) / (as.numeric(dims$n) - 1)
  pri <- .make_stblr_bayesr_priors(
    vy, dims$m, h2, nub, nue, pi, mixture_var, dims$trait_names, alpha)
  if (is.null(comp_init)) comp_init <- lapply(seq_len(dims$nt), function(i) rep(0, dims$m))
  if (is.null(r_init)) r_init <- stats$wy
  theta_init <- .stblr_logvar_theta_init(
    theta_init, ncol(annotation_info$X), dims$nt,
    colnames(annotation_info$X), dims$trait_names)
  raw <- stblr_cpg_omp_csr_logvar_bayesr(
    stats$wy, stats$ww, stats$yy,
    lapply(seq_len(dims$nt), function(i) rep(0, dims$m)), comp_init,
    use_comp_init, r_init, use_r_init, rebuild_r_before_updateE, ld_prefix,
    pri$B, pri$E, pri$ssb_prior_list, pri$sse_prior_list, pi, mixture_var,
    alpha, nub, nue, updateB, updateE, updatePi, adjE, as.integer(dims$n),
    as.integer(nit), as.integer(nburn), as.integer(nthin), as.integer(ncores),
    as.integer(seed), chain$nchains, chain$keep_chains, chain$chain_seeds,
    as.integer(updateE_start), as.integer(updateE_every), updateLDswap,
    ld_swap_prob, ld_swap_r2, as.integer(ld_swap_max_friends),
    as.integer(ld_swap_moves), annotation_info$X, theta_init, theta_prior_sd,
    updateTheta, .convergence_spec$markers %||% integer(),
    isTRUE(.convergence_spec$probability), isTRUE(.convergence_spec$b),
    isTRUE(.convergence_spec$d), isTRUE(.convergence_spec$component))
  fit <- .as_stblr_fit(raw, dims$trait_names, dims$variable_names)
  fit$input <- list(
    n = dims$n, m = dims$m, nt = dims$nt, h2 = h2, nub = nub, nue = nue,
    nit = nit, nburn = nburn, nthin = nthin, ncores = ncores, seed = seed,
    nchains = chain$nchains, keep_chains = chain$keep_chains,
    chain_seeds = if (length(chain$chain_seeds)) chain$chain_seeds else NULL,
    updateB = updateB, updateE = updateE, updatePi = updatePi,
    mixture_var = mixture_var, alpha = alpha, updateLDswap = updateLDswap,
    ld_prefix = ld_prefix)
  .stblr_attach_logvar_output(
    fit, raw, annotation_info$transform, dims$variable_names,
    dims$trait_names, "sbayesr", theta_prior_sd, updateTheta)
}
