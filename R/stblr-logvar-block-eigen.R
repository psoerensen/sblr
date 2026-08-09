.stblr_logvar_block_metadata <- function(
    fit, raw, bed, representation, eigen_policy, eigen_prop, eigen_filter,
    eigen_tau, eigen_eta, low_rank_residual_rebuild_every) {
  fit$input$ld_backend <- "block_eigen"
  fit$input$operator <- "block_eigen"
  fit$input$operator_contract <- if (identical(representation, "low_rank"))
    "block_low_rank_v1" else "block_dense_reconstructed_v1"
  fit$input$operator_representation <- representation
  fit$input$operator_scale_contract <- if (identical(representation, "low_rank"))
    "general_cross_product" else "reconstructed_cross_product"
  fit$input$eigen_policy <- eigen_policy
  fit$input$eigen_prop <- if (identical(representation, "low_rank"))
    eigen_prop else NULL
  fit$input$eigen_filter <- eigen_filter
  fit$input$eigen_tau <- eigen_tau
  fit$input$eigen_eta <- eigen_eta
  fit$input$eigen_blocks <- bed$block_start
  fit$input$eigen_diagnostics <- raw$diagnostics$block_eigen
  fit$input$low_rank_residual_rebuild_every <-
    low_rank_residual_rebuild_every
  fit$diagnostics$block_eigen <- raw$diagnostics$block_eigen
  fit
}

.stblr_block_logvar_bayesc <- function(
    stats, Glist, block_start, annotation_info, theta_prior_sd, theta_init,
    updateTheta, representation = "low_rank",
    eigen_policy = "cumulative_positive_mass", eigen_prop = 0.995,
    eigen_filter = "hard_truncate", eigen_tau = 0.01, eigen_eta = 0,
    low_rank_residual_rebuild_every = 100L,
    h2 = 0.3, adjE = 0.9, nit = 1000, nburn = 100, nthin = 1,
    ncores = 1L, seed = 1L, nchains = 1L, keep_chains = FALSE,
    chain_seeds = NULL, updateB = TRUE, updateE = TRUE, updatePi = TRUE,
    nub = 4, nue = 4, pi_init = NULL, pi_vb_init = NULL,
    pi_prior_mean = NULL, pi_prior_strength = NULL, pi_prior_a = NULL,
    pi_prior_b = NULL, b_init = NULL, d_init = NULL, use_d_init = FALSE,
    r_init = NULL, use_r_init = FALSE, rebuild_r_before_updateE = FALSE,
    .convergence_spec = NULL) {
  chain <- .stblr_validate_annotation_chain_args(
    nchains, keep_chains, chain_seeds)
  dims <- .stblr_get_nt_m_names(stats)
  .stblr_validate_stats(stats, dims$nt, dims$m)
  bed <- .stblr_csr_block_eigen_inputs(stats, Glist, block_start)
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
  raw <- stblr_cpg_omp_block_eigen_logvar_bayesc(
    stats$wy, stats$ww, stats$yy, state$b_init, state$d_init, use_d_init,
    residual, use_r_init, rebuild_r_before_updateE, pri$B, pri$E,
    pri$ssb_prior_list, pri$sse_prior_list, arch$pi, nub, nue, updateB,
    updateE, updatePi, adjE, rep(as.integer(dims$n), dims$nt),
    as.integer(nit), as.integer(nburn), as.integer(nthin), arch$pi_prior_a,
    arch$pi_prior_b, as.integer(ncores), as.integer(seed), chain$nchains,
    chain$keep_chains, chain$chain_seeds,
    .convergence_spec$markers %||% integer(),
    isTRUE(.convergence_spec$b), isTRUE(.convergence_spec$d),
    bed$bed_files, bed$n_bed, bed$cls, bed$rows, bed$af, bed$block_start,
    eigen_filter, eigen_tau, eigen_eta, representation, eigen_prop,
    as.integer(low_rank_residual_rebuild_every), annotation_info$X,
    theta_init, theta_prior_sd, updateTheta)
  fit <- .as_stblr_fit(raw, dims$trait_names, dims$variable_names)
  fit$input <- c(list(
    n = dims$n, m = dims$m, nt = dims$nt, h2 = h2, nub = nub, nue = nue,
    nit = nit, nburn = nburn, nthin = nthin, ncores = ncores, seed = seed,
    nchains = chain$nchains, keep_chains = chain$keep_chains,
    chain_seeds = if (length(chain$chain_seeds)) chain$chain_seeds else NULL,
    updateB = updateB, updateE = updateE, updatePi = updatePi,
    updateLDswap = FALSE), arch)
  fit <- .stblr_attach_logvar_output(
    fit, raw, annotation_info$transform, dims$variable_names,
    dims$trait_names, "sbayesc", theta_prior_sd, updateTheta)
  fit$input$backend <- "block_eigen_logvar_bayesc"
  .stblr_logvar_block_metadata(
    fit, raw, bed, representation, eigen_policy, eigen_prop, eigen_filter,
    eigen_tau, eigen_eta, low_rank_residual_rebuild_every)
}

.stblr_block_logvar_bayesr <- function(
    stats, Glist, block_start, annotation_info, theta_prior_sd, theta_init,
    updateTheta, representation = "low_rank",
    eigen_policy = "cumulative_positive_mass", eigen_prop = 0.995,
    eigen_filter = "hard_truncate", eigen_tau = 0.01, eigen_eta = 0,
    low_rank_residual_rebuild_every = 100L,
    residual_policy = "global_projected_legacy", block_ve_mode = "fixVe",
    resam_thresh = 1.1, minimum_ve_ratio = 0.7,
    block_ve_keep_history = FALSE,
    h2 = 0.3, adjE = 0.9, nit = 1000, nburn = 100, nthin = 1,
    ncores = 1L, seed = 1L, nchains = 1L, keep_chains = FALSE,
    chain_seeds = NULL, updateB = TRUE, updateE = TRUE, updatePi = TRUE,
    nub = 4, nue = 4, pi = NULL, mixture_var = c(0, 0.01, 0.1, 1),
    alpha = NULL, comp_init = NULL, use_comp_init = FALSE, r_init = NULL,
    use_r_init = FALSE, rebuild_r_before_updateE = FALSE,
    updateE_start = 0L, updateE_every = 1L, .convergence_spec = NULL) {
  chain <- .stblr_validate_annotation_chain_args(
    nchains, keep_chains, chain_seeds)
  dims <- .stblr_get_nt_m_names(stats)
  .stblr_validate_stats(stats, dims$nt, dims$m)
  bed <- .stblr_csr_block_eigen_inputs(stats, Glist, block_start)
  if (!is.numeric(mixture_var) || length(mixture_var) < 2L ||
      mixture_var[1L] != 0 || any(!is.finite(mixture_var)) ||
      any(mixture_var[-1L] <= 0)) {
    stop("mixture_var must start with 0 and have positive non-null components.")
  }
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
  if (is.null(comp_init)) {
    comp_init <- lapply(seq_len(dims$nt), function(i) rep(0, dims$m))
  }
  if (is.null(r_init)) r_init <- stats$wy
  theta_init <- .stblr_logvar_theta_init(
    theta_init, ncol(annotation_info$X), dims$nt,
    colnames(annotation_info$X), dims$trait_names)
  raw <- stblr_cpg_omp_block_eigen_logvar_bayesr(
    stats$wy, stats$ww, stats$yy,
    lapply(seq_len(dims$nt), function(i) rep(0, dims$m)), comp_init,
    use_comp_init, r_init, use_r_init, rebuild_r_before_updateE,
    pri$B, pri$E, pri$ssb_prior_list, pri$sse_prior_list, pi, mixture_var,
    alpha, nub, nue, updateB, updateE, updatePi, adjE, as.integer(dims$n),
    as.integer(nit), as.integer(nburn), as.integer(nthin), as.integer(ncores),
    as.integer(seed), chain$nchains, chain$keep_chains, chain$chain_seeds,
    as.integer(updateE_start), as.integer(updateE_every),
    .convergence_spec$markers %||% integer(),
    isTRUE(.convergence_spec$probability), isTRUE(.convergence_spec$b),
    isTRUE(.convergence_spec$d), isTRUE(.convergence_spec$component),
    bed$bed_files, bed$n_bed, bed$cls, bed$rows, bed$af, bed$block_start,
    eigen_filter, eigen_tau, eigen_eta, representation, eigen_prop,
    as.integer(low_rank_residual_rebuild_every), list(
      phenotype_variance = as.numeric(vy), residual_policy = residual_policy,
      block_ve_mode = block_ve_mode, resam_thresh = resam_thresh,
      minimum_ve_ratio = minimum_ve_ratio,
      block_ve_keep_history = block_ve_keep_history), annotation_info$X,
    theta_init, theta_prior_sd, updateTheta)
  fit <- .as_stblr_fit(raw, dims$trait_names, dims$variable_names)
  fit$input <- list(
    n = dims$n, m = dims$m, nt = dims$nt, h2 = h2, nub = nub, nue = nue,
    nit = nit, nburn = nburn, nthin = nthin, ncores = ncores, seed = seed,
    nchains = chain$nchains, keep_chains = chain$keep_chains,
    chain_seeds = if (length(chain$chain_seeds)) chain$chain_seeds else NULL,
    updateB = updateB, updateE = updateE, updatePi = updatePi,
    mixture_var = mixture_var, alpha = alpha, updateLDswap = FALSE,
    residual_policy = residual_policy, block_ve_mode = block_ve_mode,
    resam_thresh = resam_thresh, minimum_ve_ratio = minimum_ve_ratio)
  fit <- .stblr_attach_logvar_output(
    fit, raw, annotation_info$transform, dims$variable_names,
    dims$trait_names, "sbayesr", theta_prior_sd, updateTheta)
  fit$input$backend <- "block_eigen_logvar_bayesr"
  .stblr_logvar_block_metadata(
    fit, raw, bed, representation, eigen_policy, eigen_prop, eigen_filter,
    eigen_tau, eigen_eta, low_rank_residual_rebuild_every)
}
