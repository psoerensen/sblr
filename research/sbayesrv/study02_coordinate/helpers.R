`%||%` <- function(x, y) if (is.null(x)) y else x

.study02_stop <- function(label, ...) {
  condition <- structure(
    list(message = paste0(...), call = NULL, label = label),
    class = c("study02_coordinate_error", "error", "condition"))
  stop(condition)
}

.study02_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.study02_ids_sha256 <- function(ids) {
  digest::digest(paste(ids, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

.study02_gmean <- function(x) exp(mean(log(x)))

.study02_config <- function() {
  stages <- data.frame(
    stage = c("q50_h030", "q500_h030", "q500_h050"),
    n_causal = c(50L, 500L, 500L), h2 = c(.30, .30, .50),
    architecture_seed = c(20260901L, 20260911L, 20260921L),
    training_residual_seed = c(20260902L, 20260912L, 20260922L),
    validation_residual_seed = c(20260903L, 20260913L, 20260923L),
    stringsAsFactors = FALSE)
  list(
    dataset_id = "human_independent",
    dataset_path = file.path("simulated_human_data", "human"),
    raw_n = 5000L, raw_m = 50000L, retained_m = 37991L,
    chromosome = 1L, split_seed = 3101L, train_fraction = .70,
    qc = list(excludeMAF = .05, excludeMISS = .05, excludeCGAT = TRUE,
      excludeINDEL = TRUE, excludeDUPS = TRUE, excludeHWE = 1e-12,
      excludeMHC = FALSE),
    ld = list(max_distance_bp = 0, max_distance_variants = 1000L,
      r2_threshold = .001, block_size = 1024L, nthreads = 1L),
    qgdata_head = "b840029c277703b6cd82560aac62a4377def4399",
    qgdata_checksums = c(
      bed = "831b33866e4e6382e7ce80ce99e85bbc4d3ada6e00fabd3da97c2db16a8fc74e",
      bim = "e64ba3e68e0dcb5d06a7e4af520fb9a2966696116d674735674a90eaa0a23913",
      fam = "3ec886498fc333334e5d33f259f17767a012df8c119b023b338e81b6b5fa40da"),
    ld_cache_identity =
      "8e990e54be49212386cb47b4f14023048d8ae142f211a9c3b549e031abb15565",
    annotation_seed = 20260821L,
    theta = c(binary_informative = log(4),
      continuous_informative = log(2), correlated_proxy = 0, noise = 0),
    gamma = c(null = 0, small = .01, medium = .1, large = 1),
    nit = 2000L, nburn = 250L, nthin = 1L, nchains = 4L, ncores = 4L,
    fit_seed = 30104L,
    chain_seeds = c(130104L, 230104L, 330104L, 430104L),
    theta_prior_sd = .7, stages = stages)
}

.study02_resolve_roots <- function(repo_root) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  parent <- dirname(repo_root)
  bench <- normalizePath(file.path(parent, "sblrbench"), winslash = "/",
                         mustWork = TRUE)
  qgdata <- normalizePath(file.path(parent, "qgdata"), winslash = "/",
                          mustWork = TRUE)
  list(repo = repo_root, bench = bench, qgdata = qgdata)
}

.study02_git_head <- function(path) {
  trimws(system2("git", c("-C", path, "rev-parse", "HEAD"), stdout = TRUE))
}

.study02_parse_meta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  bits <- strsplit(lines, "=", fixed = TRUE)
  keys <- vapply(bits, `[[`, character(1), 1L)
  values <- vapply(bits, function(x) paste(x[-1L], collapse = "="),
                   character(1))
  stats::setNames(values, keys)
}

.study02_validate_resource <- function(roots, config, output_dir) {
  cache_dir <- file.path(roots$repo, "research", "cache", "ld",
                         config$ld_cache_identity)
  prefix <- file.path(cache_dir, "training_ld")
  ld_cache <- file.path(cache_dir, "training_ld_glist.rds")
  manifest_path <- file.path(cache_dir, "manifest.csv")
  resources_path <- file.path(cache_dir, "resource_files.csv")
  marker_path <- file.path(cache_dir, "marker_ids.txt")
  training_path <- file.path(cache_dir, "training_ids.txt")
  validation_path <- file.path(cache_dir, "validation_ids.txt")
  af_path <- file.path(cache_dir, "training_af.csv")
  required <- c(manifest_path, resources_path, ld_cache, marker_path,
    training_path, validation_path, af_path,
    paste0(prefix, c(".meta.txt", ".row_ptr.u64.bin",
      ".col_idx.u32.0based.bin", ".values.f32.bin")))
  if (!all(file.exists(required))) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "The content-identified Study 02 LD cache is incomplete: ",
      paste(required[!file.exists(required)], collapse = ", "))
  }

  manifest_table <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  manifest <- stats::setNames(manifest_table$value, manifest_table$field)
  if (!identical(unname(manifest[["cache_identity"]]),
                 config$ld_cache_identity) ||
      !identical(unname(manifest[["qualification_status"]]),
                 "validated Study 02 coordinate resource")) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "The Study 02 LD cache identity or qualification status differs from the contract.")
  }
  resources <- utils::read.csv(resources_path, stringsAsFactors = FALSE)
  resource_paths <- file.path(cache_dir, resources$file)
  if (!all(file.exists(resource_paths)) ||
      !identical(as.numeric(file.info(resource_paths)$size),
                 as.numeric(resources$bytes)) ||
      !identical(unname(vapply(resource_paths, .study02_sha256, character(1))),
                 unname(resources$sha256))) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "A Study 02 LD cache resource is missing or fails its size/checksum contract.")
  }

  if (!identical(.study02_git_head(roots$qgdata), config$qgdata_head)) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "qgdata HEAD differs from the frozen resource contract.")
  }
  dataset_prefix <- file.path(roots$qgdata, config$dataset_path)
  genotype_files <- paste0(dataset_prefix, c(".bed", ".bim", ".fam"))
  if (!all(file.exists(genotype_files))) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "The frozen qgdata genotype files are missing.")
  }
  observed_qg_sha <- vapply(genotype_files, .study02_sha256, character(1))
  if (!identical(unname(observed_qg_sha), unname(config$qgdata_checksums))) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "The qgdata BED/BIM/FAM checksums differ from the frozen contract.")
  }

  ld_glist <- readRDS(ld_cache)
  marker_ids <- readLines(marker_path, warn = FALSE)
  training_ids <- readLines(training_path, warn = FALSE)
  validation_ids <- readLines(validation_path, warn = FALSE)
  af_table <- utils::read.csv(af_path, stringsAsFactors = FALSE)
  training_af <- stats::setNames(af_table$allele_frequency,
                                 af_table$marker_id)
  expected_split <- sblrbench::make_prediction_split(
    ld_glist$ids, config$train_fraction, config$split_seed)
  split <- list(
    train_ids = training_ids, test_ids = validation_ids,
    train_rows = match(training_ids, ld_glist$ids),
    test_rows = match(validation_ids, ld_glist$ids),
    split_seed = config$split_seed
  )
  sparse <- ld_glist$sparseLD
  if (length(ld_glist$ids) != config$raw_n ||
      length(ld_glist$rsids[[config$chromosome]]) != config$raw_m ||
      length(marker_ids) != config$retained_m || anyDuplicated(marker_ids) ||
      length(training_ids) != 3500L || length(validation_ids) != 1500L ||
      anyDuplicated(training_ids) || anyDuplicated(validation_ids) ||
      length(intersect(training_ids, validation_ids)) ||
      anyNA(split$train_rows) || anyNA(split$test_rows) ||
      !identical(training_ids, expected_split$train_ids) ||
      !identical(validation_ids, expected_split$test_ids) ||
      !identical(as.character(sparse$marker_names), marker_ids) ||
      !identical(as.integer(sparse$rows), as.integer(split$train_rows)) ||
      !identical(as.integer(sparse$reference_n), 3500L) ||
      !identical(names(training_af), marker_ids) ||
      length(training_af) != config$retained_m ||
      any(!is.finite(training_af) | training_af <= 0 | training_af >= 1)) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "Cached markers, samples, or training allele frequencies differ from Study 02.")
  }
  meta <- .study02_parse_meta(paste0(prefix, ".meta.txt"))
  expected_meta <- c(n_used = "3500", n_variants = "37991",
    r2_threshold = "0.001", max_distance_bp = "0",
    max_distance_variants = "1000", block_size = "1024", nthreads = "1")
  if (!all(meta[names(expected_meta)] == expected_meta)) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "Cached CSR metadata differs from the exact Study 02 LD policy.")
  }
  af_sha256 <- digest::digest(
    paste(sprintf("%.17g", as.numeric(sparse$af[[1L]])), collapse = "\n"),
    algo = "sha256", serialize = FALSE)
  if (!isTRUE(all.equal(unname(as.numeric(sparse$af[[1L]])),
                        unname(training_af), tolerance = 1e-15)) ||
      !identical(unname(manifest[["training_af_sha256"]]), af_sha256)) {
    .study02_stop("BLOCKED_MISSING_STUDY02_LD_RESOURCE",
      "The RDS and manifest training allele frequencies differ.")
  }
  ld_glist$bedfiles <- genotype_files[[1L]]
  ld_glist$bimfiles <- genotype_files[[2L]]
  ld_glist$famfiles <- genotype_files[[3L]]
  ld_glist$sparseLD$bed_files <- genotype_files[[1L]]
  ld_glist$sparseLD$prefix <- prefix
  resource_manifest <- data.frame(
    dataset_id = config$dataset_id, dataset_path = config$dataset_path,
    qgdata_head = config$qgdata_head, raw_samples = config$raw_n,
    raw_markers = config$raw_m, retained_markers = length(marker_ids),
    marker_ids_sha256 = .study02_ids_sha256(marker_ids),
    training = length(split$train_ids), validation = length(split$test_ids),
    training_ids_sha256 = .study02_ids_sha256(split$train_ids),
    validation_ids_sha256 = .study02_ids_sha256(split$test_ids),
    split_seed = split$split_seed, max_distance_bp = 0,
    max_distance_variants = 1000L, r2_threshold = .001,
    block_size = 1024L, nthreads = 1L,
    allele_frequency_policy = "training_subset",
    cache_identity = config$ld_cache_identity,
    ld_prefix = file.path("research", "cache", "ld",
      config$ld_cache_identity, "training_ld"),
    ld_rebuilt = FALSE, stringsAsFactors = FALSE)
  utils::write.csv(resources, file.path(output_dir,
    "reused_resource_files.csv"), row.names = FALSE)
  utils::write.csv(resource_manifest,
                   file.path(output_dir, "resource_manifest.csv"),
                   row.names = FALSE)
  list(glist = ld_glist, ld_glist = ld_glist, prefix = prefix,
    marker_ids = marker_ids, sample_ids = ld_glist$ids, split = split,
    training_af = training_af, manifest = resource_manifest)
}

.study02_load_scaled_genotypes <- function(resource, config) {
  get <- function(ids) qgg::getG(
    Glist = resource$ld_glist, chr = config$chromosome, ids = ids,
    rsids = resource$marker_ids, impute = TRUE, scale = TRUE)
  train <- get(resource$split$train_ids)
  validation <- get(resource$split$test_ids)
  if (!identical(rownames(train), resource$split$train_ids) ||
      !identical(rownames(validation), resource$split$test_ids) ||
      !identical(colnames(train), resource$marker_ids) ||
      !identical(colnames(validation), resource$marker_ids) ||
      any(!is.finite(train)) || any(!is.finite(validation))) {
    .study02_stop("BLOCKED_Q50_H030",
      "Training-scale genotype extraction lost canonical alignment.")
  }
  training_mean_error <- max(abs(colMeans(train)))
  if (!is.finite(training_mean_error) || training_mean_error > 1e-10) {
    .study02_stop("BLOCKED_Q50_H030",
      "Training-subset allele-frequency centring was not reproduced.")
  }
  list(train = train, validation = validation,
       training_mean_error = training_mean_error)
}

.study02_annotations <- function(marker_ids, config) {
  m <- length(marker_ids)
  set.seed(config$annotation_seed)
  binary <- numeric(m)
  binary[sample.int(m, as.integer(round(.30 * m)))] <- 1
  bounded <- seq(-1, 1, length.out = m)
  continuous <- sample(bounded, m, replace = FALSE)
  auxiliary <- sample(bounded, m, replace = FALSE)
  noise <- sample(bounded, m, replace = FALSE)
  proxy <- .70 * as.numeric(scale(continuous)) +
    sqrt(1 - .70^2) * as.numeric(scale(auxiliary))
  raw <- cbind(binary_informative = binary,
    continuous_informative = continuous, correlated_proxy = proxy,
    noise = noise)
  rownames(raw) <- marker_ids
  processed <- sbayesrv_preprocess_annotations(raw, marker_ids)
  eta_q <- sbayesrv_eta_q(config$theta, processed$X)
  q <- stats::setNames(as.numeric(eta_q$q), marker_ids)
  if (!identical(rownames(processed$X), marker_ids) ||
      any(!is.finite(q) | q <= 0) ||
      abs(log(.study02_gmean(q))) > 1e-12) {
    .study02_stop("BLOCKED_Q50_H030",
      "The frozen annotation/q surface failed its contract.")
  }
  list(raw = raw, X = processed$X, q = q,
    proxy_correlation = stats::cor(processed$X[, "continuous_informative"],
      processed$X[, "correlated_proxy"]))
}

.study02_stage_pi <- function(n_causal, m) {
  if (n_causal == 50L) {
    c(null = .99, small = .01 / 3, medium = .01 / 3, large = .01 / 3)
  } else {
    active <- n_causal / m
    c(null = 1 - active, small = .60 * active,
      medium = .30 * active, large = .10 * active)
  }
}

.study02_standard_residual <- function(n, seed, variance, ids) {
  set.seed(seed)
  e <- stats::rnorm(n)
  e <- (e - mean(e)) / stats::sd(e) * sqrt(variance)
  stats::setNames(e, ids)
}

.study02_simulate_stage <- function(stage, W, annotations, config) {
  pi <- .study02_stage_pi(stage$n_causal, ncol(W$train))
  vg <- stage$h2 / (1 - stage$h2)
  common <- list(W = W$train, A = annotations$raw,
    architecture = "bayesr", nt = 1L, h2 = stage$h2, vg = vg,
    pi = pi, mixture_variances = config$gamma,
    n_causal = stage$n_causal, annotation_model = "none",
    seed = stage$architecture_seed, standardize_W = FALSE,
    scale_effects = TRUE, return_genotypes = FALSE,
    return_marker_probabilities = TRUE, compute_sumstats = FALSE)
  unit <- tryCatch(
    do.call(gsim, c(common, list(marker_multipliers = NULL))),
    error = function(e) .study02_stop(
      paste0("BLOCKED_", toupper(stage$stage)),
      "Unit-q companion simulation failed: ", conditionMessage(e)))
  unit_rng <- .Random.seed
  truth <- tryCatch(
    do.call(gsim, c(common, list(marker_multipliers = annotations$q))),
    error = function(e) .study02_stop(
      paste0("BLOCKED_", toupper(stage$stage)),
      "Variance-modulated simulation failed: ", conditionMessage(e)))
  q_rng <- .Random.seed

  active <- which(truth$component != 1L)
  null <- which(truth$component == 1L)
  unit_pre <- unit$B[, 1L] / unit$settings$effect_scale[[1L]]
  q_pre <- truth$B[, 1L] / truth$settings$effect_scale[[1L]]
  expected_q_pre <- unit_pre * sqrt(annotations$q)
  factor_error <- max(abs(q_pre - expected_q_pre))
  rng_equal <- identical(unit_rng, q_rng)
  state_equal <- identical(unit$component, truth$component) &&
    identical(unit$causal_rsids, truth$causal_rsids)
  z_unit <- unit_pre[active] /
    sqrt(config$gamma[unit$component[active]])
  z_q <- q_pre[active] /
    sqrt(config$gamma[truth$component[active]] * annotations$q[active])
  z_error <- max(abs(z_unit - z_q))
  probability_constant <- max(abs(truth$marker_probabilities -
    matrix(pi, nrow(truth$marker_probabilities), length(pi), byrow = TRUE)))

  effects <- stats::setNames(truth$B[, 1L], colnames(W$train))
  G_train <- as.numeric(W$train %*% effects)
  names(G_train) <- rownames(W$train)
  g_error <- max(abs(G_train - as.numeric(truth$G[, 1L])))
  residual_variance <- vg * (1 - stage$h2) / stage$h2
  E_train <- .study02_standard_residual(nrow(W$train),
    stage$training_residual_seed, residual_variance, rownames(W$train))
  Y_train <- G_train + E_train
  G_validation <- as.numeric(W$validation %*% effects)
  names(G_validation) <- rownames(W$validation)
  E_validation <- .study02_standard_residual(nrow(W$validation),
    stage$validation_residual_seed, residual_variance,
    rownames(W$validation))
  Y_validation <- G_validation + E_validation
  observed_vg <- stats::var(G_train)
  observed_h2 <- observed_vg / (observed_vg + stats::var(E_train))
  y_error <- max(abs(Y_train - G_train - E_train))

  checks <- data.frame(
    stage = stage$stage,
    requested_causal = stage$n_causal,
    observed_causal = length(active),
    null_effect_max_abs = if (length(null)) max(abs(effects[null])) else 0,
    component_probability_sum_error = abs(sum(pi) - 1),
    probability_surface_max_error = probability_constant,
    marker_multiplier_geometric_mean = .study02_gmean(annotations$q),
    marker_multiplier_factor_max_error = factor_error,
    underlying_z_max_error = z_error,
    unit_q_rng_identical = rng_equal,
    unit_q_states_identical = state_equal,
    genetic_value_max_error = g_error,
    phenotype_identity_max_error = y_error,
    target_vg = vg, realized_training_vg = observed_vg,
    target_h2 = stage$h2, realized_training_h2 = observed_h2,
    train_validation_overlap = length(intersect(rownames(W$train),
                                                 rownames(W$validation))),
    stringsAsFactors = FALSE)
  passed <- length(active) == stage$n_causal &&
    identical(names(truth$component), colnames(W$train)) &&
    identical(names(truth$marker_multipliers), colnames(W$train)) &&
    max(abs(effects[null])) == 0 && abs(sum(pi) - 1) < 1e-14 &&
    probability_constant < 1e-14 && abs(log(.study02_gmean(annotations$q))) < 1e-12 &&
    factor_error < 1e-12 && z_error < 1e-12 && rng_equal && state_equal &&
    g_error < 1e-9 && y_error < 1e-12 &&
    abs(observed_vg - vg) < 1e-10 &&
    abs(observed_h2 - stage$h2) < 1e-12 &&
    checks$train_validation_overlap == 0
  if (!passed) {
    .study02_stop(paste0("BLOCKED_", toupper(stage$stage)),
      "A mandatory simulation identity failed before fitting.")
  }
  component <- stats::setNames(as.integer(truth$component), colnames(W$train))
  list(stage = stage, pi = pi, effects = effects, component = component,
    causal = effects != 0, G_train = G_train, E_train = E_train,
    Y_train = Y_train, G_validation = G_validation,
    E_validation = E_validation, Y_validation = Y_validation,
    checks = checks, effect_scale = truth$settings$effect_scale[[1L]])
}

.study02_make_stats <- function(simulation, resource, config) {
  y <- matrix(simulation$Y_train, ncol = 1L,
    dimnames = list(resource$split$train_ids, "trait1"))
  stats <- make_summary_stats(Glist = resource$ld_glist, y = y,
    chr = config$chromosome, rows = resource$split$train_rows,
    scale = TRUE, nthreads = 1L)
  reference_af <- resource$training_af
  if (!identical(stats$marker_names, resource$marker_ids) ||
      !identical(stats$trait_names, "trait1") ||
      !identical(as.integer(stats$n), 3500L) ||
      !isTRUE(all.equal(unname(stats$af[[1L]]), unname(reference_af),
                        tolerance = 0))) {
    .study02_stop(paste0("BLOCKED_", toupper(simulation$stage$stage)),
      "Phenotype-specific summary statistics differ from the Study 02 scale.")
  }
  stats
}

.study02_extract_fit <- function(fit, marker_ids, method, stage) {
  effect <- as.matrix(extract_posterior(fit, "realised_effects"))
  pip <- as.matrix(extract_posterior(fit, "pips"))
  component <- extract_posterior(fit, "component_probabilities")
  if (is.list(component) && !is.array(component)) component <- component[[1L]]
  component <- as.matrix(component)
  posterior_summary <- summarise_posterior(fit)
  diagnostics <- extract_diagnostics(fit)
  consistency <- check_stblr_consistency(fit, require_chains = TRUE,
                                         verbose = FALSE)
  residual <- posterior_summary[posterior_summary$parameter == "ve", ,
                                drop = FALSE]
  aligned <- identical(rownames(effect), marker_ids) &&
    identical(rownames(pip), marker_ids) &&
    identical(rownames(component), marker_ids)
  valid <- isTRUE(consistency$ok) && aligned && ncol(effect) == 1L &&
    ncol(pip) == 1L && all(is.finite(effect)) && all(is.finite(pip)) &&
    all(is.finite(component)) && max(abs(rowSums(component) - 1)) < 1e-8 &&
    nrow(residual) == 1L && is.finite(residual$min) && residual$min > 0 &&
    is.list(diagnostics)
  if (!valid) {
    .study02_stop(paste0("BLOCKED_", toupper(stage)),
      method, " failed canonical posterior or diagnostic validation.")
  }
  list(effect = stats::setNames(effect[, 1L], marker_ids),
    pip = stats::setNames(pip[, 1L], marker_ids), component = component,
    summary = posterior_summary, diagnostics = diagnostics,
    consistency = consistency,
    minimum_residual_variance = residual$min)
}

.study02_fit_stage <- function(simulation, stats, resource, annotations,
                               config) {
  common <- list(stats = stats, Glist = resource$ld_glist,
    ld_prefix = resource$prefix, nit = config$nit, nburn = config$nburn,
    nthin = config$nthin, nchains = config$nchains, ncores = config$ncores,
    seed = config$fit_seed, chain_seeds = config$chain_seeds,
    keep_chains = TRUE, convergence = "core",
    convergence_control = list(warn = FALSE, keep_traces = TRUE),
    h2 = simulation$stage$h2, pi = simulation$pi,
    mixture_var = config$gamma, updateLDswap = FALSE, verbose = FALSE)
  start <- proc.time()[["elapsed"]]
  ordinary_fit <- tryCatch(
    do.call(stblr_csr, c(common, list(method = "sbayesr"))),
    error = function(e) .study02_stop(
      paste0("BLOCKED_", toupper(simulation$stage$stage)),
      "Ordinary SBayesR failed: ", conditionMessage(e)))
  ordinary_seconds <- proc.time()[["elapsed"]] - start
  ordinary <- .study02_extract_fit(ordinary_fit, resource$marker_ids,
    "ordinary SBayesR", simulation$stage$stage)

  start <- proc.time()[["elapsed"]]
  joint_fit <- tryCatch(
    do.call(stblr_csr_annot, c(common, list(
      annotations = annotations$raw, annotation_model = "log_variance",
      method = "sbayesr", theta_prior_sd = config$theta_prior_sd,
      theta_init = NULL, updateTheta = TRUE))),
    error = function(e) .study02_stop(
      paste0("BLOCKED_", toupper(simulation$stage$stage)),
      "SBayesRV failed: ", conditionMessage(e)))
  joint_seconds <- proc.time()[["elapsed"]] - start
  joint <- .study02_extract_fit(joint_fit, resource$marker_ids,
    "SBayesRV", simulation$stage$stage)
  required <- c("theta_summary", "theta_trace", "marker_prior_scale")
  if (!all(required %in% names(joint_fit)) ||
      !identical(rownames(joint_fit$marker_prior_scale), resource$marker_ids) ||
      !all(config$theta |> names() %in% joint_fit$theta_summary$annotation)) {
    .study02_stop(paste0("BLOCKED_", toupper(simulation$stage$stage)),
      "SBayesRV documented theta/q output is incomplete or misaligned.")
  }
  joint$theta_summary <- joint_fit$theta_summary
  joint$theta_trace <- joint_fit$theta_trace
  joint$marker_prior_scale <- stats::setNames(
    as.numeric(joint_fit$marker_prior_scale[, 1L]), resource$marker_ids)
  list(ordinary = ordinary, joint = joint,
    runtime = data.frame(stage = simulation$stage$stage,
      method = c("ordinary_sbayesr", "sbayesrv"),
      seconds = c(ordinary_seconds, joint_seconds), stringsAsFactors = FALSE))
}

.study02_average_precision <- function(score, truth) {
  ordering <- order(score, decreasing = TRUE)
  hit <- as.logical(truth[ordering])
  mean(cumsum(hit)[hit] / which(hit))
}

.study02_method_metrics <- function(stage, method, extracted, simulation, W) {
  prediction <- as.numeric(W$validation %*% extracted$effect)
  truth_g <- simulation$G_validation
  causal <- simulation$causal
  fit_calibration <- stats::lm(truth_g ~ prediction)
  ranks <- rank(-extracted$pip, ties.method = "min")
  metrics <- c(
    genetic_value_prediction_correlation = stats::cor(prediction, truth_g),
    genetic_value_prediction_mse = mean((prediction - truth_g)^2),
    genetic_value_prediction_nmse = mean((prediction - truth_g)^2) /
      stats::var(truth_g),
    phenotype_prediction_correlation = stats::cor(prediction,
      simulation$Y_validation),
    calibration_intercept = stats::coef(fit_calibration)[[1L]],
    calibration_slope = stats::coef(fit_calibration)[[2L]],
    effect_rmse = sqrt(mean((extracted$effect - simulation$effects)^2)),
    pip_brier_score = mean((extracted$pip - causal)^2),
    average_precision = .study02_average_precision(extracted$pip, causal),
    top_50_causal_recall = mean(ranks[causal] <= 50),
    top_100_causal_recall = mean(ranks[causal] <= 100),
    top_500_causal_recall = mean(ranks[causal] <= 500),
    top_1000_causal_recall = mean(ranks[causal] <= 1000))
  list(metrics = data.frame(stage = stage, method = method,
      metric = names(metrics), value = as.numeric(metrics),
      stringsAsFactors = FALSE),
    causal_ranks = data.frame(stage = stage, method = method,
      marker_id = names(ranks)[causal], pip = extracted$pip[causal],
      rank = as.numeric(ranks[causal]),
      true_q = NA_real_, component = simulation$component[causal],
      stringsAsFactors = FALSE))
}

.study02_theta_q <- function(joint, annotations, simulation, config) {
  theta_summary <- joint$theta_summary
  theta_summary$truth <- config$theta[match(theta_summary$annotation,
                                             names(config$theta))]
  theta_mean <- stats::setNames(theta_summary$mean, theta_summary$annotation)
  log_q_mean <- as.numeric(annotations$X %*% theta_mean[colnames(annotations$X)])
  q_from_mean_theta <- exp(log_q_mean)
  trace <- joint$theta_trace
  draws <- seq.int(config$nburn + 1L, dim(trace)[1L])
  theta_draws <- do.call(rbind, lapply(seq_len(dim(trace)[3L]), function(chain)
    trace[draws, , chain, drop = FALSE][, , 1L]))
  colnames(theta_draws) <- colnames(annotations$X)
  q_sum <- numeric(nrow(annotations$X))
  chunk <- 50L
  for (start in seq.int(1L, nrow(theta_draws), by = chunk)) {
    take <- start:min(start + chunk - 1L, nrow(theta_draws))
    q_sum <- q_sum + rowSums(exp(annotations$X %*%
      t(theta_draws[take, , drop = FALSE])))
  }
  q_posterior_mean <- q_sum / nrow(theta_draws)
  true_log_q <- log(annotations$q)
  causal <- simulation$causal
  q_summary <- data.frame(
    stage = simulation$stage$stage,
    log_q_correlation = stats::cor(true_log_q, log_q_mean),
    log_q_rmse = sqrt(mean((true_log_q - log_q_mean)^2)),
    causal_log_q_correlation = stats::cor(true_log_q[causal],
      log_q_mean[causal]),
    causal_log_q_rmse = sqrt(mean((true_log_q[causal] -
      log_q_mean[causal])^2)),
    posterior_mean_q_min = min(q_posterior_mean),
    posterior_mean_q_max = max(q_posterior_mean),
    posterior_mean_q_geometric_mean = .study02_gmean(q_posterior_mean),
    q_at_mean_theta_geometric_mean = .study02_gmean(q_from_mean_theta),
    draw_q_geometric_mean_min = exp(min(rowSums(sweep(theta_draws, 2L,
      colMeans(annotations$X), "*")))),
    draw_q_geometric_mean_max = exp(max(rowSums(sweep(theta_draws, 2L,
      colMeans(annotations$X), "*")))),
    stringsAsFactors = FALSE)
  list(theta = theta_summary, q_summary = q_summary,
       q_posterior_mean = stats::setNames(q_posterior_mean,
                                           rownames(annotations$X)))
}

.study02_write_summary <- function(output_dir, resource, annotations,
                                    config, status, attempted) {
  lines <- c(
    "# Study 02-coordinate SBayesRV research run",
    "",
    paste0("Status: `", status, "`"),
    "",
    paste0("Reused CSR prefix: `", resource$manifest$ld_prefix, "`."),
    "LD was validated and reused from the content-identified research cache;",
    "it was not rebuilt or copied.",
    paste0("Annotation proxy correlation: ",
      format(annotations$proxy_correlation, digits = 6), "."),
    paste0("True q range: ", format(min(annotations$q), digits = 6),
      " to ", format(max(annotations$q), digits = 6),
      "; geometric mean ", format(.study02_gmean(annotations$q), digits = 8),
      "."),
    paste0("Attempted stages: ", paste(attempted, collapse = ", "), "."),
    "",
    "This is one research replicate. Mechanical success and descriptive paired",
    "differences are not benchmark evidence or a general performance claim.")
  writeLines(lines, file.path(output_dir, "summary.md"))
}

run_study02_coordinate <- function(repo_root = ".") {
  config <- .study02_config()
  roots <- .study02_resolve_roots(repo_root)
  output_dir <- file.path(roots$repo, "research", "sbayesrv",
    "study02_coordinate", "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  resource <- .study02_validate_resource(roots, config, output_dir)
  W <- .study02_load_scaled_genotypes(resource, config)
  annotations <- .study02_annotations(resource$marker_ids, config)
  utils::write.csv(data.frame(
    annotation = colnames(annotations$X), theta_true = config$theta,
    processed_mean = colMeans(annotations$X),
    stringsAsFactors = FALSE), file.path(output_dir, "annotation_truth.csv"),
    row.names = FALSE)
  utils::write.csv(data.frame(training_column_mean_max_abs =
    W$training_mean_error, proxy_correlation = annotations$proxy_correlation,
    q_min = min(annotations$q), q_median = stats::median(annotations$q),
    q_max = max(annotations$q), q_geometric_mean = .study02_gmean(annotations$q),
    stringsAsFactors = FALSE), file.path(output_dir,
      "annotation_surface_summary.csv"), row.names = FALSE)

  truth_checks <- metrics <- ranks <- fit_status <- runtimes <-
    theta <- q_recovery <- component_q <- list()
  attempted <- character()
  status <- "SBayesRV_PREDICTION_SIGNAL_OBSERVED_AT_STUDY02_COORDINATE"
  utils::write.csv(config$stages, file.path(output_dir, "stage_design.csv"),
                   row.names = FALSE)
  for (i in seq_len(nrow(config$stages))) {
    stage <- config$stages[i, , drop = FALSE]
    attempted <- c(attempted, stage$stage)
    message("[", stage$stage, "] validating simulation identities")
    simulation <- .study02_simulate_stage(stage, W, annotations, config)
    truth_checks[[i]] <- simulation$checks
    utils::write.csv(do.call(rbind, truth_checks), file.path(output_dir,
      "truth_checks.csv"), row.names = FALSE)
    stats <- .study02_make_stats(simulation, resource, config)
    message("[", stage$stage, "] fitting ordinary SBayesR")
    fits <- .study02_fit_stage(simulation, stats, resource, annotations, config)
    message("[", stage$stage, "] both fits validated; evaluating validation set")
    ordinary <- .study02_method_metrics(stage$stage, "ordinary_sbayesr",
      fits$ordinary, simulation, W)
    joint <- .study02_method_metrics(stage$stage, "sbayesrv", fits$joint,
      simulation, W)
    ordinary$causal_ranks$true_q <- annotations$q[ordinary$causal_ranks$marker_id]
    joint$causal_ranks$true_q <- annotations$q[joint$causal_ranks$marker_id]
    metrics[[i]] <- rbind(ordinary$metrics, joint$metrics)
    ranks[[i]] <- rbind(ordinary$causal_ranks, joint$causal_ranks)
    runtimes[[i]] <- fits$runtime
    tq <- .study02_theta_q(fits$joint, annotations, simulation, config)
    tq$theta$stage <- stage$stage
    theta[[i]] <- tq$theta
    q_recovery[[i]] <- tq$q_summary
    causal_component <- split(which(simulation$causal),
      simulation$component[simulation$causal])
    component_q[[i]] <- do.call(rbind, lapply(names(causal_component), function(k) {
      idx <- causal_component[[k]]
      data.frame(stage = stage$stage, component = as.integer(k), n = length(idx),
        true_log_q_mean = mean(log(annotations$q[idx])),
        estimated_log_q_mean = mean(log(tq$q_posterior_mean[idx])),
        stringsAsFactors = FALSE)
    }))
    fit_status[[i]] <- data.frame(stage = stage$stage,
      method = c("ordinary_sbayesr", "sbayesrv"), status = "ok",
      minimum_residual_variance = c(fits$ordinary$minimum_residual_variance,
        fits$joint$minimum_residual_variance),
      canonical_consistency = c(fits$ordinary$consistency$ok,
        fits$joint$consistency$ok), stringsAsFactors = FALSE)
    utils::write.csv(do.call(rbind, metrics), file.path(output_dir,
      "metrics.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, ranks), file.path(output_dir,
      "causal_ranks.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, runtimes), file.path(output_dir,
      "runtime.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, theta), file.path(output_dir,
      "theta_summary.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, q_recovery), file.path(output_dir,
      "q_recovery.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, component_q), file.path(output_dir,
      "q_recovery_by_component.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, fit_status), file.path(output_dir,
      "fit_status.csv"), row.names = FALSE)
    differences <- do.call(rbind, lapply(split(do.call(rbind, metrics),
      do.call(rbind, metrics)$stage), function(x) {
        wide <- reshape(x, idvar = c("stage", "metric"), timevar = "method",
                        direction = "wide")
        data.frame(stage = wide$stage, metric = wide$metric,
          sbayesrv_minus_ordinary = wide$value.sbayesrv -
            wide$value.ordinary_sbayesr,
          orientation = ifelse(grepl("mse|rmse|brier", wide$metric),
            "negative_favours_sbayesrv", "positive_favours_sbayesrv"),
          stringsAsFactors = FALSE)
      }))
    utils::write.csv(differences, file.path(output_dir,
      "paired_differences.csv"), row.names = FALSE)
    rm(simulation, stats, fits); gc(verbose = FALSE)
  }
  utils::write.csv(data.frame(status = status, stages_attempted =
    paste(attempted, collapse = ";"), ld_rebuilt = FALSE,
    baseline_refitted = FALSE, stringsAsFactors = FALSE),
    file.path(output_dir, "status.csv"), row.names = FALSE)
  .study02_write_summary(output_dir, resource, annotations, config, status,
                         attempted)
  status
}
