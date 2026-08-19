# Study 06-specific annotation construction and annotation-driven simulation.

annotation_zscore <- function(x) {
  z <- (x - mean(x)) / stats::sd(x)
  if (any(!is.finite(z)))
    stop("Annotation standardization failed.", call. = FALSE)
  as.numeric(z)
}

select_study06_v2_blocks <- function(marker_ids, positions, spec) {
  marker_ids <- as.character(marker_ids)
  positions <- as.numeric(positions)
  design <- spec$data$block_design
  if (length(marker_ids) != length(positions) || anyNA(marker_ids) ||
      anyDuplicated(marker_ids) || any(!is.finite(positions)))
    stop("QC marker IDs and positions must be unique, finite, and aligned.",
      call. = FALSE)
  genomic_order <- order(positions, marker_ids)
  marker_ids <- marker_ids[genomic_order]
  positions <- positions[genomic_order]
  block_count <- as.integer(design$block_count)
  block_size <- as.integer(design$block_size)
  if (length(marker_ids) < block_count * block_size)
    stop("Too few QC-qualified markers for the prespecified v2 blocks.",
      call. = FALSE)
  boundaries <- floor(seq(0, length(marker_ids), length.out = block_count + 1L))
  selected <- lapply(seq_len(block_count), function(block) {
    lo <- boundaries[block] + 1L
    hi <- boundaries[block + 1L]
    available <- hi - lo + 1L
    if (available < block_size)
      stop("A genomic stratum cannot supply the prespecified block size.",
        call. = FALSE)
    start <- lo + floor((available - block_size) / 2)
    index <- seq.int(start, length.out = block_size)
    data.frame(marker_id = marker_ids[index], position_bp = positions[index],
      qc_index = genomic_order[index], block_id = block,
      within_block_index = seq_len(block_size),
      stringsAsFactors = FALSE)
  })
  panel <- do.call(rbind, selected)
  rownames(panel) <- NULL
  if (nrow(panel) != design$marker_target || anyDuplicated(panel$marker_id) ||
      !identical(as.integer(table(panel$block_id)),
        rep(block_size, block_count)))
    stop("Deterministic v2 block selection violated its marker contract.",
      call. = FALSE)
  panel
}

audit_study06_v2_blocks <- function(panel, scaled_train,
                                     method_marker_ids = list()) {
  if (!is.data.frame(panel) || !all(c("marker_id", "position_bp", "block_id",
      "within_block_index") %in% names(panel)) ||
      !is.matrix(scaled_train) ||
      !identical(colnames(scaled_train), panel$marker_id))
    stop("Block audit inputs are incomplete or marker-misaligned.", call. = FALSE)
  block_rows <- lapply(split(seq_len(nrow(panel)), panel$block_id), function(idx) {
    correlation <- stats::cor(scaled_train[, idx, drop = FALSE])
    eigenvalue <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
    positive <- eigenvalue[eigenvalue > max(eigenvalue) * .Machine$double.eps *
      length(eigenvalue)]
    retained <- if (length(positive))
      which(cumsum(positive) / sum(positive) >= .995)[1L] else 0L
    data.frame(block_id = panel$block_id[idx[1L]], marker_count = length(idx),
      first_marker = panel$marker_id[idx[1L]], last_marker = panel$marker_id[idx[length(idx)]],
      start_bp = min(panel$position_bp[idx]), end_bp = max(panel$position_bp[idx]),
      physical_span_bp = diff(range(panel$position_bp[idx])),
      within_block_rank = qr(scaled_train[, idx, drop = FALSE])$rank,
      minimum_positive_eigenvalue = if (length(positive)) min(positive) else NA_real_,
      maximum_positive_eigenvalue = if (length(positive)) max(positive) else NA_real_,
      retained_eigen_count_0_995 = retained,
      retained_positive_mass = if (retained) sum(positive[seq_len(retained)]) /
        sum(positive) else NA_real_, stringsAsFactors = FALSE)
  })
  block_table <- do.call(rbind, block_rows)
  all_correlation <- abs(stats::cor(scaled_train))
  block <- panel$block_id
  cross <- outer(block, block, `!=`) & upper.tri(all_correlation)
  equality <- if (!length(method_marker_ids)) TRUE else
    all(vapply(method_marker_ids, identical, logical(1), panel$marker_id))
  list(blocks = block_table,
    summary = data.frame(marker_count = nrow(panel),
      block_count = length(unique(panel$block_id)),
      maximum_cross_block_absolute_correlation = max(all_correlation[cross]),
      method_marker_identity_and_order_equal = equality,
      eigen_prop = .995, stringsAsFactors = FALSE))
}

validate_annotation_design <- function(A, marker_ids, spec) {
  required <- spec$annotation_design$columns
  if (!is.matrix(A) || !is.numeric(A))
    stop("Annotation design must be a numeric matrix.", call. = FALSE)
  if (!identical(rownames(A), marker_ids))
    stop("Annotation marker IDs or order differ from the benchmark markers.",
      call. = FALSE)
  if (!identical(colnames(A), required))
    stop("Annotation columns must preserve this exact order: ",
      paste(required, collapse = ", "), ".", call. = FALSE)
  if (any(!is.finite(A)))
    stop("Annotation design contains missing or non-finite values.",
      call. = FALSE)
  if (any(A[, "Intercept"] != 1))
    stop("The annotation intercept must equal one for every marker.",
      call. = FALSE)
  if (!all(A[, "enriched_binary"] %in% c(0, 1)))
    stop("enriched_binary must contain only zero and one.", call. = FALSE)
  expected <- max(1L, round(length(marker_ids) *
    spec$annotation_design$enriched_fraction))
  if (sum(A[, "enriched_binary"]) != expected)
    stop("The enriched binary group size differs from the audited design.",
      call. = FALSE)
  continuous <- A[, c("continuous_signal", "null_annotation"), drop = FALSE]
  if (any(abs(colMeans(continuous)) > 1e-10) ||
      any(abs(apply(continuous, 2L, stats::sd) - 1) > 1e-10))
    stop("Continuous annotations must use the audited sample z-score rule.",
      call. = FALSE)
  if (qr(A)$rank != ncol(A))
    stop("Annotation design is not full column rank.", call. = FALSE)
  invisible(TRUE)
}

construct_annotation_design <- function(marker_ids, spec) {
  if (!is.character(marker_ids) || anyNA(marker_ids) ||
      anyDuplicated(marker_ids))
    stop("marker_ids must be unique non-missing strings.", call. = FALSE)
  set.seed(as.integer(spec$annotation_design$seed))
  m <- length(marker_ids)
  selected <- sort(sample.int(m, max(1L, round(m *
    spec$annotation_design$enriched_fraction))))
  enriched <- integer(m)
  enriched[selected] <- 1L
  A <- cbind(Intercept = 1,
    enriched_binary = enriched,
    continuous_signal = annotation_zscore(stats::rnorm(m)),
    null_annotation = annotation_zscore(stats::rnorm(m)))
  rownames(A) <- marker_ids
  storage.mode(A) <- "double"
  validate_annotation_design(A, marker_ids, spec)
  A
}

annotation_marker_probabilities <- function(A, alpha, mixture_var) {
  out <- sblr::sbayesrc_marker_pi(A, alpha, gamma = mixture_var)
  if (!identical(dim(out), c(nrow(A), length(mixture_var))) ||
      any(!is.finite(out)) || any(out < 0 | out > 1) ||
      any(abs(rowSums(out) - 1) > 1e-12))
    stop("Invalid annotation-implied component probabilities.",
      call. = FALSE)
  rownames(out) <- rownames(A)
  colnames(out) <- paste0("component_", seq_len(ncol(out)) - 1L)
  out
}

annotation_reverse_sticks <- function(component_probability) {
  p <- as.numeric(component_probability)
  if (length(p) < 2L || any(!is.finite(p)) || any(p <= 0) ||
      abs(sum(p) - 1) > 1e-10)
    stop("Invalid marginal component probabilities.", call. = FALSE)
  remaining <- rev(cumsum(rev(p)))
  stats::qnorm(vapply(seq_len(length(p) - 1L),
    function(j) remaining[j + 1L] / remaining[j], numeric(1)))
}

construct_annotation_truth <- function(A, spec) {
  simulation <- spec$controls$simulation
  target_active <- simulation$target_expected_nonnull / nrow(A)
  target_marginal <- c(1 - target_active,
    target_active * simulation$active_component_weights)
  nonintercept <- spec$annotation_design$informative_nonintercept_alpha
  objective <- function(intercept) {
    marginal <- colMeans(annotation_marker_probabilities(A,
      rbind(Intercept = intercept, nonintercept), simulation$mixture_var))
    sum((marginal - target_marginal)^2)
  }
  initial <- annotation_reverse_sticks(target_marginal)
  calibrated <- stats::optim(initial, objective, method = "BFGS",
    control = list(reltol = 1e-14, maxit = 2000L))
  if (calibrated$convergence != 0L || calibrated$value > 1e-12)
    stop("Informative marginal component calibration failed.", call. = FALSE)
  informative <- rbind(Intercept = calibrated$par, nonintercept)
  informative_pi <- annotation_marker_probabilities(A, informative,
    simulation$mixture_var)
  marginal <- colMeans(informative_pi)
  uninformative <- matrix(0, nrow = ncol(A),
    ncol = length(simulation$mixture_var) - 1L,
    dimnames = list(colnames(A), paste0("step_", 1:3)))
  uninformative["Intercept", ] <- annotation_reverse_sticks(marginal)
  uninformative_pi <- annotation_marker_probabilities(A, uninformative,
    simulation$mixture_var)
  if (max(abs(colMeans(uninformative_pi) - marginal)) > 1e-10)
    stop("Uninformative marginal calibration failed.", call. = FALSE)
  enriched <- A[, "enriched_binary"] == 1
  expected_nonnull <- 1 - informative_pi[, 1L]
  share <- sum(expected_nonnull[enriched]) / sum(expected_nonnull)
  share_range <- spec$annotation_design$enriched_expected_nonnull_share_range
  if (abs(sum(expected_nonnull) - simulation$target_expected_nonnull) > .01 ||
      share < share_range[1L] || share > share_range[2L] ||
      any(uninformative[-1L, , drop = FALSE] != 0))
    stop("Annotation truth does not satisfy the audited design.",
      call. = FALSE)
  list(informative_annotations = informative,
    uninformative_annotations = uninformative,
    marginal_component_probability = marginal,
    target_marginal_component_probability = target_marginal,
    expected_nonnull = sum(expected_nonnull),
    enriched_expected_nonnull_share = share)
}

annotation_design_summary <- function(A, truth, spec) {
  do.call(rbind, lapply(names(spec$scenarios), function(scenario) {
    probability <- annotation_marker_probabilities(A, truth[[scenario]],
      spec$controls$simulation$mixture_var)
    enriched <- A[, "enriched_binary"] == 1
    out <- data.frame(study_version = spec$study_version,
      scenario = scenario, marker_count = nrow(A),
      annotation_count = ncol(A), enriched_count = sum(enriched),
      enriched_prevalence = mean(enriched),
      continuous_signal_mean = mean(A[, "continuous_signal"]),
      continuous_signal_sd = stats::sd(A[, "continuous_signal"]),
      null_annotation_mean = mean(A[, "null_annotation"]),
      null_annotation_sd = stats::sd(A[, "null_annotation"]),
      expected_nonnull = sum(1 - probability[, 1L]),
      enriched_expected_nonnull_probability = mean((1 - probability[, 1L])[enriched]),
      unenriched_expected_nonnull_probability = mean((1 - probability[, 1L])[!enriched]),
      enriched_expected_nonnull_share =
        sum((1 - probability[, 1L])[enriched]) / sum(1 - probability[, 1L]),
      annotation_rank = qr(A)$rank,
      stringsAsFactors = FALSE)
    for (j in seq_len(ncol(probability)))
      out[[paste0("marginal_component_", j - 1L)]] <- mean(probability[, j])
    correlations <- stats::cor(A[, -1L, drop = FALSE])
    out$maximum_absolute_annotation_correlation <-
      max(abs(correlations[upper.tri(correlations)]))
    out$cor_enriched_continuous <- correlations["enriched_binary",
      "continuous_signal"]
    out$cor_enriched_null <- correlations["enriched_binary",
      "null_annotation"]
    out$cor_continuous_null <- correlations["continuous_signal",
      "null_annotation"]
    out
  }))
}

sample_annotation_components <- function(probability, seed) {
  set.seed(as.integer(seed))
  u <- stats::runif(nrow(probability))
  cumulative <- t(apply(probability, 1L, cumsum))
  component <- rowSums(cumulative < u) + 1L
  if (any(component < 1L | component > ncol(probability)))
    stop("Annotation component allocation failed.", call. = FALSE)
  component
}

validate_annotation_truth_realization <- function(component, A, spec) {
  simulation <- spec$controls$simulation
  counts <- tabulate(component, nbins = length(simulation$mixture_var))
  active <- sum(counts[-1L])
  if (active < simulation$nonnull_sanity_range[1L] ||
      active > simulation$nonnull_sanity_range[2L])
    stop("Realized active-marker count is outside the prespecified v2 range.",
      call. = FALSE)
  if (any(counts[-1L] < simulation$minimum_realized_active_counts))
    stop("A realized active component is below its prespecified minimum.",
      call. = FALSE)
  rows <- lapply(seq_len(length(counts) - 1L), function(stick) {
    eligible <- component >= stick
    outcome <- component > stick
    n_eligible <- sum(eligible)
    ones <- sum(outcome[eligible])
    zeros <- n_eligible - ones
    if (n_eligible < simulation$minimum_stick_eligible[stick] ||
        min(ones, zeros) < simulation$minimum_stick_binary_outcome[stick])
      stop("A probit stick has insufficient eligible markers or binary outcomes.",
        call. = FALSE)
    cell <- table(factor(A[eligible, "enriched_binary"], levels = 0:1),
      factor(outcome[eligible], levels = c(FALSE, TRUE)))
    if (any(cell == 0L))
      stop("A required enriched-annotation by stick-outcome cell is empty.",
        call. = FALSE)
    if (qr(A[eligible, , drop = FALSE])$rank != ncol(A))
      stop("The annotation design is rank deficient on a required stick subset.",
        call. = FALSE)
    data.frame(stick = stick, eligible_count = n_eligible,
      outcome_zero_count = zeros, outcome_one_count = ones,
      enriched_outcome_cell_minimum = min(cell), stringsAsFactors = FALSE)
  })
  list(component_counts = counts, stick_counts = do.call(rbind, rows))
}

simulate_annotation_architecture <- function(coordinate, scaled_genotypes,
                                              train_rows, A, truth, spec) {
  scenario <- as.character(coordinate$scenario)
  if (!scenario %in% names(spec$scenarios))
    stop("Unknown annotation scenario: ", scenario, call. = FALSE)
  if (!identical(colnames(scaled_genotypes), rownames(A)))
    stop("Simulation annotation-marker alignment failed.", call. = FALSE)
  mixture <- spec$controls$simulation$mixture_var
  probability <- annotation_marker_probabilities(A, truth[[scenario]], mixture)
  component <- sample_annotation_components(probability,
    coordinate$component_seed)
  truth_validation <- validate_annotation_truth_realization(component, A, spec)
  set.seed(as.integer(coordinate$effect_seed))
  raw_effect <- numeric(ncol(scaled_genotypes))
  active <- component > 1L
  raw_effect[active] <- stats::rnorm(sum(active),
    sd = sqrt(mixture[component[active]]))
  genetic_raw <- as.numeric(scaled_genotypes %*% raw_effect)
  h2 <- spec$controls$simulation$h2
  scale <- sqrt((h2 / (1 - h2)) /
    stats::var(genetic_raw[train_rows]))
  effects <- raw_effect * scale
  genetic <- as.numeric(scaled_genotypes %*% effects)
  set.seed(as.integer(coordinate$residual_seed))
  residual <- stats::rnorm(nrow(scaled_genotypes))
  residual <- residual - mean(residual[train_rows])
  residual <- residual / stats::sd(residual[train_rows])
  phenotype <- genetic + residual
  realized_variance <- c(raw_genetic = stats::var(genetic_raw[train_rows]),
    genetic = stats::var(genetic[train_rows]),
    residual = stats::var(residual[train_rows]),
    phenotype = stats::var(phenotype[train_rows]))
  realized_h2 <- realized_variance[["genetic"]] /
    (realized_variance[["genetic"]] + realized_variance[["residual"]])
  if (!is.finite(realized_h2) ||
      abs(realized_h2 - h2) > spec$controls$simulation$realized_h2_tolerance)
    stop("Realized phenotype heritability failed the v2 tolerance.", call. = FALSE)
  trait <- spec$data$trait
  marker_ids <- colnames(scaled_genotypes)
  sample_ids <- rownames(scaled_genotypes)
  effect_matrix <- matrix(effects, ncol = 1L,
    dimnames = list(marker_ids, trait))
  genetic_matrix <- matrix(genetic, ncol = 1L,
    dimnames = list(sample_ids, trait))
  residual_matrix <- matrix(residual, ncol = 1L,
    dimnames = list(sample_ids, trait))
  phenotype_matrix <- matrix(phenotype, ncol = 1L,
    dimnames = list(sample_ids, trait))
  causal_ids <- marker_ids[active]
  raw <- list(y = phenotype_matrix, W = scaled_genotypes, B = effect_matrix,
    G = genetic_matrix, E = residual_matrix,
    causal = list(shared = causal_ids,
      specific = stats::setNames(list(character()), trait), all = causal_ids),
    rsids = marker_ids, ids = sample_ids, h2_target = h2,
    h2_observed = realized_h2,
    shared_idx = which(active),
    specific_idx = stats::setNames(list(integer()), trait),
    causal_rsids = causal_ids)
  simulation <- as_sblrbench_simulation(raw, study = spec$study,
    architecture = scenario, replicate = as.integer(coordinate$replicate),
    seed = as.integer(coordinate$component_seed))
  simulation$data$train_ids <- sample_ids[train_rows]
  simulation$data$test_ids <- sample_ids[-train_rows]
  marker_truth <- data.frame(marker_id = marker_ids,
    component_index = component - 1L,
    component_variance = mixture[component], raw_effect = raw_effect,
    effect = effects, true_nonnull = active,
    enriched_binary = as.integer(A[, "enriched_binary"]),
    stringsAsFactors = FALSE)
  for (j in seq_len(ncol(probability)))
    marker_truth[[paste0("true_prior_component_", j - 1L)]] <- probability[, j]
  simulation$extras$annotation_matrix <- A
  simulation$extras$true_alpha <- truth[[scenario]]
  simulation$extras$true_marker_prior <- probability
  simulation$extras$marker_truth <- marker_truth
  simulation$extras$component_counts <- tabulate(component,
    nbins = length(mixture))
  simulation$extras$stick_counts <- truth_validation$stick_counts
  simulation$extras$effect_scale <- scale
  simulation$extras$study_version <- spec$study_version
  simulation$extras$target_variance <- c(heritability = h2,
    residual = 1, genetic = h2 / (1 - h2))
  simulation$extras$realized_variance <- realized_variance
  simulation$extras$realized_heritability <- realized_h2
  validate_sblrbench_simulation(simulation)
  simulation
}

prepare_annotation_simulations <- function(spec, profile, data,
                                           mode = c("final", "qualification")) {
  mode <- match.arg(mode)
  A <- construct_annotation_design(data$markers$marker_ids, spec)
  truth <- construct_annotation_truth(A, spec)
  coordinates <- unique(benchmark_annotation_seeds(spec, profile,
    mode = mode)[c("scenario", "replicate", "component_seed",
      "effect_seed", "residual_seed")])
  bundles <- lapply(seq_len(nrow(coordinates)), function(i) {
    coordinate <- as.list(coordinates[i, , drop = FALSE])
    simulation <- simulate_annotation_architecture(coordinate,
      data$scaled$all, data$split$train_rows, A, truth, spec)
    oracle <- check_oracle_genetic_values(simulation,
      tolerance = spec$validation$oracle_tolerance)
    if (!isTRUE(oracle$ok))
      stop("Study 06 simulation oracle failed.", call. = FALSE)
    test_simulation <- subset_sblrbench_simulation_samples(simulation,
      data$split$test_ids)
    stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
      spec$data)
    list(coordinate = coordinate, simulation = simulation,
      test_simulation = test_simulation, stats = stats, oracle = oracle,
      annotations = A, annotation_truth = truth,
      marker_truth = simulation$extras$marker_truth)
  })
  attr(bundles, "annotations") <- A
  attr(bundles, "annotation_truth") <- truth
  bundles
}
