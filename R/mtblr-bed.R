.mtblr_bed_center_tolerance <- function(y) {
  1e-10 * pmax(1, sqrt(colMeans(y^2)))
}

.mtblr_bed_matrix <- function(x, m, nt, name, mode = c("numeric", "state")) {
  mode <- match.arg(mode)
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) != nt || any(lengths(x) != m)) {
      stop(name, " must contain one length-m vector per trait.", call. = FALSE)
    }
    x <- do.call(cbind, x)
  }
  x <- as.matrix(x)
  if (!identical(dim(x), c(m, nt))) {
    stop(name, " must be an m by nt matrix or trait list.", call. = FALSE)
  }
  if (mode == "numeric") {
    storage.mode(x) <- "double"
    if (any(!is.finite(x))) stop(name, " must be finite.", call. = FALSE)
  } else {
    if (anyNA(x) || any(!x %in% 0:1)) {
      stop("state must contain only binary values.", call. = FALSE)
    }
    storage.mode(x) <- "integer"
  }
  x
}

.mtblr_bed_initialization <- function(beta, b, state, models, m, nt) {
  supplied <- c(beta = !is.null(beta), b = !is.null(b),
                state = !is.null(state))
  if (!is.null(beta)) beta <- .mtblr_bed_matrix(beta, m, nt, "beta")
  if (!is.null(b)) b <- .mtblr_bed_matrix(b, m, nt, "b")
  if (!is.null(state)) {
    state <- .mtblr_bed_matrix(state, m, nt, "state", "state")
  }
  if (is.null(b)) b <- matrix(0, m, nt)
  if (is.null(state)) state <- matrix(as.integer(b != 0), m, nt)
  if (is.null(beta)) beta <- b
  pattern_key <- function(x) apply(x, 1L, paste, collapse = "_")
  if (any(!pattern_key(state) %in% pattern_key(models))) {
    stop("Every state row must equal one supplied model pattern.",
         call. = FALSE)
  }
  if (any(b[state == 0L] != 0)) {
    stop("Effective effects must be exactly zero for inactive states.",
         call. = FALSE)
  }
  active <- state == 1L
  if (any(active) &&
      !isTRUE(all.equal(unname(b[active]), unname(beta[active]),
                       tolerance = 1e-12))) {
    stop("Active effective effects must equal their latent effects.",
         call. = FALSE)
  }
  policy <- if (!any(supplied)) "all_zero_defaults" else
    paste0(names(supplied)[supplied], collapse = "_")
  list(beta = beta, b = b, state = state, policy = policy)
}

.mtblr_bed_memory_estimate <- function(
  n, m, nt, nmodels, trace_length,
  nchains = 1L, ncores = 1L, keep_chains = FALSE,
  used_workers = NULL, convergence_memory = NULL
) {
  bytes_per_marker <- ceiling(n / 4)
  stride <- 64 * ceiling(bytes_per_marker / 64)
  legacy_components <- c(
    packed_genotype = m * stride,
    phenotype = 8 * n * nt,
    sample_residual = 8 * n * nt,
    effective_effects = 8 * m * nt,
    latent_effects = 8 * m * nt,
    state = 4 * m * nt,
    decoded_marker_workspace = 8 * n,
    marker_maps = 5 * 8 * m,
    marker_order_and_sets = 8 * m,
    covariance_work = 6 * 8 * nt * nt,
    model_work = 8 * nmodels * (nt * nt + 2 * nt + 2),
    marker_wy = 8 * m * nt,
    final_marker_r = 8 * m * nt,
    traces = 8 * trace_length * 5 * nt
  )
  shared_names <- c("packed_genotype", "phenotype", "marker_maps",
                    "marker_order_and_sets", "marker_wy")
  private_names <- c("sample_residual", "decoded_marker_workspace",
                     "covariance_work", "model_work")
  result_names <- c("effective_effects", "latent_effects", "state", "traces")
  pooled_names <- "final_marker_r"
  shared_components <- legacy_components[shared_names]
  private_components <- legacy_components[private_names]
  result_components <- legacy_components[result_names]
  pooled_components <- legacy_components[pooled_names]
  retained_components <- c(
    marker_bm_dm_b = 3 * 8 * m * nt,
    marker_state = 4 * m * nt,
    traces = 8 * trace_length * 5 * nt,
    covariance_matrices = 6 * 8 * nt * nt,
    pi_final_and_mean = 2 * 8 * nmodels,
    diagnostics = 8 * 8
  )
  requested_workers <- min(ncores, nchains)
  if (is.null(used_workers)) used_workers <- requested_workers
  shared_bytes <- sum(shared_components)
  private_bytes <- sum(private_components)
  result_bytes <- sum(result_components)
  retained_bytes <- sum(retained_components)
  pooled_bytes <- sum(pooled_components)
  retained_total <- if (keep_chains) nchains * retained_bytes else 0
  convergence_total <- if (is.null(convergence_memory)) 0 else
    convergence_memory$estimated_total_bytes
  total <- shared_bytes + requested_workers * private_bytes +
    nchains * result_bytes + retained_total + pooled_bytes
  execution_total <- shared_bytes + used_workers * private_bytes +
    nchains * result_bytes + retained_total + pooled_bytes
  total <- total + convergence_total
  execution_total <- execution_total + convergence_total
  components <- if (nchains == 1L && requested_workers == 1L &&
                    !keep_chains) legacy_components else c(
    shared_components,
    private_worker_copies = requested_workers * private_bytes,
    chain_result_copies = nchains * result_bytes,
    retained_chain_copies = retained_total,
    pooled_output = pooled_bytes,
    convergence = convergence_total
  )
  if (convergence_total > 0 && !"convergence" %in% names(components)) {
    components <- c(components, convergence = convergence_total)
  }
  convergence_components <- if (is.null(convergence_memory)) {
    c(trace_capture = 0, workspace = 0, summary_output = 0,
      retained_traces = 0)
  } else c(
    trace_capture = convergence_memory$trace_capture_bytes,
    workspace = convergence_memory$maximum_workspace_bytes,
    summary_output = convergence_memory$summary_output_bytes,
    retained_traces = convergence_memory$retained_trace_bytes)
  list(
    label = "analytical working-memory estimate",
    estimate_kind = "analytical upper-bound estimate",
    measured_rss = FALSE,
    measured_peak_rss = FALSE,
    components_bytes = components,
    shared_components_bytes = shared_components,
    shared_bytes = shared_bytes,
    private_components_bytes = private_components,
    private_state_bytes_per_worker = private_bytes,
    result_components_bytes = result_components,
    result_bytes_per_chain = result_bytes,
    retained_components_bytes = retained_components,
    retained_chain_bytes_per_chain = retained_bytes,
    pooled_output_bytes = pooled_bytes,
    nchains = as.integer(nchains),
    requested_cores = as.integer(ncores),
    requested_worker_count = as.integer(requested_workers),
    used_workers = as.integer(used_workers),
    keep_chains = keep_chains,
    estimated_concurrent_bytes = shared_bytes + requested_workers * private_bytes,
    estimated_chain_results_bytes = nchains * result_bytes,
    estimated_retained_output_bytes = retained_total,
    estimated_total_bytes = total,
    estimated_total_gib = total / 1024^3,
    execution_estimated_concurrent_bytes = shared_bytes + used_workers * private_bytes,
    execution_estimated_total_bytes = execution_total,
    execution_estimated_total_gib = execution_total / 1024^3,
    convergence_requested = if (is.null(convergence_memory)) FALSE else
      isTRUE(convergence_memory$requested),
    convergence_trace_capture = if (is.null(convergence_memory)) FALSE else
      isTRUE(convergence_memory$trace_capture),
    convergence_keep_traces = if (is.null(convergence_memory)) FALSE else
      isTRUE(convergence_memory$keep_traces),
    convergence_components_bytes = convergence_components,
    convergence_trace_capture_bytes = unname(convergence_components[1L]),
    convergence_workspace_bytes = unname(convergence_components[2L]),
    convergence_summary_output_bytes = unname(convergence_components[3L]),
    convergence_retained_trace_bytes = unname(convergence_components[4L]),
    convergence_estimated_total_bytes = convergence_total,
    convergence_estimated_total_gib = convergence_total / 1024^3
  )
}

.mtblr_bed_marker_metadata <- function(dat, Glist) {
  marker_ids <- as.character(dat$variable_names)
  metadata <- data.frame(
    marker_id = marker_ids,
    chromosome_or_file = rep(dat$chr, lengths(dat$cls)),
    bed_column = unlist(dat$cls, use.names = FALSE),
    allele_frequency = unlist(dat$af, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  explicit <- Glist$marker_metadata
  if (is.list(explicit) && length(explicit) >= max(dat$chr) &&
      all(vapply(dat$chr, function(cc) is.data.frame(explicit[[cc]]),
                 logical(1)))) {
    selected <- do.call(rbind, Map(function(cc, cl) {
      explicit[[cc]][cl, , drop = FALSE]
    }, dat$chr, dat$cls))
    if (!is.null(selected$marker_id) &&
        identical(as.character(selected$marker_id), marker_ids)) {
      for (field in c("effect_allele", "other_allele")) {
        if (!is.null(selected[[field]])) metadata[[field]] <- selected[[field]]
      }
    }
  }
  .mtblr_marker_metadata(marker_ids, metadata)
}

.mtblr_bed_trait_metadata <- function(metadata, trait_names, n,
                                      preprocessing,
                                      residual_covariance) {
  out <- if (is.null(metadata)) {
    data.frame(trait_id = trait_names, stringsAsFactors = FALSE)
  } else as.data.frame(metadata, stringsAsFactors = FALSE)
  if (is.null(out$trait_id)) out$trait_id <- trait_names
  if (nrow(out) != length(trait_names) ||
      !identical(as.character(out$trait_id), trait_names) ||
      anyNA(out$trait_id) || any(!nzchar(out$trait_id)) ||
      anyDuplicated(out$trait_id)) {
    stop("trait_metadata trait_id must uniquely match trait order.",
         call. = FALSE)
  }
  out$sample_size <- rep.int(n, length(trait_names))
  out$phenotype_mean_before <- preprocessing$mean_before
  out$phenotype_mean_after <- preprocessing$mean_after
  out$phenotype_variance <- preprocessing$variance_after
  out$phenotype_centering <- preprocessing$centering_status
  out$residual_covariance <- residual_covariance
  out$data_level <- "individual"
  out
}

.mtblr_bed_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.mtblr_bed_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be a positive integer-compatible scalar within the R integer range.",
         call. = FALSE)
  }
  as.integer(x)
}

.mtblr_bed_chain_controls <- function(nchains, ncores, chain_seeds,
                                      keep_chains) {
  nchains <- .mtblr_bed_positive_integer(nchains, "nchains")
  ncores <- .mtblr_bed_positive_integer(ncores, "ncores")
  keep_chains <- .mtblr_bed_logical(keep_chains, "keep_chains")
  requested <- chain_seeds
  if (is.null(chain_seeds)) {
    native <- integer()
  } else {
    if (!is.numeric(chain_seeds) || length(chain_seeds) != nchains ||
        any(!is.finite(chain_seeds)) || any(chain_seeds != floor(chain_seeds)) ||
        any(chain_seeds < -2147483648 | chain_seeds > 2147483647)) {
      stop(paste0(
        "chain_seeds must be NULL or a length-nchains integer-compatible ",
        "numeric vector within the signed 32-bit range."), call. = FALSE)
    }
    native <- unname(chain_seeds)
  }
  list(nchains = nchains, ncores = ncores, requested = requested,
       native = native, keep_chains = keep_chains)
}

.mtblr_bed_blas_environment <- function() {
  variables <- c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS",
                 "OMP_THREAD_LIMIT")
  as.list(Sys.getenv(variables, unset = NA_character_))
}

.mtblr_bed_convergence_memory <- function(convergence, controls, nchains,
                                           nit, nt) {
  estimate <- .blr_convergence_memory_estimate(
    nchains, nit, nt, keep_traces = controls$keep_traces)
  if (identical(convergence, "none")) {
    estimate$trace_capture_bytes <- 0
    estimate$maximum_workspace_bytes <- 0
    estimate$workspace_bytes_per_quantity <- 0
    estimate$summary_output_bytes <- 0
    estimate$retained_trace_bytes <- 0
  } else if (!isTRUE(controls$trace_route_required)) {
    estimate$trace_capture_bytes <- 0
    estimate$maximum_workspace_bytes <- 0
    estimate$workspace_bytes_per_quantity <- 0
    estimate$retained_trace_bytes <- 0
  }
  estimate$estimated_total_bytes <- sum(c(
    estimate$trace_capture_bytes, estimate$maximum_workspace_bytes,
    estimate$summary_output_bytes, estimate$retained_trace_bytes))
  estimate$estimated_total_gib <- estimate$estimated_total_bytes / 1024^3
  estimate$requested <- !identical(convergence, "none")
  estimate$trace_capture <- isTRUE(controls$trace_route_required)
  estimate$keep_traces <- isTRUE(controls$keep_traces)
  estimate
}

#' Fit joint multivariate BayesC, BayesR, and BayesRC models from PLINK BED
#'
#' Fits one joint individual-level multivariate BayesC likelihood using a
#' shared set of individuals and standardized genotypes from one BED-backed
#' genotype list. Unlike [stblr_bed()], which currently fits traits as separate
#' scalar likelihoods, `mtblr_bed()` uses joint inclusion patterns and a joint
#' marker-effect covariance. Full residual covariance is the default;
#' diagonal covariance is an explicit reduction model.
#'
#' Phenotypes are centered by default but never scaled. Covariates are not
#' fitted, and the complete aligned phenotype matrix must be finite. The
#' returned memory estimate is analytical working memory, not measured RSS or
#' measured peak RSS.
#'
#' @param y Numeric phenotype vector, matrix, or numeric data frame.
#' @param Glist One BED-backed genotype list.
#' @param covar Must be `NULL`; pass pre-adjusted phenotypes when required.
#' @param chr Optional BED file/chromosome indices.
#' @param cls Optional one-based marker columns, one vector per selected file.
#' @param rows Optional one-based BED rows in phenotype order.
#' @param scale Must be `TRUE`; genotypes are standardized with supplied
#'   selected allele frequencies.
#' @param center Center aligned phenotype columns in R. If `FALSE`, columns
#'   must already satisfy the native centering tolerance.
#' @param residual_covariance Either `"full"` or `"diagonal"`.
#' @param method One of `"bayesc"`, `"bayesr"`, or `"bayesrc"`; packed BED is individual
#'   level and therefore rejects summary-statistics `s` model names.
#' @param trait_metadata Optional data frame with one row per trait.
#' @param sets Optional disjoint complete list of one-based marker sets.
#' @param block_size Block size used for default sets within one BED file.
#' @param beta,b,state Optional latent effects, effective effects, and binary
#'   inclusion states as marker-by-trait matrices or trait lists.
#' @param h2 Heritability scalar or one value per trait, strictly in `(0,1)`.
#' @param pi Initial non-null probability or full pattern-probability vector.
#' @param models,pimodels Joint model patterns and probabilities.
#' @param mixture_var Fixed BayesR component-variance multipliers: one leading
#'   zero followed by unique, strictly increasing positive values.
#' @param joint_pi Optional initial probability vector over the deterministic
#'   joint pattern-by-component states.
#' @param joint_pi_prior Optional positive Dirichlet prior over joint states.
#' @param component Optional zero-based component initialization, one value per
#'   marker and consistent with `state`.
#' @param annotations Required marker-by-annotation numeric matrix or data
#'   frame for `"bayesrc"`. Explicit unique marker IDs are matched to selected
#'   BED markers.
#' @param add_intercept Add one intercept when none is supplied.
#' @param standardize_annotations Standardize eligible non-intercept columns.
#' @param center_binary_annotations Center and scale binary annotations when
#'   standardization is enabled.
#' @param alpha_init Optional processed-annotation-by-stick coefficient matrix.
#' @param sigmaSqAlpha_init Optional positive variance initialization per stick.
#' @param intercept_flat Use a flat prior for the first intercept coefficient.
#' @param sigmaSqAlpha_a,sigmaSqAlpha_b Positive annotation-variance prior
#'   hyperparameters.
#' @param pi_floor Probability floor used by probit stick-breaking.
#' @param alpha_update_every Positive iteration interval for coefficient updates.
#' @param updateAlpha Update annotation coefficients and their variances.
#' @param selection_s Optional fixed scalar MAF-S exponent for BayesR; the BED
#'   model name remains `bayesr` because the data are individual-level.
#' @param selection_maf Optional allele frequencies aligned to the selected
#'   marker order for the independent `selection_s` scale policy.
#' @param allow_reference_maf_for_selection_s Logical fallback control. Packed
#'   BED uses analysis-genotype frequencies by construction when this is NULL.
#' @param estimate_selection_s Logical; sampled MT S is currently unsupported.
#' @param selection_s_init,selection_s_prior,selection_s_proposal_sd Reserved
#'   sampled-S controls, rejected while `estimate_selection_s` is unsupported.
#' @param vg,vb,ve Initial genetic, marker-effect, and residual covariance.
#' @param ssb_prior,sse_prior Covariance prior scale matrices.
#' @param updateB,updateE,updatePi Scalar logical update controls.
#' @param nub,nue Covariance prior degrees of freedom.
#' @param nit,nburn,nthin MCMC controls applied identically to every chain.
#' @param seed Base fit-local native RNG seed when `chain_seeds` is `NULL`.
#' @param nchains Number of complete joint-MT chains.
#' @param ncores Requested package OpenMP chain workers. Used workers are capped
#'   by `nchains`; this does not control BLAS threads.
#' @param chain_seeds Optional signed 32-bit integer-compatible seed vector of
#'   length `nchains`, retained in supplied order. Negative values represent
#'   unsigned 32-bit seeds above the signed integer maximum.
#' @param keep_chains Retain compact per-chain posterior records.
#' @param convergence Convergence-diagnostic mode: `"auto"` preserves core-only
#'   automatic behavior, `"none"` disables capture, `"core"` explicitly requests
#'   the five core quantities, and `"extended"` adds applicable Tier 2 groups.
#' @param convergence_control Optional uniquely named list controlling warnings,
#'   thresholds, retention, extended groups, explicit selected-marker quantities,
#'   full probability-state opt-in, and hard diagnostic trace-memory guards.
#' @param memory_warning_gb Positive warning threshold in GiB, or `Inf`.
#' @param verbose Print resolved execution metadata.
#' @return An object of class `mtblr_fit` with BED diagnostics, phenotype
#'   preprocessing, alignment provenance, and an analytical memory estimate.
#'
#' @section Multichain execution:
#' Each chain fits one complete joint multivariate model. Chain-level OpenMP
#' dispatch is static, and `ncores` is capped by `nchains`. If OpenMP is not
#' available, a multi-core request warns once and runs serially. With
#' `chain_seeds = NULL`, chain zero uses `seed` and later chains use the native
#' modulo-2^32 `seed + 9176*c` policy; explicit signed seeds are used in their
#' supplied order.
#'
#' Posterior marker means, covariance means, and model-probability means pool
#' retained samples across chains. Traces are iterationwise chain means, while
#' final effects, states, residual scores, covariance matrices, and final model
#' probabilities come from primary chain 1. The `*_sd`, `*_min`, and `*_max`
#' fields summarize per-chain posterior means; they are not posterior standard
#' deviations, credible intervals, R-hat, ESS, or MCSE. Multiple chains alone
#' do not establish convergence.
#'
#' Compact retained chains omit shared BED data, phenotypes, marker order,
#' marker residuals, sample residuals, and genetic values. Timing is diagnostic
#' and nondeterministic. The package never changes global BLAS settings; users
#' should normally configure BLAS to one thread when running concurrent chain
#' workers to avoid oversubscription.
#'
#' @section Convergence diagnostics:
#' Core diagnostics use post-burn, unthinned per-chain `vbs`, `vgs`, `ves`,
#' `vle`, and `vld` traces. They report rank-normalized split and folded R-hat (using
#' their maximum), bulk and two-tail ESS, mean ESS, posterior SD, and mean
#' MCSE. R-hat requires at least two chains and four post-burn draws; ESS and
#' MCSE require at least six draws. Two or three chains carry an advisory that
#' four chains are generally recommended. Fixed `B` or `E` quantities are
#' marked `not_updated` rather than diagnosed.
#'
#' Diagnostics are independent of `keep_chains`. Optional convergence traces
#' contain only these post-burn trait-level traces, use no extra thinning, and increase
#' the analytical memory estimate. Threshold warnings are advisory and are
#' aggregated at most once per fit; multiple chains or passing thresholds do
#' not prove convergence. Extended mode adds applicable low-dimensional and
#' explicitly selected-marker scalar diagnostics. Matrix/simplex convergence,
#' ESS-for-SD, and quantile/median MCSE remain outside the formal engine.
#' @export
mtblr_bed <- function(
  y, Glist, covar = NULL, chr = NULL, cls = NULL, rows = NULL,
  scale = TRUE, center = TRUE,
  residual_covariance = c("full", "diagonal"), method = "bayesc",
  trait_metadata = NULL, sets = NULL, block_size = 1000,
  beta = NULL, b = NULL, state = NULL, h2 = 0.5, pi = 0.001,
  models = NULL, pimodels = NULL, mixture_var = NULL, joint_pi = NULL,
  joint_pi_prior = NULL, component = NULL,
  annotations = NULL, add_intercept = TRUE,
  standardize_annotations = TRUE, center_binary_annotations = FALSE,
  alpha_init = NULL, sigmaSqAlpha_init = NULL, intercept_flat = TRUE,
  sigmaSqAlpha_a = 2, sigmaSqAlpha_b = 2, pi_floor = 1e-12,
  alpha_update_every = 1L, updateAlpha = TRUE, selection_s = NULL,
  selection_maf = NULL, allow_reference_maf_for_selection_s = FALSE,
  estimate_selection_s = FALSE, selection_s_init = NULL,
  selection_s_prior = NULL, selection_s_proposal_sd = NULL,
  vg = NULL, vb = NULL, ve = NULL,
  ssb_prior = NULL, sse_prior = NULL,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE, nub = 4, nue = 4,
  nit = 1000, nburn = 500, nthin = 1, seed = 1,
  nchains = 1L, ncores = 1L, chain_seeds = NULL, keep_chains = FALSE,
  convergence = c("auto", "none", "core", "extended"), convergence_control = NULL,
  memory_warning_gb = 8, verbose = FALSE
) {
  semantics <- .mtblr_resolve_public_method(method, "packed_bed")
  if (missing(Glist) || is.null(Glist)) {
    stop("One BED-backed Glist is required by mtblr_bed().", call. = FALSE)
  }
  if (is.list(Glist) && is.null(Glist$bedfiles) && length(Glist) &&
      all(vapply(Glist, function(x) is.list(x) && !is.null(x$bedfiles),
                 logical(1)))) {
    stop("mtblr_bed() accepts one Glist, not one Glist per trait.",
         call. = FALSE)
  }
  if (!is.list(Glist) || is.null(Glist$bedfiles)) {
    stop("Glist must contain BED file information in Glist$bedfiles.",
         call. = FALSE)
  }
  if (!is.null(covar)) {
    stop(paste(
      "mtblr_bed() does not currently fit or project covariates.",
      "Pass phenotypes that have already been adjusted for the desired covariates."),
      call. = FALSE)
  }
  if (!isTRUE(scale) || length(scale) != 1L) {
    stop("mtblr_bed() requires scale = TRUE for standardized genotypes.",
         call. = FALSE)
  }
  center <- .mtblr_bed_logical(center, "center")
  residual_covariance <- match.arg(residual_covariance)
  convergence <- match.arg(convergence)
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size <= 0 ||
      block_size != as.integer(block_size)) {
    stop("block_size must be a positive integer-compatible scalar.",
         call. = FALSE)
  }
  if (is.data.frame(y) && !all(vapply(y, is.numeric, logical(1)))) {
    stop("A phenotype data frame must contain only numeric columns.",
         call. = FALSE)
  }
  numeric_y <- if (is.data.frame(y)) {
    all(vapply(y, is.numeric, logical(1)))
  } else {
    is.numeric(y)
  }
  if (!numeric_y || (length(y) == 0L && is.null(dim(y)))) {
    stop("y must be a nonempty numeric vector, matrix, or data frame.",
         call. = FALSE)
  }
  input_ids <- if (is.null(dim(y))) names(y) else rownames(y)
  input_sample_count <- if (is.null(dim(y))) length(y) else nrow(y)
  explicit_rows <- !is.null(rows)
  explicit_cls <- !is.null(cls)
  dat <- .make_bed_marker_data(
    Glist = Glist, y = y, chr = chr, cls = cls,
    block_size = as.integer(block_size), rows = rows)
  Y <- as.matrix(dat$y)
  if (nrow(Y) <= 1L || ncol(Y) <= 0L || any(!is.finite(Y))) {
    stop("The aligned phenotype must be a complete finite matrix with more than one row.",
         call. = FALSE)
  }
  trait_names <- colnames(Y)
  if (is.null(trait_names)) trait_names <- paste0("T", seq_len(ncol(Y)))
  if (anyNA(trait_names) || any(!nzchar(trait_names)) ||
      anyDuplicated(trait_names)) {
    stop("Trait names must be unique, nonempty, and non-missing.",
         call. = FALSE)
  }
  colnames(Y) <- trait_names
  mean_before <- colMeans(Y)
  if (center) Y <- sweep(Y, 2L, mean_before, "-")
  mean_after <- colMeans(Y)
  tolerance <- .mtblr_bed_center_tolerance(Y)
  if (any(abs(mean_after) > tolerance)) {
    stop("center = FALSE requires phenotype columns already centered to the Phase 17O tolerance.",
         call. = FALSE)
  }
  variance_after <- apply(Y, 2L, stats::var)
  if (any(!is.finite(variance_after)) || any(variance_after <= 0)) {
    stop("Every aligned phenotype trait must have positive finite variance.",
         call. = FALSE)
  }
  centering_status <- if (center) "centered_by_adapter" else
    "verified_precentered"
  preprocessing <- list(
    mean_before = unname(mean_before), mean_after = unname(mean_after),
    variance_after = unname(variance_after), center_requested = center,
    center_applied = center, centering_tolerance = unname(tolerance),
    centering_status = rep(centering_status, dat$nt),
    phenotype_scaling = "not_performed",
    missing_phenotype_policy = "complete_matrix_required")

  frequencies <- unlist(dat$af, use.names = FALSE)
  if (length(frequencies) != dat$m || any(!is.finite(frequencies)) ||
      any(frequencies <= 0 | frequencies >= 1)) {
    stop("Selected Glist allele frequencies must be finite and strictly inside (0, 1).",
         call. = FALSE)
  }
  marker_metadata <- .mtblr_bed_marker_metadata(dat, Glist)
  pattern_spec <- .mtblr_models(models, pimodels, pi, dat$nt)
  maf_info <- .mtblr_resolve_selection_maf(
    selection_maf, !is.null(selection_s) || isTRUE(estimate_selection_s),
    dat$m, analysis_frequency = frequencies,
    allow_reference_maf_for_selection_s =
      allow_reference_maf_for_selection_s)
  mixture <- .mtblr_bayesr_spec(
    semantics$prior_kernel,pattern_spec,maf_info$values,dat$m,mixture_var,joint_pi,
    joint_pi_prior,component,selection_s,estimate_selection_s,
    selection_s_init,selection_s_prior,selection_s_proposal_sd)
  if (identical(semantics$prior_kernel, "bayesrc")) {
    if (!is.null(joint_pi) || !is.null(joint_pi_prior))
      stop("joint_pi and joint_pi_prior are not BayesRC controls; use pimodels for conditional pattern initialization.", call. = FALSE)
    mixture$method_code <- 6L
  }
  bayesrc <- .mtblr_bayesrc_controls(
    semantics$prior_kernel, annotations, marker_metadata$marker_id,
    pattern_spec, mixture, add_intercept, standardize_annotations,
    center_binary_annotations, alpha_init, sigmaSqAlpha_init,
    intercept_flat, sigmaSqAlpha_a, sigmaSqAlpha_b, pi_floor,
    alpha_update_every, updateAlpha)
  model_spec <- mixture$patterns
  null_index <- which(rowSums(model_spec$matrix) == 0L)
  p_active <- 1 - sum(model_spec$probabilities[null_index])
  if (!is.finite(p_active) || p_active <= 0) {
    stop("Model probabilities must assign positive mass to non-null patterns.",
         call. = FALSE)
  }
  if (is.null(sets)) {
    labels <- unique(dat$sets)
    default_sets <- lapply(labels, function(label) which(dat$sets == label))
    set_spec <- .mtblr_sets(default_sets, dat$m)
    set_source <- if (length(dat$chr) > 1L) "chromosome_or_file" else
      "block_size"
  } else {
    set_spec <- .mtblr_sets(sets, dat$m)
    set_source <- "explicit_sets"
  }
  h2 <- as.numeric(h2)
  if (!(length(h2) %in% c(1L, dat$nt)) || any(!is.finite(h2)) ||
      any(h2 <= 0 | h2 >= 1)) {
    stop("h2 must be a finite scalar or length-nt vector in (0, 1).",
         call. = FALSE)
  }
  if (length(h2) == 1L) h2 <- rep(h2, dat$nt)
  for (name in c("nub", "nue")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= max(2, dat$nt - 1L)) {
      stop(name, " must be finite and greater than max(2, nt - 1).",
           call. = FALSE)
    }
  }
  vy <- colSums(Y^2) / (nrow(Y) - 1)
  vg0 <- diag(vy * h2, dat$nt)
  ve0 <- diag(vy * (1 - h2), dat$nt)
  vb0 <- vg0 / (dat$m * p_active)
  vg <- .mtblr_cov(vg, vg0, "vg", dat$nt)
  vb <- .mtblr_cov(vb, vb0, "vb", dat$nt)
  if (residual_covariance == "diagonal" && !is.null(ve) &&
      any(as.matrix(ve)[row(as.matrix(ve)) != col(as.matrix(ve))] != 0)) {
    stop("ve must be exactly diagonal when residual_covariance = 'diagonal'.",
         call. = FALSE)
  }
  ve <- .mtblr_cov(ve, ve0, "ve", dat$nt)
  ssb0 <- ((nub - 2) / nub) * vg / (dat$m * p_active)
  ssb_prior <- .mtblr_cov(ssb_prior, ssb0, "ssb_prior", dat$nt)
  if (residual_covariance == "diagonal" && !is.null(sse_prior) &&
      any(as.matrix(sse_prior)[row(as.matrix(sse_prior)) !=
                               col(as.matrix(sse_prior))] != 0)) {
    stop("sse_prior must be exactly diagonal when residual_covariance = 'diagonal'.",
         call. = FALSE)
  }
  sse0 <- ((nue - 2) / nue) * ve
  sse_prior <- .mtblr_cov(sse_prior, sse0, "sse_prior", dat$nt)
  initialization <- if (semantics$prior_kernel == "bayesc") {
    .mtblr_bed_initialization(beta,b,state,model_spec$matrix,dat$m,dat$nt)
  } else {
    .mtblr_bayesr_initialization(
      beta,b,state,mixture$component_init,pattern_spec$matrix,dat$m,dat$nt,
      semantics$prior_kernel)
  }
  updateB <- .mtblr_bed_logical(updateB, "updateB")
  updateE <- .mtblr_bed_logical(updateE, "updateE")
  updatePi <- .mtblr_bed_logical(updatePi, "updatePi")
  for (name in c("nit", "nburn", "nthin", "seed")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != as.integer(value) ||
        (name == "nburn" && value < 0) ||
        (name != "nburn" && value <= 0)) {
      stop(name, " must be an integer-compatible scalar in its valid range.",
           call. = FALSE)
    }
  }
  chain_control <- .mtblr_bed_chain_controls(
    nchains, ncores, chain_seeds, keep_chains)
  nchains <- chain_control$nchains
  ncores <- chain_control$ncores
  keep_chains <- chain_control$keep_chains
  native_chain_seeds <- chain_control$native
  convergence_controls <- .mtblr_bed_convergence_controls(
    convergence, convergence_control, nchains)
  extended_plan <- .blr_mtblr_extended_plan(
    convergence_controls, marker_metadata$marker_id, trait_names, method,
    mixture, bayesrc, updateB, updateE, updatePi, residual_covariance,
    nchains, as.integer(nit))
  convergence_controls$selected_markers_resolved <- extended_plan$selected
  convergence_memory <- .mtblr_bed_convergence_memory(
    convergence, convergence_controls, nchains, as.integer(nit), dat$nt)
  if (!is.numeric(memory_warning_gb) || length(memory_warning_gb) != 1L ||
      is.na(memory_warning_gb) || memory_warning_gb <= 0) {
    stop("memory_warning_gb must be a positive finite scalar or Inf.",
         call. = FALSE)
  }
  memory <- .mtblr_bed_memory_estimate(
    nrow(Y), dat$m, dat$nt, nrow(model_spec$matrix), nit + nburn,
    nchains, ncores, keep_chains,
    convergence_memory = convergence_memory)
  memory <- .mtblr_bayesr_memory(
    memory,method,dat$m,nchains,ncores,nit+nburn,nrow(model_spec$matrix),
    mixture$component_count)
  memory <- .mtblr_bayesrc_memory(
    memory,bayesrc,dat$m,nchains,ncores,nit+nburn)
  memory <- .blr_add_extended_memory(memory, extended_plan)
  if (memory$estimated_total_gib > memory_warning_gb) {
    warning(sprintf(
      paste0("mtblr_bed analytical upper-bound estimate (not measured RSS; ",
             "not measured peak RSS): n=%d, m=%d, nt=%d, models=%d, ",
             "nchains=%d, ncores=%d, requested workers=%d, keep_chains=%s, ",
             "convergence=%s, convergence trace capture=%s, ",
             "convergence trace retention=%s, convergence=%.6f GiB, ",
             "estimated %.6f GiB exceeds threshold %.6f GiB."),
      nrow(Y), dat$m, dat$nt, nrow(model_spec$matrix),
      nchains, ncores, memory$requested_worker_count, keep_chains,
      convergence, convergence_controls$trace_route_required,
      convergence_controls$keep_traces,
      memory$convergence_estimated_total_gib,
      memory$estimated_total_gib, memory_warning_gb), call. = FALSE)
  }
  normalized_bed_files <- normalizePath(
    dat$bed_files, winslash = "/", mustWork = TRUE)
  unmatched <- if (!is.null(input_ids)) {
    setdiff(as.character(input_ids), rownames(Y) %||% character())
  } else character()
  row_status <- if (explicit_rows) "explicit_rows" else if (!is.null(input_ids))
    "matched_by_id" else "all_glist_rows"
  selected_rows <- if (is.null(dat$rows)) seq_len(dat$n_total) else dat$rows
  alignment <- list(
    individual_policy = "shared_individual_level",
    input_sample_count = input_sample_count,
    selected_sample_count = dat$n_used,
    row_selection_status = row_status,
    selected_rows = selected_rows,
    sample_order_status = "phenotype_order_preserved",
    unmatched_input_ids = unmatched,
    unmatched_input_count = length(unmatched),
    duplicate_policy = "error",
    phenotype_missingness_status = "complete",
    phenotype_centering_status = centering_status,
    marker_selection_status = if (explicit_cls) "explicit_cls" else
      "default_rsidsLD",
    marker_order_status = "selected_glist_order_preserved",
    genotype_orientation_status = "by_construction_same_glist",
    genotype_scale_status = "standardized_genotype",
    covariate_policy = "pre_adjusted_or_unadjusted_as_supplied")
  trait_metadata <- .mtblr_bed_trait_metadata(
    trait_metadata, trait_names, dat$n_used, preprocessing,
    residual_covariance)
  blas_thread_environment <- .mtblr_bed_blas_environment()
  blas_policy <- "package_does_not_modify_blas_threads"

  native_arguments <- list(
    bed_files = dat$bed_files, n_bed = dat$n_total, cls = dat$cls,
    rows = dat$rows, af = frequencies, Y = Y,
    beta_init = lapply(seq_len(dat$nt), function(t) initialization$beta[, t]),
    b_init = lapply(seq_len(dat$nt), function(t) initialization$b[, t]),
    state_init = lapply(seq_len(dat$nt), function(t) initialization$state[, t]),
    sets = set_spec$native, B = vb, E = ve,
    ssb_prior = lapply(seq_len(dat$nt), function(t) ssb_prior[t, ]),
    sse_prior = lapply(seq_len(dat$nt), function(t) sse_prior[t, ]),
    models = model_spec$native, pi = model_spec$probabilities,
    nub = nub, nue = nue, updateB = updateB, updateE = updateE,
    updatePi = updatePi, residual_covariance = residual_covariance,
    nit = as.integer(nit), nburn = as.integer(nburn),
    nthin = as.integer(nthin), seed = as.integer(seed),
    method = mixture$method_code,
    nchains = nchains, ncores = ncores,
    chain_seeds = native_chain_seeds, keep_chains = keep_chains,
    joint_component = mixture$joint_component,
    joint_multiplier = mixture$joint_multiplier,
    joint_names = mixture$joint_names,
    component_count = mixture$component_count,
    marker_scale = mixture$marker_scale,
    pi_prior = mixture$pi_prior,
    component_init = initialization$component %||% integer(),
    annotations = bayesrc$annotations, alpha_init = bayesrc$alpha_init,
    sigma_alpha_init = bayesrc$sigma_alpha_init,
    pattern_pi_init = bayesrc$pattern_pi_init,
    pattern_pi_prior = bayesrc$pattern_pi_prior,
    updateAlpha = bayesrc$updateAlpha,
    intercept_flat = bayesrc$intercept_flat,
    sigma_alpha_a = bayesrc$sigma_alpha_a,
    sigma_alpha_b = bayesrc$sigma_alpha_b,
    pi_floor = bayesrc$pi_floor,
    alpha_update_every = bayesrc$alpha_update_every,
    convergence_covariance = extended_plan$native$convergence_covariance,
    convergence_probability = extended_plan$native$convergence_probability,
    convergence_annotations = extended_plan$native$convergence_annotations,
    convergence_full_probability =
      extended_plan$native$convergence_full_probability,
    convergence_markers = extended_plan$native$convergence_markers,
    convergence_b = extended_plan$native$convergence_b,
    convergence_d = extended_plan$native$convergence_d,
    convergence_component = extended_plan$native$convergence_component)
  native_route <- if (isTRUE(convergence_controls$trace_route_required)) {
    mtblr_bed_convergence_trace_internal
  } else {
    mtblr_bed_chains_internal
  }
  native_result <- do.call(native_route, native_arguments)
  if (!is.null(bayesrc$model_parameters)) {
    if (isTRUE(convergence_controls$trace_route_required))
      native_result$raw <- .mtblr_bayesrc_enrich_raw(
        native_result$raw, bayesrc, method, updatePi)
    else native_result <- .mtblr_bayesrc_enrich_raw(
      native_result, bayesrc, method, updatePi)
  }
  if (isTRUE(convergence_controls$trace_route_required)) {
    native_result$raw <- .mtblr_bayesr_enrich_raw(
      native_result$raw, method, mixture$model_parameters)
  } else {
    native_result <- .mtblr_bayesr_enrich_raw(
      native_result, method, mixture$model_parameters)
  }
  convergence_traces <- NULL
  if (isTRUE(convergence_controls$trace_route_required)) {
    diagnostic_result <- .mtblr_bed_convergence_internal(
      native_result = native_result, trait_names = trait_names,
      updateB = updateB, updateE = updateE,
      control = convergence_controls$thresholds,
      keep_traces = convergence_controls$keep_traces,
      extended_plan = extended_plan, model = method)
    raw <- diagnostic_result$raw
    convergence_traces <- diagnostic_result$convergence_traces
  } else {
    raw <- native_result
    raw$diagnostics$convergence <- if (identical(convergence, "none")) {
      .blr_convergence_not_requested(
        trait_names, updateB, updateE, nchains, as.integer(nit),
        convergence_controls$thresholds)
    } else {
      .blr_convergence_unavailable(
        trait_names, updateB, updateE, nchains, as.integer(nit),
        convergence_controls$thresholds)
    }
  }
  convergence_result <- raw$diagnostics$convergence
  convergence_result$warning_messages <- if (identical(convergence, "none") ||
                                                (identical(convergence, "auto") &&
                                                 nchains < 2L)) {
    character()
  } else {
    .blr_convergence_warning_messages(convergence_result, convergence)
  }
  raw$diagnostics$convergence <-
    .blr_validate_convergence_result(convergence_result)
  raw <- .validate_mtblr_raw(raw)
  bed_diagnostics <- raw$diagnostics$mt_bed
  alignment$nchains <- nchains
  alignment$ncores_requested <- ncores
  alignment$used_workers <- as.integer(bed_diagnostics$used_workers)
  alignment$openmp_available <- bed_diagnostics$openmp_available
  alignment$chain_seeds <- bed_diagnostics$chain_seeds
  alignment$keep_chains <- keep_chains
  alignment$chain_topology <- "one_complete_joint_mt_model_per_chain"
  memory <- .mtblr_bed_memory_estimate(
    nrow(Y), dat$m, dat$nt, nrow(model_spec$matrix), nit + nburn,
    nchains, ncores, keep_chains,
    used_workers = as.integer(bed_diagnostics$used_workers),
    convergence_memory = convergence_memory)
  memory <- .mtblr_bayesr_memory(
    memory,method,dat$m,nchains,ncores,nit+nburn,nrow(model_spec$matrix),
    mixture$component_count)
  memory <- .mtblr_bayesrc_memory(
    memory,bayesrc,dat$m,nchains,ncores,nit+nburn)
  memory <- .blr_add_extended_memory(memory, extended_plan)
  raw$model$names <- model_spec$names
  raw$pi$names <- model_spec$names
  raw$data <- list(
    marker_metadata = marker_metadata, trait_metadata = trait_metadata,
    n_total = dat$n_total, n_used = dat$n_used, m = dat$m, nt = dat$nt,
    bed_files = normalized_bed_files, chr = dat$chr, cls = dat$cls,
    selected_rows = selected_rows, allele_frequencies = frequencies,
    genotype_scale = "standardized_genotype",
    data_level = "individual",
    selection_maf_source = maf_info$selection_maf_source,
    selection_maf_population = maf_info$selection_maf_population,
    selection_maf_alignment_status = maf_info$selection_maf_alignment_status,
    selection_maf_fallback_used = maf_info$selection_maf_fallback_used,
    annotation_source = bayesrc$metadata$annotation_source %||% "not_applicable",
    annotation_marker_alignment_status =
      bayesrc$metadata$annotation_marker_alignment_status %||% "not_applicable",
    missing_genotype_policy = "mean_imputed_after_centering",
    phenotype_centering = preprocessing,
    phenotype_units = "retained_not_scaled",
    residual_covariance = residual_covariance,
    sets = set_spec$public, set_source = set_source,
    memory_estimate = memory,
    nchains = nchains, ncores_requested = ncores,
    used_workers = as.integer(bed_diagnostics$used_workers),
    openmp_available = bed_diagnostics$openmp_available,
    chain_seeds = bed_diagnostics$chain_seeds,
    keep_chains = keep_chains, blas_policy = blas_policy,
    blas_thread_environment = blas_thread_environment,
    convergence = convergence,
    convergence_requested = convergence_controls$diagnostic_requested,
    convergence_scope = if (identical(convergence, "none")) "none" else convergence,
    convergence_keep_traces = convergence_controls$keep_traces,
    convergence_control = convergence_controls[c(
      "warn", "rhat_threshold", "ess_per_chain_threshold",
      "mcse_mean_over_sd_threshold", "keep_traces",
      "extended_groups_requested", "extended_groups_resolved",
      "selected_markers", "selected_marker_quantities",
      "selected_markers_resolved",
      "full_probability_states", "max_trace_gb", "allow_large_traces")],
    convergence_memory_estimate = convergence_memory)
  raw$alignment <- alignment
  input <- list(
    method = method, model = method,
    backend = paste0("mt_bed_", semantics$prior_kernel),
    prior_kernel = semantics$prior_kernel,
    data_level = "individual",
    effect_scale_policy = if (!is.null(selection_s))
      if (semantics$prior_kernel %in% c("bayesr", "bayesrc")) "component_maf_s" else "maf_s"
      else if (semantics$prior_kernel %in% c("bayesr", "bayesrc")) "component" else "unit",
    model_semantics_version = 2L,
    model_semantics = "s_prefix_means_summary_statistics",
    residual_covariance = residual_covariance,
    genotype_scale = "standardized_genotype",
    phenotype_centering = centering_status,
    phenotype_scaling = "not_performed",
    covariate_policy = "pre_adjusted_or_unadjusted_as_supplied",
    covariates_fitted = FALSE,
    missing_phenotype_policy = "complete_matrix_required",
    cpo = "unsupported", le_ld = "trait_diagonal_decomposition",
    sample_residual_returned = FALSE, genetic_values_returned = FALSE,
    n = dat$n_used, n_total = dat$n_total, n_used = dat$n_used,
    m = dat$m, nt = dat$nt, chr = dat$chr, cls = dat$cls,
    rows = selected_rows, block_size = as.integer(block_size),
    sets = set_spec$public, set_source = set_source, h2 = h2, vy = vy,
    vg = vg, vb = vb, ve = ve, ssb_prior = ssb_prior,
    sse_prior = sse_prior, nub = nub, nue = nue,
    models = model_spec$matrix, model_names = model_spec$names,
    pimodels = model_spec$probabilities, mixture_var = mixture$mixture_var,
    joint_pi_prior = mixture$pi_prior, selection_s = mixture$selection_s,
    estimate_selection_s = estimate_selection_s,
    selection_maf_source = maf_info$selection_maf_source,
    selection_maf_population = maf_info$selection_maf_population,
    selection_maf_alignment_status = maf_info$selection_maf_alignment_status,
    selection_maf_fallback_used = maf_info$selection_maf_fallback_used,
    initialization_policy = initialization$policy,
    updateB = updateB, updateE = updateE, updatePi = updatePi,
    updateAlpha = bayesrc$updateAlpha,
    annotation_policy = if (is.null(bayesrc$model_parameters)) "global" else
      "annotation_probit_stick",
    nit = as.integer(nit), nburn = as.integer(nburn),
    nthin = as.integer(nthin), seed = as.integer(seed),
    nchains = nchains, ncores = ncores, ncores_requested = ncores,
    used_workers = as.integer(bed_diagnostics$used_workers),
    openmp_available = bed_diagnostics$openmp_available,
    base_seed = as.integer(seed),
    chain_seeds_requested = chain_control$requested,
    chain_seeds_resolved = bed_diagnostics$chain_seeds,
    chain_seed_policy = "native_uint32_base_plus_9176_or_explicit_signed_seeds",
    keep_chains = keep_chains,
    primary_chain = 1L, final_state_policy = "primary_chain",
    posterior_summary_policy = "pooled_retained_samples",
    trace_policy = "iterationwise_chain_mean",
    blas_policy = blas_policy,
    blas_thread_environment = blas_thread_environment,
    trait_metadata = trait_metadata, marker_metadata = marker_metadata,
    alignment = alignment, memory_estimate = memory,
    memory_warning_gb = memory_warning_gb,
    convergence = convergence,
    convergence_requested = convergence_controls$diagnostic_requested,
    convergence_scope = if (identical(convergence, "none")) "none" else convergence,
    convergence_trace_route = if (convergence_controls$trace_route_required)
      "mtblr_bed_convergence_trace_internal" else
      "mtblr_bed_chains_internal",
    convergence_warning_enabled = convergence_controls$warn,
    convergence_warning_emitted = convergence_controls$warn &&
      length(raw$diagnostics$convergence$warning_messages) > 0L,
    convergence_keep_traces = convergence_controls$keep_traces,
    convergence_control = convergence_controls[c(
      "warn", "rhat_threshold", "ess_per_chain_threshold",
      "mcse_mean_over_sd_threshold", "keep_traces",
      "extended_groups_requested", "extended_groups_resolved",
      "selected_markers", "selected_marker_quantities",
      "selected_markers_resolved",
      "full_probability_states", "max_trace_gb", "allow_large_traces")],
    convergence_thresholds = convergence_controls$thresholds,
    convergence_status = raw$diagnostics$convergence$overall_status,
    convergence_computed = raw$diagnostics$convergence$computed,
    convergence_memory_estimate = convergence_memory)
  if (isTRUE(verbose)) {
    print(input[c("backend", "n_used", "m", "nt",
                  "residual_covariance", "seed", "nchains", "ncores",
                  "used_workers")])
  }
  fit <- .as_mtblr_fit(
    raw, marker_metadata$marker_id, trait_names, marker_metadata,
    trait_metadata, alignment, input)
  fit$bed_diagnostics <- raw$diagnostics$mt_bed
  fit$phenotype_preprocessing <- preprocessing
  fit$memory_estimate <- memory
  fit <- .mtblr_bayesr_format_fit(fit, mixture$model_parameters)
  fit <- .mtblr_bayesrc_format_fit(fit, raw$annotations, bayesrc)
  if (isTRUE(bayesrc$maf_annotation_overlap) && !is.null(selection_s))
    warning("MAF-derived annotations and selection_s are both active; MAF may influence component probabilities and effect-size variance.", call. = FALSE)
  fit["convergence_traces"] <- list(convergence_traces)
  if (isTRUE(input$convergence_warning_emitted)) {
    warning(raw$diagnostics$convergence$warning_messages[1L], call. = FALSE)
  }
  .blr_finalize_fit(
    fit, "mtblr", method, "packed_bed", data = raw$data,
    diagnostics = raw$diagnostics, memory_estimate = memory)
}
