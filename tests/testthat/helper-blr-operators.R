# ---- consolidated from tests/testthat/helper-mtblr-bed-contract.R ----
blr_bed_contract_bed_code <- function(dosage) {
  ifelse(is.na(dosage), 1L,
    ifelse(dosage == 2, 0L, ifelse(dosage == 1, 2L, 3L)))
}

blr_bed_contract_write_bed <- function(path, dosage) {
  dosage <- as.matrix(dosage)
  n <- nrow(dosage)
  bytes_per_marker <- ceiling(n / 4)
  payload <- raw(ncol(dosage) * bytes_per_marker)
  for (j in seq_len(ncol(dosage))) {
    codes <- blr_bed_contract_bed_code(dosage[, j])
    for (i in seq_len(n)) {
      byte <- (j - 1L) * bytes_per_marker + (i - 1L) %/% 4L + 1L
      shift <- 2L * ((i - 1L) %% 4L)
      payload[byte] <- as.raw(
        bitwOr(as.integer(payload[byte]), bitwShiftL(codes[i], shift))
      )
    }
  }
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  writeBin(payload, con)
  invisible(path)
}

blr_bed_contract_decode_bed <- function(bed_files, n_bed, cls, rows = NULL) {
  if (is.null(rows)) rows <- seq_len(n_bed)
  stopifnot(length(bed_files) == length(cls))
  out <- vector("list", sum(lengths(cls)))
  at <- 1L
  bytes_per_marker <- ceiling(n_bed / 4)
  dosage <- c(2, NA_real_, 1, 0)
  for (f in seq_along(bed_files)) {
    bytes <- readBin(bed_files[f], "raw", n = file.info(bed_files[f])$size)
    stopifnot(identical(bytes[1:3], as.raw(c(0x6c, 0x1b, 0x01))))
    for (column in cls[[f]]) {
      first <- 4L + (column - 1L) * bytes_per_marker
      packed <- bytes[first:(first + bytes_per_marker - 1L)]
      codes <- vapply(rows, function(row) {
        byte <- packed[(row - 1L) %/% 4L + 1L]
        bitwAnd(bitwShiftR(as.integer(byte), 2L * ((row - 1L) %% 4L)), 3L)
      }, integer(1))
      out[[at]] <- dosage[codes + 1L]
      at <- at + 1L
    }
  }
  do.call(cbind, out)
}

blr_bed_contract_transform_genotypes <- function(dosage, af, scale = TRUE) {
  dosage <- as.matrix(dosage)
  stopifnot(ncol(dosage) == length(af))
  out <- dosage
  for (j in seq_len(ncol(out))) {
    if (scale) {
      denom <- sqrt(2 * af[j] * (1 - af[j]))
      out[, j] <- (out[, j] - 2 * af[j]) / denom
      out[is.na(out[, j]), j] <- 0
    } else {
      out[is.na(out[, j]), j] <- 2 * af[j]
    }
  }
  out
}

blr_bed_contract_fixture <- function() {
  file1 <- rbind(
    c(0, 1, 2, NA), c(1, 2, 0, 1), c(2, NA, 1, 0),
    c(0, 0, 1, 2), c(1, 2, NA, 1), c(2, 1, 0, 2), c(NA, 0, 2, 1)
  )
  file2 <- rbind(
    c(2, 0, 1), c(1, NA, 2), c(0, 1, 1), c(2, 2, 0),
    c(NA, 1, 2), c(1, 0, NA), c(0, 2, 1)
  )
  paths <- c(tempfile(fileext = ".bed"), tempfile(fileext = ".bed"))
  blr_bed_contract_write_bed(paths[1], file1)
  blr_bed_contract_write_bed(paths[2], file2)
  cls <- list(c(4L, 2L, 1L), c(3L, 1L))
  rows <- c(7L, 2L, 5L, 1L, 6L)
  af <- c(.31, .43, .27, .38, .46)
  dosage <- blr_bed_contract_decode_bed(paths, 7L, cls, rows)
  list(paths = paths, n_bed = 7L, cls = cls, rows = rows, af = af,
       dosage = dosage, X = blr_bed_contract_transform_genotypes(dosage, af, TRUE))
}

blr_bed_contract_rebuild_residual <- function(Y, X, effective) {
  as.matrix(Y) - as.matrix(X) %*% as.matrix(effective)
}

blr_bed_contract_marker_score <- function(X, residual, marker, current_effect) {
  x <- X[, marker]
  drop(crossprod(x, residual) + drop(crossprod(x)) * current_effect)
}

blr_bed_contract_update_residual <- function(residual, x, delta) {
  as.matrix(residual) - tcrossprod(x, delta)
}

blr_bed_contract_genetic_values <- function(X, effective) {
  as.matrix(X) %*% as.matrix(effective)
}

blr_bed_contract_genetic_covariance <- function(X, effective) {
  U <- blr_bed_contract_genetic_values(X, effective)
  crossprod(U) / nrow(U)
}

blr_bed_contract_marker_conditional <- function(score, w, B, E, models, pi) {
  score <- as.numeric(score)
  models <- as.matrix(models)
  P <- solve(B)
  Omega <- solve(E)
  out <- lapply(seq_len(nrow(models)), function(k) {
    D <- diag(models[k, ], nrow = length(score))
    C <- P + w * D %*% Omega %*% D
    rhs <- D %*% Omega %*% score
    covariance <- solve(C)
    mean <- drop(covariance %*% rhs)
    log_weight <- log(pi[k]) - determinant(C, logarithm = TRUE)$modulus / 2 +
      drop(crossprod(rhs, mean)) / 2
    list(C = C, rhs = drop(rhs), mean = mean, covariance = covariance,
         log_weight = as.numeric(log_weight))
  })
  log_weights <- vapply(out, `[[`, numeric(1), "log_weight")
  probabilities <- exp(log_weights - max(log_weights))
  probabilities <- probabilities / sum(probabilities)
  list(models = out, log_weights = log_weights, probabilities = probabilities)
}

blr_bed_contract_memory_formula <- function(n, m, nt, nmodels = 2^nt, retained = 1L) {
  bytes_per_marker <- ceiling(n / 4)
  stride <- 64 * ceiling(bytes_per_marker / 64)
  list(
    packed_owner = m * stride,
    phenotype = 8 * n * nt,
    sample_residual = 8 * n * nt,
    effective_effect = 8 * m * nt,
    latent_effect = 8 * m * nt,
    state = 4 * m * nt,
    decoded_marker = 8 * n,
    marker_map = 5 * 8 * m,
    covariance_work = 8 * nt * nt * 6,
    model_work = 8 * nmodels * (nt * nt + 2 * nt + 2),
    trace_minimum = 8 * retained * (3 * nt + nmodels)
  )
}

blr_bed_contract_contract <- list(
  implementation_status = "audit_only_no_sampler",
  owner = "PackedBedMatrix",
  view = "BedPackedGenotypeView",
  phenotype = "complete_finite_centered_pre_adjusted_same_rows",
  phenotype_scaling = "not_performed",
  missing_phenotypes = "unsupported",
  covariates = "pre_adjusted_no_native_argument",
  genotype_scale = "standardized_only",
  residual_layout = "arma_mat_n_by_nt_column_major",
  residual_covariance = "full_canonical_diagonal_reduction",
  marker_decode = "one_reusable_double_workspace",
  marker_order = "summary_mt_marginal_score_stable",
  sets = "explicit_disjoint_complete",
  cpo = "unsupported",
  le_ld = "unsupported_initially",
  raw_schema = "mtblr_raw_version_1",
  raw_backend = "mt_bed_bayesc",
  raw_data_level = "individual",
  marker_wy = "X_transpose_Y",
  marker_r = "X_transpose_R_final",
  sample_outputs = "internal_only",
  execution = "serial_one_chain_fit_local_mt19937"
)

# ---- consolidated from tests/testthat/helper-mtblr-bed-internal.R ----
blr_bed_internal_prior_list <- function(x) {
  lapply(seq_len(nrow(x)), function(i) as.numeric(x[i, ]))
}

blr_bed_internal_models <- function(nt) {
  grid <- as.matrix(expand.grid(rep(list(0:1), nt)))
  lapply(seq_len(nrow(grid)), function(i) as.integer(grid[i, ]))
}

blr_bed_internal_case <- function(nt = 2L, residual_covariance = "diagonal",
                          updates = FALSE, multiple_sets = FALSE,
                          nonzero = FALSE) {
  fixture <- blr_bed_contract_fixture()
  n <- nrow(fixture$X)
  m <- ncol(fixture$X)
  Y <- vapply(seq_len(nt), function(trait) {
    value <- sin(seq_len(n) * (.41 + trait / 13)) +
      cos(seq_len(n) * (.23 + trait / 17))
    value - mean(value)
  }, numeric(n))
  colnames(Y) <- paste0("T", seq_len(nt))
  models <- blr_bed_internal_models(nt)
  pi <- seq_along(models)
  pi <- pi / sum(pi)
  B <- diag(seq(.55, .85, length.out = nt), nrow = nt)
  if (nt > 1L) B[upper.tri(B)] <- B[lower.tri(B)] <- .08
  E <- diag(seq(.8, 1.1, length.out = nt), nrow = nt)
  if (residual_covariance == "full" && nt > 1L) {
    E[upper.tri(E)] <- E[lower.tri(E)] <- .12
  }
  ssb <- diag(.25, nt)
  sse <- diag(.35, nt)
  if (residual_covariance == "full" && nt > 1L) {
    sse[upper.tri(sse)] <- sse[lower.tri(sse)] <- .04
  }
  beta <- replicate(nt, rep(0, m), simplify = FALSE)
  effective <- replicate(nt, rep(0, m), simplify = FALSE)
  state <- replicate(nt, rep(0L, m), simplify = FALSE)
  if (nonzero) {
    state[[1L]][1L] <- 1L
    beta[[1L]][1L] <- effective[[1L]][1L] <- .07
    if (nt > 1L) beta[[2L]][2L] <- -.03
  }
  sets <- if (multiple_sets) {
    list(as.integer(seq(0L, m - 1L, by = 2L)),
         as.integer(seq(1L, m - 1L, by = 2L)))
  } else {
    list(0:(m - 1L))
  }
  list(
    fixture = fixture, Y = Y, beta = beta, b = effective, state = state,
    sets = sets, B = B, E = E, ssb = ssb, sse = sse,
    models = models, pi = pi, nub = max(4, nt + 2),
    nue = max(4, nt + 2), updateB = updates, updateE = updates,
    updatePi = updates, residual_covariance = residual_covariance,
    nit = 5L, nburn = 2L, nthin = 1L, seed = 17015L, method = 4L
  )
}

blr_bed_internal_args <- function(case) {
  f <- case$fixture
  list(
    bed_files = f$paths, n_bed = f$n_bed, cls = f$cls, rows = f$rows,
    af = f$af, Y = case$Y, beta_init = case$beta, b_init = case$b,
    state_init = case$state, sets = case$sets, B = case$B, E = case$E,
    ssb_prior = blr_bed_internal_prior_list(case$ssb),
    sse_prior = blr_bed_internal_prior_list(case$sse), models = case$models,
    pi = case$pi, nub = case$nub, nue = case$nue,
    updateB = case$updateB, updateE = case$updateE,
    updatePi = case$updatePi,
    residual_covariance = case$residual_covariance,
    nit = case$nit, nburn = case$nburn, nthin = case$nthin,
    seed = case$seed, method = case$method
  )
}

blr_bed_internal_call <- function(case) {
  do.call(sblr:::mtblr_bed_internal, blr_bed_internal_args(case))
}

blr_bed_internal_dense_args <- function(case) {
  X <- case$fixture$X
  Y <- case$Y
  XX <- crossprod(X)
  nt <- ncol(Y)
  m <- ncol(X)
  indices <- replicate(m, 0:(m - 1L), simplify = FALSE)
  values <- lapply(seq_len(nt), function(trait) {
    lapply(seq_len(m), function(marker) as.numeric(XX[marker, ]))
  })
  list(
    wy = lapply(seq_len(nt), function(t) drop(crossprod(X, Y[, t]))),
    ww = replicate(nt, diag(XX), simplify = FALSE),
    yy = colSums(Y^2), b = case$b, XXvalues = values,
    XXindices = indices, sets = case$sets, B = case$B, E = case$E,
    ssb_prior = blr_bed_internal_prior_list(case$ssb),
    sse_prior = blr_bed_internal_prior_list(case$sse), models = case$models,
    pi = case$pi, nub = case$nub, nue = case$nue,
    updateB = case$updateB, updateE = case$updateE,
    updatePi = case$updatePi, n = rep.int(nrow(X), nt),
    nit = case$nit, nburn = case$nburn, nthin = case$nthin,
    seed = case$seed, method = case$method
  )
}

blr_bed_internal_legacy_matrix <- function(legacy, field) {
  do.call(cbind, legacy[[field]])
}

blr_bed_internal_compare_dense <- function(case, tolerance = 1e-10) {
  bed <- blr_bed_internal_call(case)
  dense <- do.call(sblr:::mtblr, blr_bed_internal_dense_args(case))
  fields <- list(
    bm = 1L, dm = 2L, wy = 3L, r = 4L, b = 5L, state = 6L,
    vbs = 8L, vgs = 9L, ves = 10L, covb = 11L, covg = 12L,
    cove = 13L, vb = 14L, vg = 15L, ve = 16L
  )
  actual <- lapply(names(fields), function(name) {
    bed[[if (name %in% c("bm", "dm", "wy", "r", "b", "state"))
      "marker" else if (name %in% c("vbs", "vgs", "ves"))
      "trace" else "variance"]][[name]]
  })
  names(actual) <- names(fields)
  actual$order <- bed$marker$order
  actual$pi <- bed$pi
  expected <- lapply(fields, function(field) {
    blr_bed_internal_legacy_matrix(dense, field)
  })
  storage.mode(expected$state) <- "integer"
  expected$order <- as.integer(dense[[7L]][[1L]])
  expected$pi <- list(final = dense[[17L]][[1L]],
                      mean = dense[[18L]][[1L]])
  testthat::expect_equal(actual, expected, tolerance = tolerance)
  testthat::expect_identical(actual$state, expected$state)
  testthat::expect_identical(actual$order, expected$order)
  invisible(bed)
}

blr_bed_internal_cleanup <- function(case) {
  unlink(case$fixture$paths)
  invisible()
}

# ---- consolidated from tests/testthat/helper-mtblr-bed-public.R ----
blr_bed_public_case <- function(nt = 2L, use_all_rows = FALSE,
                          matched_ids = FALSE, uncentered = FALSE) {
  fixture <- blr_bed_contract_fixture()
  full_ids <- c(paste0("a", 1:4), paste0("b", 1:3))
  Glist <- list(
    n = fixture$n_bed, ids = paste0("id", seq_len(fixture$n_bed)),
    bedfiles = fixture$paths,
    rsids = list(full_ids[1:4], full_ids[5:7]),
    rsidsLD = list(c("a4", "a2", "a1"), c("b3", "b1")),
    af = list(c(.27, .43, .20, .31), c(.46, .20, .38))
  )
  rows <- if (use_all_rows) NULL else fixture$rows
  n <- if (use_all_rows) fixture$n_bed else length(fixture$rows)
  Y <- vapply(seq_len(nt), function(trait) {
    x <- sin(seq_len(n) * (.31 + trait / 11)) +
      cos(seq_len(n) * (.19 + trait / 17))
    if (!uncentered) x <- x - mean(x)
    x
  }, numeric(n))
  colnames(Y) <- paste0("T", seq_len(nt))
  if (matched_ids) {
    stopifnot(!use_all_rows)
    rownames(Y) <- Glist$ids[fixture$rows]
    rows <- NULL
  }
  list(fixture = fixture, Glist = Glist, Y = Y, rows = rows)
}

blr_bed_public_cleanup <- function(case) {
  unlink(case$fixture$paths)
  invisible()
}

blr_bed_public_or <- function(x, y) if (is.null(x)) y else x

blr_bed_public_public_args <- function(case, residual_covariance = "full",
                                  updates = FALSE, center = FALSE, ...) {
  args <- list(
    y = case$Y, Glist = case$Glist, rows = case$rows,
    residual_covariance = residual_covariance, center = center,
    updateB = updates, updateE = updates, updatePi = updates,
    nit = 5L, nburn = 2L, nthin = 1L, seed = 17016L,
    memory_warning_gb = Inf)
  extra <- list(...)
  args[names(extra)] <- extra
  args
}

blr_bed_public_native_args <- function(public_args) {
  y <- public_args$y
  dat <- sblr:::.make_bed_marker_data(
    public_args$Glist, y, blr_bed_public_or(public_args$chr, NULL),
    blr_bed_public_or(public_args$cls, NULL),
    blr_bed_public_or(public_args$block_size, 1000L),
    blr_bed_public_or(public_args$rows, NULL))
  Y <- as.matrix(dat$y)
  if (isTRUE(blr_bed_public_or(public_args$center, TRUE))) {
    Y <- sweep(Y, 2L, colMeans(Y), "-")
  }
  nt <- ncol(Y); m <- dat$m
  mod <- sblr:::.mtblr_models(
    blr_bed_public_or(public_args$models, NULL),
    blr_bed_public_or(public_args$pimodels, NULL),
    blr_bed_public_or(public_args$pi, .001), nt)
  null <- which(rowSums(mod$matrix) == 0L)
  p_active <- 1 - sum(mod$probabilities[null])
  labels <- unique(dat$sets)
  defaults <- lapply(labels, function(label) which(dat$sets == label))
  sets <- sblr:::.mtblr_sets(blr_bed_public_or(public_args$sets, defaults), m)
  h2 <- blr_bed_public_or(public_args$h2, .5)
  if (length(h2) == 1L) h2 <- rep(h2, nt)
  vy <- colSums(Y^2) / (nrow(Y) - 1)
  vg <- blr_bed_public_or(public_args$vg, diag(vy * h2, nt))
  ve <- blr_bed_public_or(public_args$ve, diag(vy * (1 - h2), nt))
  vb <- blr_bed_public_or(public_args$vb, vg / (m * p_active))
  nub <- blr_bed_public_or(public_args$nub, 4)
  nue <- blr_bed_public_or(public_args$nue, 4)
  ssb <- blr_bed_public_or(public_args$ssb_prior,
                     ((nub - 2) / nub) * vg / (m * p_active))
  sse <- blr_bed_public_or(public_args$sse_prior,
                     ((nue - 2) / nue) * ve)
  init <- sblr:::.mtblr_bed_initialization(
    blr_bed_public_or(public_args$beta, NULL),
    blr_bed_public_or(public_args$b, NULL),
    blr_bed_public_or(public_args$state, NULL), mod$matrix, m, nt)
  list(
    bed_files = dat$bed_files, n_bed = dat$n_total, cls = dat$cls,
    rows = dat$rows, af = unlist(dat$af, use.names = FALSE), Y = Y,
    beta_init = lapply(seq_len(nt), function(t) init$beta[, t]),
    b_init = lapply(seq_len(nt), function(t) init$b[, t]),
    state_init = lapply(seq_len(nt), function(t) init$state[, t]),
    sets = sets$native, B = vb, E = ve,
    ssb_prior = lapply(seq_len(nt), function(t) ssb[t, ]),
    sse_prior = lapply(seq_len(nt), function(t) sse[t, ]),
    models = mod$native, pi = mod$probabilities, nub = nub, nue = nue,
    updateB = blr_bed_public_or(public_args$updateB, TRUE),
    updateE = blr_bed_public_or(public_args$updateE, TRUE),
    updatePi = blr_bed_public_or(public_args$updatePi, TRUE),
    residual_covariance = blr_bed_public_or(public_args$residual_covariance, "full"),
    nit = blr_bed_public_or(public_args$nit, 1000L),
    nburn = blr_bed_public_or(public_args$nburn, 500L),
    nthin = blr_bed_public_or(public_args$nthin, 1L),
    seed = blr_bed_public_or(public_args$seed, 1L),
    method = 4L)
}

blr_bed_public_compare_public_internal <- function(args, tolerance = 1e-12) {
  fit <- do.call(mtblr_bed, args)
  raw <- do.call(sblr:::mtblr_bed_internal, blr_bed_public_native_args(args))
  marker_fields <- c("bm", "dm", "wy", "r", "b")
  trace_fields <- c("vbs", "vgs", "ves")
  variance_fields <- c(
    covb = "cov_b_mean", covg = "cov_g_mean", cove = "cov_e_mean",
    vb = "cov_b_final", vg = "cov_g_final", ve = "cov_e_final")
  actual <- c(
    lapply(marker_fields, function(field) unname(fit[[field]])),
    lapply(trace_fields, function(field) unname(fit[[field]])),
    lapply(unname(variance_fields), function(field) unname(fit[[field]])),
    list(pi = unname(fit$pi_final), pim = unname(fit$pi_mean)))
  names(actual) <- c(marker_fields, trace_fields, names(variance_fields),
                     "pi", "pim")
  expected <- c(
    raw$marker[marker_fields], raw$trace[trace_fields],
    raw$variance[names(variance_fields)],
    list(pi = raw$pi$final, pim = raw$pi$mean))
  testthat::expect_equal(actual, expected, tolerance = tolerance)
  testthat::expect_identical(unname(fit$d), raw$marker$state)
  testthat::expect_identical(fit$marker_order, raw$marker$order)
  testthat::expect_equal(
    list(fit$raw_schema_version, fit$input$backend, fit$input$data_level),
    list(1L, "mt_bed_bayesc", "individual"))
  serial_diagnostics <- names(raw$diagnostics$mt_bed)
  testthat::expect_equal(fit$diagnostics$mt_bed[serial_diagnostics],
                         raw$diagnostics$mt_bed)
  invisible(fit)
}

# ---- consolidated from tests/testthat/helper-mtblr-bed-chains-internal.R ----
blr_bed_chains_seed_oracle <- function(seed, nchains) {
  base <- as.double(seed) %% 2^32
  (base + 9176 * (seq_len(nchains) - 1)) %% 2^32
}

blr_bed_chains_args <- function(case, nchains = 1L, ncores = 1L,
                          chain_seeds = integer(), keep_chains = FALSE) {
  c(blr_bed_internal_args(case), list(
    nchains = as.integer(nchains), ncores = as.integer(ncores),
    chain_seeds = as.integer(chain_seeds), keep_chains = keep_chains,
    joint_component = integer(), joint_multiplier = numeric(),
    joint_names = character(), component_count = 0L,
    marker_scale = numeric(), pi_prior = numeric(),
    component_init = integer(),
    annotations = matrix(numeric(), 0L, 0L),
    alpha_init = matrix(numeric(), 0L, 0L),
    sigma_alpha_init = numeric(), pattern_pi_init = numeric(),
    pattern_pi_prior = numeric(), updateAlpha = FALSE,
    intercept_prior_resolved = matrix(numeric(), 0L, 0L),
    sigma_alpha_a = 2, sigma_alpha_b = 2,
    pi_floor = 1e-12, alpha_update_every = 1L
  ))
}

blr_bed_chains_call <- function(case, nchains = 1L, ncores = 1L,
                          chain_seeds = integer(), keep_chains = FALSE) {
  do.call(sblr:::mtblr_bed_chains_internal,
          blr_bed_chains_args(case, nchains, ncores, chain_seeds, keep_chains))
}

blr_bed_chains_numerical_raw <- function(raw) {
  list(
    marker = raw$marker[c("bm", "dm", "wy", "r", "b", "state", "order")],
    trace = raw$trace,
    variance = raw$variance,
    pi = raw$pi,
    diagnostics = raw$diagnostics[c("marker", "covb", "covg", "cove", "pi")]
  )
}

blr_bed_chains_without_timing <- function(raw) {
  bed <- raw$diagnostics$mt_bed
  bed[c("requested_cores", "used_workers", "chain_seconds", "seconds_mean",
        "seconds_max", "dispatch_seconds")] <- NULL
  raw$diagnostics$mt_bed <- bed
  if (!is.null(raw$chains)) {
    raw$chains <- lapply(raw$chains, function(chain) {
      chain$diagnostics$seconds <- NULL
      chain
    })
  }
  raw
}

# ---- consolidated from tests/testthat/helper-mtblr-bed-public-chains.R ----
blr_bed_public_chains_public_args <- function(case, nchains = 1L, ncores = 1L,
                                 chain_seeds = NULL, keep_chains = FALSE,
                                 residual_covariance = "full",
                                 updates = FALSE,
                                 convergence = "none", ...) {
  args <- blr_bed_public_public_args(
    case, residual_covariance, updates, center = FALSE,
    nchains = nchains, ncores = ncores, chain_seeds = chain_seeds,
    keep_chains = keep_chains, convergence = convergence, ...)
  args
}

blr_bed_public_chains_internal_args <- function(public_args) {
  args <- blr_bed_public_native_args(public_args)
  args$nchains <- as.integer(blr_bed_public_or(public_args$nchains, 1L))
  args$ncores <- as.integer(blr_bed_public_or(public_args$ncores, 1L))
  seeds <- blr_bed_public_or(public_args$chain_seeds, NULL)
  args$chain_seeds <- if (is.null(seeds)) integer() else seeds
  args$keep_chains <- blr_bed_public_or(public_args$keep_chains, FALSE)
  args$joint_component <- integer()
  args$joint_multiplier <- numeric()
  args$joint_names <- character()
  args$component_count <- 0L
  args$marker_scale <- numeric()
  args$pi_prior <- numeric()
  args$component_init <- integer()
  args$annotations <- matrix(numeric(), 0L, 0L)
  args$alpha_init <- matrix(numeric(), 0L, 0L)
  args$sigma_alpha_init <- numeric()
  args$pattern_pi_init <- numeric()
  args$pattern_pi_prior <- numeric()
  args$updateAlpha <- FALSE
  args$intercept_prior_resolved <- matrix(numeric(), 0L, 0L)
  args$sigma_alpha_a <- 2
  args$sigma_alpha_b <- 2
  args$pi_floor <- 1e-12
  args$alpha_update_every <- 1L
  args
}

blr_bed_public_chains_internal <- function(public_args) {
  do.call(sblr:::mtblr_bed_chains_internal,
          blr_bed_public_chains_internal_args(public_args))
}

blr_bed_public_chains_fit_numerics <- function(fit) {
  mapping <- c(
    bm = "bm", dm = "dm", wy = "wy", r = "r", b = "b", d = "d",
    marker_order = "marker_order", vbs = "vbs", vgs = "vgs", ves = "ves",
    covb = "cov_b_mean", covg = "cov_g_mean", cove = "cov_e_mean",
    vb = "cov_b_final", vg = "cov_g_final", ve = "cov_e_final",
    pi = "pi_final", pim = "pi_mean",
    bm_sd = "bm_chain_mean_sd", bm_min = "bm_chain_mean_min",
    bm_max = "bm_chain_mean_max", dm_sd = "dm_chain_mean_sd",
    dm_min = "dm_chain_mean_min", dm_max = "dm_chain_mean_max")
  setNames(lapply(unname(mapping), function(field) unname(fit[[field]])),
           names(mapping))
}

blr_bed_public_chains_raw_numerics <- function(raw) {
  values <- list(
    bm = raw$marker$bm, dm = raw$marker$dm, wy = raw$marker$wy,
    r = raw$marker$r, b = raw$marker$b, d = raw$marker$state,
    marker_order = raw$marker$order,
    vbs = raw$trace$vbs, vgs = raw$trace$vgs, ves = raw$trace$ves,
    covb = raw$variance$covb, covg = raw$variance$covg,
    cove = raw$variance$cove, vb = raw$variance$vb,
    vg = raw$variance$vg, ve = raw$variance$ve,
    pi = raw$pi$final, pim = raw$pi$mean,
    bm_sd = raw$marker$bm_sd, bm_min = raw$marker$bm_min,
    bm_max = raw$marker$bm_max, dm_sd = raw$marker$dm_sd,
    dm_min = raw$marker$dm_min, dm_max = raw$marker$dm_max)
  lapply(values, unname)
}

blr_bed_public_chains_compare_public_internal <- function(args, tolerance = 1e-12) {
  fit <- do.call(mtblr_bed, args)
  raw <- blr_bed_public_chains_internal(args)
  testthat::expect_equal(blr_bed_public_chains_fit_numerics(fit),
                         blr_bed_public_chains_raw_numerics(raw), tolerance = tolerance)
  testthat::expect_identical(fit$nchains, raw$meta$nchains)
  testthat::expect_identical(fit$chain_seeds,
                             raw$diagnostics$mt_bed$chain_seeds)
  testthat::expect_identical(fit$chain_diagnostics$used_workers,
                             raw$diagnostics$mt_bed$used_workers)
  testthat::expect_identical(is.null(fit$chains), is.null(raw$chains))
  if (!is.null(fit$chains)) {
    testthat::expect_identical(names(fit$chains), names(raw$chains))
  }
  invisible(list(fit = fit, raw = raw))
}

blr_bed_public_chains_without_timing <- function(fit) {
  fit$diagnostics$mt_bed[c("chain_seconds", "seconds_mean", "seconds_max",
                           "dispatch_seconds", "requested_cores",
                           "used_workers")] <- NULL
  fit$chain_diagnostics[c("chain_seconds", "seconds_mean", "seconds_max",
                          "dispatch_seconds", "requested_cores",
                          "used_workers")] <- NULL
  if (!is.null(fit$chains)) {
    fit$chains <- lapply(fit$chains, function(chain) {
      chain$diagnostics$seconds <- NULL
      chain
    })
  }
  fit$input$used_workers <- NULL
  fit$input$ncores <- NULL
  fit$input$ncores_requested <- NULL
  fit$input$memory_estimate <- NULL
  fit$input$memory_warning_gb <- NULL
  fit$memory_estimate <- NULL
  fit
}

blr_bed_public_chains_preexisting_fields <- function(fit) {
  fields <- c("bm", "dm", "wy", "r", "b", "d", "marker_order",
              "vbs", "vgs", "ves", "covb", "covg", "cove",
              "vb", "vg", "ve", "pi", "pim", "rb", "rg", "re",
              "bed_diagnostics", "phenotype_preprocessing")
  fit[fields]
}

