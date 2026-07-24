phase17n_bed_code <- function(dosage) {
  ifelse(is.na(dosage), 1L,
    ifelse(dosage == 2, 0L, ifelse(dosage == 1, 2L, 3L)))
}

phase17n_write_bed <- function(path, dosage) {
  dosage <- as.matrix(dosage)
  n <- nrow(dosage)
  bytes_per_marker <- ceiling(n / 4)
  payload <- raw(ncol(dosage) * bytes_per_marker)
  for (j in seq_len(ncol(dosage))) {
    codes <- phase17n_bed_code(dosage[, j])
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

phase17n_decode_bed <- function(bed_files, n_bed, cls, rows = NULL) {
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

phase17n_transform_genotypes <- function(dosage, af, scale = TRUE) {
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

phase17n_fixture <- function() {
  file1 <- rbind(
    c(0, 1, 2, NA), c(1, 2, 0, 1), c(2, NA, 1, 0),
    c(0, 0, 1, 2), c(1, 2, NA, 1), c(2, 1, 0, 2), c(NA, 0, 2, 1)
  )
  file2 <- rbind(
    c(2, 0, 1), c(1, NA, 2), c(0, 1, 1), c(2, 2, 0),
    c(NA, 1, 2), c(1, 0, NA), c(0, 2, 1)
  )
  paths <- c(tempfile(fileext = ".bed"), tempfile(fileext = ".bed"))
  phase17n_write_bed(paths[1], file1)
  phase17n_write_bed(paths[2], file2)
  cls <- list(c(4L, 2L, 1L), c(3L, 1L))
  rows <- c(7L, 2L, 5L, 1L, 6L)
  af <- c(.31, .43, .27, .38, .46)
  dosage <- phase17n_decode_bed(paths, 7L, cls, rows)
  list(paths = paths, n_bed = 7L, cls = cls, rows = rows, af = af,
       dosage = dosage, X = phase17n_transform_genotypes(dosage, af, TRUE))
}

phase17n_rebuild_residual <- function(Y, X, effective) {
  as.matrix(Y) - as.matrix(X) %*% as.matrix(effective)
}

phase17n_marker_score <- function(X, residual, marker, current_effect) {
  x <- X[, marker]
  drop(crossprod(x, residual) + drop(crossprod(x)) * current_effect)
}

phase17n_update_residual <- function(residual, x, delta) {
  as.matrix(residual) - tcrossprod(x, delta)
}

phase17n_genetic_values <- function(X, effective) {
  as.matrix(X) %*% as.matrix(effective)
}

phase17n_genetic_covariance <- function(X, effective) {
  U <- phase17n_genetic_values(X, effective)
  crossprod(U) / nrow(U)
}

phase17n_marker_conditional <- function(score, w, B, E, models, pi) {
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

phase17n_memory_formula <- function(n, m, nt, nmodels = 2^nt, retained = 1L) {
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

phase17n_contract <- list(
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
