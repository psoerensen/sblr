phase17o_prior_list <- function(x) {
  lapply(seq_len(nrow(x)), function(i) as.numeric(x[i, ]))
}

phase17o_models <- function(nt) {
  grid <- as.matrix(expand.grid(rep(list(0:1), nt)))
  lapply(seq_len(nrow(grid)), function(i) as.integer(grid[i, ]))
}

phase17o_case <- function(nt = 2L, residual_covariance = "diagonal",
                          updates = FALSE, multiple_sets = FALSE,
                          nonzero = FALSE) {
  fixture <- phase17n_fixture()
  n <- nrow(fixture$X)
  m <- ncol(fixture$X)
  Y <- vapply(seq_len(nt), function(trait) {
    value <- sin(seq_len(n) * (.41 + trait / 13)) +
      cos(seq_len(n) * (.23 + trait / 17))
    value - mean(value)
  }, numeric(n))
  colnames(Y) <- paste0("T", seq_len(nt))
  models <- phase17o_models(nt)
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

phase17o_args <- function(case) {
  f <- case$fixture
  list(
    bed_files = f$paths, n_bed = f$n_bed, cls = f$cls, rows = f$rows,
    af = f$af, Y = case$Y, beta_init = case$beta, b_init = case$b,
    state_init = case$state, sets = case$sets, B = case$B, E = case$E,
    ssb_prior = phase17o_prior_list(case$ssb),
    sse_prior = phase17o_prior_list(case$sse), models = case$models,
    pi = case$pi, nub = case$nub, nue = case$nue,
    updateB = case$updateB, updateE = case$updateE,
    updatePi = case$updatePi,
    residual_covariance = case$residual_covariance,
    nit = case$nit, nburn = case$nburn, nthin = case$nthin,
    seed = case$seed, method = case$method
  )
}

phase17o_call <- function(case) {
  do.call(sblr:::mtblr_bed_internal, phase17o_args(case))
}

phase17o_dense_args <- function(case) {
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
    ssb_prior = phase17o_prior_list(case$ssb),
    sse_prior = phase17o_prior_list(case$sse), models = case$models,
    pi = case$pi, nub = case$nub, nue = case$nue,
    updateB = case$updateB, updateE = case$updateE,
    updatePi = case$updatePi, n = rep.int(nrow(X), nt),
    nit = case$nit, nburn = case$nburn, nthin = case$nthin,
    seed = case$seed, method = case$method
  )
}

phase17o_legacy_matrix <- function(legacy, field) {
  do.call(cbind, legacy[[field]])
}

phase17o_compare_dense <- function(case, tolerance = 1e-10) {
  bed <- phase17o_call(case)
  dense <- do.call(sblr:::mtblr, phase17o_dense_args(case))
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
    phase17o_legacy_matrix(dense, field)
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

phase17o_cleanup <- function(case) {
  unlink(case$fixture$paths)
  invisible()
}
