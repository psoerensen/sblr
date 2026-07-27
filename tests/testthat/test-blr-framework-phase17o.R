test_that("Phase 17O production marker kernel matches the independent oracle", {
  configurations <- list(
    list(B = matrix(c(1, .2, .2, .8), 2),
         E = matrix(c(1, .3, .3, 1.2), 2),
         score = c(.7, -.2), w = 3.4),
    list(B = matrix(c(1, .1, .04, .1, .8, .06, .04, .06, .7), 3),
         E = matrix(c(1, .2, .08, .2, 1.1, .05, .08, .05, .9), 3),
         score = c(.7, -.2, .4), w = 2.7),
    list(B = matrix(c(.9, .18, .18, .7), 2),
         E = diag(c(.8, 1.3)), score = c(-.4, .9), w = 4.1)
  )
  for (configuration in configurations) {
    nt <- length(configuration$score)
    model_matrix <- as.matrix(expand.grid(rep(list(0:1), nt)))
    models <- lapply(seq_len(nrow(model_matrix)),
                     function(i) as.integer(model_matrix[i, ]))
    pi <- seq_len(length(models))
    pi[2L] <- 0
    pi <- pi / sum(pi)
    native <- sblr:::mtblr_bed_marker_contract_internal(
      configuration$score, configuration$w, configuration$B,
      configuration$E, models, pi)
    oracle <- phase17n_marker_conditional(
      configuration$score, configuration$w, configuration$B,
      configuration$E, model_matrix, pi)
    expect_equal(native$probability, oracle$probabilities, tolerance = 1e-12)
    expect_equal(native$log_weight, oracle$log_weights, tolerance = 1e-12)
    expect_equal(native$C, lapply(oracle$models, `[[`, "C"),
                 tolerance = 1e-12)
    expect_equal(native$rhs, lapply(oracle$models, `[[`, "rhs"),
                 tolerance = 1e-12)
    expect_equal(native$mean, lapply(oracle$models, `[[`, "mean"),
                 tolerance = 1e-12)
    expect_equal(native$covariance,
                 lapply(oracle$models, `[[`, "covariance"),
                 tolerance = 1e-12)
  }
})

test_that("Phase 17O diagonal mode reduces exactly to corrected dense MT", {
  configurations <- list(
    fixed = c(FALSE, FALSE, FALSE),
    B = c(TRUE, FALSE, FALSE),
    E = c(FALSE, TRUE, FALSE),
    pi = c(FALSE, FALSE, TRUE),
    all = c(TRUE, TRUE, TRUE)
  )
  for (flags in configurations) {
    case <- phase17o_case()
    on.exit(phase17o_cleanup(case), add = TRUE)
    case$updateB <- flags[1L]
    case$updateE <- flags[2L]
    case$updatePi <- flags[3L]
    phase17o_compare_dense(case)
  }
  for (case in list(
    phase17o_case(updates = TRUE, multiple_sets = TRUE),
    phase17o_case(nt = 3L, updates = TRUE),
    phase17o_case(nt = 1L, updates = TRUE)
  )) {
    on.exit(phase17o_cleanup(case), add = TRUE)
    phase17o_compare_dense(case)
  }
})

test_that("Phase 17O full E execution satisfies sample-space identities", {
  for (case in list(
    phase17o_case(residual_covariance = "full", updates = FALSE),
    phase17o_case(residual_covariance = "full", updates = TRUE),
    phase17o_case(nt = 3L, residual_covariance = "full", updates = TRUE),
    phase17o_case(residual_covariance = "full", updates = TRUE,
                  multiple_sets = TRUE, nonzero = TRUE)
  )) {
    on.exit(phase17o_cleanup(case), add = TRUE)
    raw <- phase17o_call(case)
    X <- case$fixture$X
    residual <- case$Y - X %*% raw$marker$b
    genetic <- X %*% raw$marker$b
    expect_identical(.validate_mtblr_raw(raw), raw)
    expect_identical(raw$meta$backend, "mt_bed_bayesc")
    expect_identical(raw$meta$data_level, "individual")
    expect_equal(raw$marker$wy, unname(crossprod(X, case$Y)),
                 tolerance = 1e-12)
    expect_equal(raw$marker$r, unname(crossprod(X, residual)),
                 tolerance = 1e-12)
    expect_equal(raw$variance$vg, crossprod(genetic) / nrow(X),
                 tolerance = 1e-12)
    expect_true(all(raw$marker$state %in% 0:1))
    expect_equal(sum(raw$pi$final), 1, tolerance = 1e-12)
    expect_true(min(eigen(raw$variance$vb, symmetric = TRUE,
                          only.values = TRUE)$values) > 0)
    expect_true(min(eigen(raw$variance$ve, symmetric = TRUE,
                          only.values = TRUE)$values) > 0)
    expect_true(min(eigen(raw$variance$vg, symmetric = TRUE,
                          only.values = TRUE)$values) > -1e-10)
    expect_identical(raw$diagnostics$mt_bed$owner, "PackedBedMatrix")
    expect_identical(raw$diagnostics$mt_bed$view,
                     "BedPackedGenotypeView")
  }
})

test_that("Phase 17O preserves supplied full E when its update is disabled", {
  case <- phase17o_case(residual_covariance = "full")
  on.exit(phase17o_cleanup(case), add = TRUE)
  raw <- phase17o_call(case)
  expect_equal(raw$variance$ve, case$E, tolerance = 0)
  expect_true(any(raw$variance$ve[row(raw$variance$ve) !=
                                    col(raw$variance$ve)] != 0))
  expect_identical(raw$diagnostics$mt_bed$full_e_updates, 0)
  expect_identical(raw$diagnostics$mt_bed$diagonal_e_updates, 0)
})

test_that("Phase 17O validates initial latent, effective, and state values", {
  case <- phase17o_case(nonzero = TRUE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  expect_silent(phase17o_call(case))
  expect_true(case$beta[[2L]][2L] != 0 && case$b[[2L]][2L] == 0)

  args <- phase17o_args(case)
  bad <- args
  bad$state_init[[1L]][1L] <- 2L
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "binary")
  bad <- args
  bad$b_init[[2L]][2L] <- .1
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "inactive")
  bad <- args
  bad$b_init[[1L]][1L] <- .08
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "must agree")
  bad <- args
  bad$state_init[[1L]][3L] <- 1L
  bad$models <- list(c(0L, 0L), c(1L, 1L))
  bad$pi <- c(.5, .5)
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "supplied model")
  bad <- args
  bad$beta_init[[1L]][1L] <- Inf
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "finite")
})

test_that("Phase 17O binding rejects malformed data and execution contracts", {
  case <- phase17o_case()
  on.exit(phase17o_cleanup(case), add = TRUE)
  args <- phase17o_args(case)
  mutations <- list(
    list(name = "bed_files", value = character(), message = "nonempty"),
    list(name = "n_bed", value = 1L, message = "greater than one"),
    list(name = "rows", value = integer(), message = "nonempty"),
    list(name = "rows", value = c(1L, 1L), message = "duplicates"),
    list(name = "af", value = c(0, args$af[-1L]), message = "in \\(0, 1\\)"),
    list(name = "residual_covariance", value = "mask", message = "invalid"),
    list(name = "method", value = 3L, message = "invalid"),
    list(name = "nit", value = 0L, message = "invalid")
  )
  for (mutation in mutations) {
    bad <- args
    bad[[mutation$name]] <- mutation$value
    expect_error(do.call(sblr:::mtblr_bed_internal, bad), mutation$message)
  }
  bad <- args
  bad$Y[1L, 1L] <- bad$Y[1L, 1L] + 1
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "centered")
  bad <- args
  colnames(bad$Y) <- c("T", "T")
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "unique")
  bad <- args
  bad$Y[1L, 1L] <- NA_real_
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "finite")
  bad <- args
  bad$sets <- list(c(0L, 1L), c(1L, 2L, 3L, 4L))
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "disjoint")
  bad <- args
  bad$models <- bad$models[-1L]
  bad$pi <- rep(1 / length(bad$models), length(bad$models))
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "null")
  bad <- args
  bad$E[1L, 2L] <- bad$E[2L, 1L] <- .1
  expect_error(do.call(sblr:::mtblr_bed_internal, bad), "diagonal")
})

test_that("Phase 17O execution is fit-local and seed reproducible", {
  case <- phase17o_case(residual_covariance = "full", updates = TRUE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  first <- phase17o_call(case)
  second <- phase17o_call(case)
  expect_identical(first, second)

  other <- case
  other$seed <- case$seed + 1L
  third <- phase17o_call(other)
  expect_false(identical(first$marker$b, third$marker$b))

  dense <- phase17o_dense_args(phase17o_case())
  invisible(do.call(sblr:::mtblr, dense))
  after_summary <- phase17o_call(case)
  expect_identical(first, after_summary)
})

test_that("Phase 17O is reproducible in a fresh process", {
  skip_if_not_installed("callr")
  case <- phase17o_case(residual_covariance = "full", updates = TRUE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  args <- phase17o_args(case)
  expected <- phase17o_call(case)
  root <- if (exists("blr_has_source_tree", mode = "function") &&
              blr_has_source_tree()) {
    blr_source_root
  } else {
    NULL
  }
  observed <- callr::r(function(args, root) {
    if (!is.null(root)) {
      pkgload::load_all(root, compile = FALSE, quiet = TRUE)
    } else {
      library(sblr)
    }
    do.call(getFromNamespace("mtblr_bed_internal", "sblr"), args)
  }, list(args = args, root = root))
  expect_identical(observed, expected)
})

test_that("Phase 17O architecture retains exactly two internal routes", {
  # Phase 17P adds only the public R adapter; these Phase 17O native routes
  # remain internal and singular.
  expect_true("mtblr_bed" %in% getNamespaceExports("sblr"))
  expect_true(exists("mtblr_bed", envir = asNamespace("sblr"),
                     inherits = FALSE))
  root <- blr_repo_path()
  skip_if(is.null(root), "source architecture requires repository source")
  mtblr <- readLines(file.path(root, "src", "mtblr.cpp"), warn = FALSE)
  core <- readLines(file.path(root, "src", "blr_mt_bed_core_impl.h"),
                    warn = FALSE)
  exports <- readLines(file.path(root, "src", "RcppExports.cpp"), warn = FALSE)
  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  expect_equal(sum(grepl("mtblr_bed_internal\\(", mtblr, fixed = FALSE)), 1L)
  expect_equal(sum(grepl("mtblr_bed_marker_contract_internal\\(",
                         mtblr, fixed = FALSE)), 1L)
  expect_equal(sum(grepl("run_mt_bed_bayesc_core\\(", core)), 1L)
  expect_equal(sum(grepl("PackedBedMatrix owner=", mtblr, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("mtblr_bed_internal", exports, fixed = TRUE)), 5L)
  expect_equal(sum(grepl("mtblr_bed_marker_contract_internal", exports,
                         fixed = TRUE)), 5L)
  expect_equal(sum(grepl("export\\(mtblr_bed\\)", namespace)), 1L)
  expect_false(any(grepl("mtblr_eigen\\(", core)))
  expect_false(any(grepl("#pragma omp|omp_get|OpenMP", core,
                         ignore.case = TRUE)))
})
