.is_stblr_raw <- getFromNamespace(".is_stblr_raw", "sblr")
check_stblr_consistency <- getFromNamespace("check_stblr_consistency", "sblr")

make_fake_stblr_raw <- function(
    nt = 1L,
    m = 4L,
    niter = 6L,
    model = "bayesc",
    backend = "csr_bayesc",
    component_names = NULL
) {
  trait_names <- paste0("T", seq_len(nt))
  marker_mat <- function() matrix(stats::rnorm(m * nt), nrow = m, ncol = nt)
  trace_mat <- function() matrix(stats::runif(niter * nt), nrow = niter, ncol = nt)
  variance_mat <- function() diag(nt) * 0.5

  raw <- list(
    schema = list(class = "stblr_raw", version = 1L),
    meta = list(
      nt = nt, m = m, n_trace = niter, model = model, backend = backend
    ),
    marker = list(
      bm = marker_mat(), dm = matrix(stats::runif(m * nt), nrow = m, ncol = nt),
      wy = marker_mat(), r = marker_mat(), b = marker_mat(),
      state = matrix(0, nrow = m, ncol = nt)
    ),
    trace = list(
      vbs = trace_mat(), vgs = trace_mat(), ves = trace_mat(),
      vle = trace_mat(), vld = trace_mat(), pis = trace_mat()
    ),
    variance = list(
      covb = variance_mat(), covg = variance_mat(), cove = variance_mat(),
      vb = variance_mat(), vg = variance_mat(), ve = variance_mat()
    ),
    pi = list(
      final = matrix(0.02, nrow = nt, ncol = 1),
      mean = matrix(0.02, nrow = nt, ncol = 1),
      names = NULL
    ),
    diagnostics = list(),
    chains = NULL,
    prior = NULL,
    group = NULL,
    annotation = NULL,
    component = NULL,
    selection = list()
  )

  if (model %in% c("bayesr", "sbayesrc")) {
    if (is.null(component_names)) {
      component_names <- c("component_0", "component_1", "component_2")
    }
    ncomp <- length(component_names)
    comp_prob <- lapply(seq_len(nt), function(tt) {
      x <- matrix(stats::runif(m * ncomp), nrow = m, ncol = ncomp)
      x / rowSums(x)
    })
    raw$pi$final <- matrix(1 / ncomp, nrow = nt, ncol = ncomp)
    raw$pi$mean <- matrix(1 / ncomp, nrow = nt, ncol = ncomp)
    raw$pi$names <- component_names
    raw$component <- list(
      names = component_names,
      prob = comp_prob
    )
  }

  raw
}

test_that(".is_stblr_raw() recognizes a well-formed schema and rejects others", {
  raw <- make_fake_stblr_raw()
  expect_true(.is_stblr_raw(raw))

  expect_false(.is_stblr_raw(list(a = 1, b = 2)))
  expect_false(.is_stblr_raw(unname(list(1, 2, 3))))

  bad_version <- raw
  bad_version$schema$version <- 2L
  expect_false(.is_stblr_raw(bad_version))

  bad_class <- raw
  bad_class$schema$class <- "something_else"
  expect_false(.is_stblr_raw(bad_class))
})

test_that(".validate_stblr_raw() passes through valid objects and stops on invalid ones", {
  raw <- make_fake_stblr_raw()
  expect_true(.validate_stblr_raw(raw))

  bad_class <- raw
  bad_class$schema$class <- "something_else"
  expect_error(.validate_stblr_raw(bad_class), "Unsupported backend output")
  expect_error(.validate_stblr_raw(bad_class, backend = "csr_bayesc"), "csr_bayesc")
})

test_that(".as_stblr_fit() maps a BayesC-shaped raw object to the canonical fit fields", {
  raw <- make_fake_stblr_raw(nt = 2L, m = 5L, niter = 8L)
  trait_names <- c("trait1", "trait2")
  variable_names <- paste0("m", seq_len(5L))

  fit <- .as_stblr_fit(raw, trait_names = trait_names, variable_names = variable_names)

  for (nm in c("bm", "dm", "wy", "r", "b", "d")) {
    expect_true(is.matrix(fit[[nm]]), info = nm)
    expect_identical(dim(fit[[nm]]), c(5L, 2L), info = nm)
    expect_identical(rownames(fit[[nm]]), variable_names, info = nm)
    expect_identical(colnames(fit[[nm]]), trait_names, info = nm)
  }

  for (nm in c("vbs", "vgs", "ves", "vle", "vld", "pis")) {
    expect_true(is.matrix(fit[[nm]]), info = nm)
    expect_identical(nrow(fit[[nm]]), 8L, info = nm)
    expect_identical(colnames(fit[[nm]]), trait_names, info = nm)
  }

  for (nm in c("covb", "covg", "cove", "vb", "vg", "ve")) {
    expect_true(is.matrix(fit[[nm]]), info = nm)
    expect_identical(dim(fit[[nm]]), c(2L, 2L), info = nm)
  }

  expect_length(fit$pi, 2L)
  expect_length(fit$pim, 2L)

  fit$input <- list(nchains = 1L)
  chk <- check_stblr_consistency(fit, verbose = FALSE)
  expect_true(chk$ok)
})

test_that(".as_stblr_fit() maps a BayesR-shaped raw object with components and pis", {
  raw <- make_fake_stblr_raw(nt = 1L, m = 6L, model = "bayesr", backend = "csr_bayesr")
  variable_names <- paste0("m", seq_len(6L))

  fit <- sblr:::.blr_finalize_fit(
    .as_stblr_fit(raw, trait_names = "trait1", variable_names = variable_names),
    "stblr", "bayesr", "csr"
  )

  expect_true(!is.null(fit$pi_trace))
  expect_true(!is.null(fit$component_probabilities))
  expect_named(fit$component_probabilities, "trait1")
  expect_identical(dim(fit$component_probabilities$trait1), c(6L, 3L))
  expect_equal(fit$dm[, "trait1"],
               1 - fit$component_probabilities$trait1[, "component_0"],
               tolerance = 1e-8)
  expect_true(!is.null(fit$dm_component_mean))

  fit$input <- list(nchains = 1L)
  chk <- check_stblr_consistency(fit, verbose = FALSE)
  expect_true(chk$ok)
})

test_that("check_stblr_consistency() flags malformed comp_prob rows", {
  raw <- make_fake_stblr_raw(nt = 1L, m = 4L, model = "bayesr", backend = "csr_bayesr")
  fit <- sblr:::.blr_finalize_fit(
    .as_stblr_fit(raw, trait_names = "trait1", variable_names = paste0("m", 1:4)),
    "stblr", "bayesr", "csr"
  )

  fit$component_probabilities$trait1[1L, ] <- c(0.5, 0.5, 0.5)
  fit$input <- list(nchains = 1L)

  chk <- check_stblr_consistency(fit, verbose = FALSE)
  expect_false(chk$ok)
  expect_false(chk$checks$ok[match(
    "component_probabilities.trait1.rowsums", chk$checks$check
  )])
})

test_that("check_stblr_consistency() validates fit$selection_s against trait count", {
  raw <- make_fake_stblr_raw(nt = 2L, m = 3L)
  fit <- .as_stblr_fit(
    raw, trait_names = c("trait1", "trait2"), variable_names = paste0("m", 1:3)
  )
  fit$selection_s <- stats::setNames(c(0.1, 0.2), c("trait1", "trait2"))
  fit$input <- list(nchains = 1L)

  chk <- check_stblr_consistency(fit, verbose = FALSE)
  expect_true(chk$ok)
  expect_true(chk$checks$ok[match("selection_s.length", chk$checks$check)])

  fit$selection_s <- c(0.1)
  chk_bad <- check_stblr_consistency(fit, verbose = FALSE)
  expect_false(chk_bad$ok)
  expect_false(chk_bad$checks$ok[match("selection_s.length", chk_bad$checks$check)])
})

test_that("check_stblr_consistency() checks backend-specific fields when fit$input$backend is known", {
  raw <- make_fake_stblr_raw(nt = 1L, m = 3L)
  fit <- .as_stblr_fit(raw, trait_names = "trait1", variable_names = paste0("m", 1:3))
  fit$input <- list(backend = "csr_group_bayesc")

  chk_missing <- check_stblr_consistency(fit, verbose = FALSE)
  expect_false(chk_missing$ok)
  expect_true(any(grepl("^backend\\.csr_group_bayesc\\.", chk_missing$checks$check)))

  fit$group <- list(group_names = "all")
  fit$group_pi <- matrix(0.1, 1, 1)
  fit$group_vb_multiplier <- matrix(1, 1, 1)
  chk_ok <- check_stblr_consistency(fit, verbose = FALSE)
  expect_true(all(chk_ok$checks$ok[grepl("^backend\\.csr_group_bayesc\\.", chk_ok$checks$check)]))
})
