test_that("make_sbayesrc_alpha_init returns deterministic valid defaults", {
 A <- cbind(
  intercept = rep(1, 4),
  binary = c(0, 1, 0, 1)
 )
 gamma <- c(0, 0.01, 0.1, 1)

 out <- sblr::make_sbayesrc_alpha_init(
  A,
  gamma = gamma,
  pi_init = 0.2,
  active_comp_weights = c(1, 2, 1)
 )

 expect_equal(dim(out$alpha_init), c(ncol(A), length(gamma) - 1L))
 expect_length(out$sigmaSqAlpha_init, length(gamma) - 1L)
 expect_equal(sum(out$active_comp_weights), 1)
 expect_equal(unname(out$active_comp_weights), c(0.25, 0.5, 0.25))
 expect_equal(sum(out$component_prob_init), 1)
 expect_identical(out$gamma, gamma)
 expect_identical(out$annotation_names, colnames(A))
 expect_identical(
  out$component_names,
  c("null", "gamma_0.01", "gamma_0.10", "gamma_1.00")
 )
 expect_true(all(is.finite(out$alpha_init)))
 expect_equal(unname(out$alpha_init["binary", ]), rep(0, length(gamma) - 1L))
})

test_that("make_sbayesrc_alpha_init accepts supplied initial values", {
 A <- cbind(intercept = rep(1, 3), binary = c(0, 1, 0))
 gamma <- c(0, 0.01, 0.1)
 alpha_init <- matrix(
  seq_len(ncol(A) * (length(gamma) - 1L)) / 10,
  nrow = ncol(A)
 )
 sigmaSqAlpha_init <- c(0.5, 1.5)

 out <- sblr::make_sbayesrc_alpha_init(
  A,
  gamma = gamma,
  alpha_init = alpha_init,
  sigmaSqAlpha_init = sigmaSqAlpha_init
 )

 expect_equal(unname(out$alpha_init), alpha_init)
 expect_equal(unname(out$sigmaSqAlpha_init), sigmaSqAlpha_init)
 expect_identical(rownames(out$alpha_init), colnames(A))
 expect_identical(colnames(out$alpha_init), c("step_1", "step_2"))
})

test_that("make_sbayesrc_alpha_init rejects invalid inputs", {
 A <- cbind(intercept = rep(1, 3), binary = c(0, 1, 0))

 expect_error(sblr::make_sbayesrc_alpha_init(A, gamma = c(0.1, 1)), "gamma\\[1\\]")
 expect_error(sblr::make_sbayesrc_alpha_init(A, gamma = c(0, 0)), "positive")
 expect_error(sblr::make_sbayesrc_alpha_init(A, pi_init = 0), "pi_init")
 expect_error(sblr::make_sbayesrc_alpha_init(A, pi_init = 1), "pi_init")
 expect_error(
  sblr::make_sbayesrc_alpha_init(A, active_comp_weights = c(1, -1, 1)),
  "active_comp_weights"
 )
 expect_error(
  sblr::make_sbayesrc_alpha_init(A, alpha_init = matrix(0, nrow = 1, ncol = 3)),
  "alpha_init"
 )
 expect_error(
  sblr::make_sbayesrc_alpha_init(A, sigmaSqAlpha_init = c(1, 0, 1)),
  "sigmaSqAlpha_init"
 )
 expect_error(
  sblr::make_sbayesrc_alpha_init(A, sigmaSqAlpha_init = c(1, 1)),
  "sigmaSqAlpha_init"
 )
})

test_that("format_sbayesrc_csr_fit formats a synthetic sampler result", {
 nt <- 2L
 m <- 4L
 n_anno <- 3L
 gamma <- c(0, 0.01, 0.1)
 nstep <- length(gamma) - 1L
 ncomp <- length(gamma)
 trait_names <- c("trait_a", "trait_b")
 variable_names <- paste0("marker_", seq_len(m))
 annotation_names <- paste0("annot_", seq_len(n_anno))

 marker_slot <- function(offset = 0) {
  lapply(seq_len(nt), function(trait) seq_len(m) + offset + trait)
 }
 trace_slot <- function(offset = 0) {
  lapply(seq_len(nt), function(trait) seq_len(3) + offset + trait)
 }
 covariance_slot <- function(scale = 1) {
  list(c(scale, 0.1), c(0.1, scale))
 }

 fit <- vector("list", 24)
 fit[1:7] <- lapply(seq_len(7), marker_slot)
 fit[8:10] <- lapply(seq_len(3), trace_slot)
 fit[11:16] <- lapply(seq_len(6), covariance_slot)
 fit[17:18] <- lapply(seq_len(2), function(i) {
  list(c(0.9, 0.1), c(0.8, 0.2))
 })
 fit[[19]] <- lapply(seq_len(nt), function(i) {
  seq_len(n_anno * nstep) / 10 + i
 })
 fit[[20]] <- lapply(seq_len(nt), function(i) seq_len(nstep) + i)
 fit[21:22] <- lapply(seq_len(2), trace_slot)
 fit[[23]] <- lapply(seq_len(nt), function(i) {
  rep(c(0.8, 0.15, 0.05), m)
 })
 fit[[24]] <- lapply(seq_len(nt), function(i) c(3, 1, 0))

 out <- sblr:::format_sbayesrc_csr_fit(
  fit,
  nt = nt,
  m = m,
  gamma = gamma,
  n_anno = n_anno,
  trait_names = trait_names,
  variable_names = variable_names,
  annotation_names = annotation_names
 )

 expect_named(
  out,
  c(
   "bm", "dm", "wy", "r", "b", "component", "marker_index",
   "vbs", "vgs", "ves", "covb", "covg", "cove", "vb", "vg", "ve",
   "pi", "pim", "alpha_flat", "alpha", "sigmaSqAlpha", "vle", "vld",
   "comp_prob", "ncomp", "rb", "rg", "re"
  )
 )
 expect_equal(dim(out$bm), c(m, nt))
 expect_equal(dim(out$dm), c(m, nt))
 expect_equal(dim(out$vbs), c(3, nt))
 expect_equal(dim(out$vgs), c(3, nt))
 expect_equal(dim(out$ves), c(3, nt))
 expect_length(out$alpha, nt)
 expect_equal(dim(out$alpha[[1]]), c(n_anno, nstep))
 expect_equal(dim(out$sigmaSqAlpha), c(nt, nstep))
 expect_length(out$comp_prob, nt)
 expect_equal(dim(out$comp_prob[[1]]), c(m, ncomp))
 expect_equal(dim(out$ncomp), c(nt, ncomp))
 expect_identical(rownames(out$bm), variable_names)
 expect_identical(colnames(out$bm), trait_names)
 expect_identical(names(out$alpha), trait_names)
 expect_identical(names(out$comp_prob), trait_names)
})
