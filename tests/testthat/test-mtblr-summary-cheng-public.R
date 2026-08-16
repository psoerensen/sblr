phase6b_alleles <- function(marker_ids) {
  data.frame(
    marker_id = marker_ids,
    effect = rep(c("A", "C"), length.out = length(marker_ids)),
    other = rep(c("G", "T"), length.out = length(marker_ids)),
    coding = "effect_allele_count", stringsAsFactors = FALSE)
}

phase6b_csr_resource <- function(id, matrix, marker_ids,
                                 approximation = "exact_declared_operator") {
  stopifnot(identical(dim(matrix), rep(length(marker_ids), 2L)))
  rows <- lapply(seq_along(marker_ids), function(row) {
    columns <- which(matrix[row, ] != 0)
    list(columns = columns, values = matrix[row, columns])
  })
  list(
    resource_id = id, marker_ids = marker_ids,
    alleles = phase6b_alleles(marker_ids),
    genotype_coding = "effect_allele_count", centering = "declared",
    standardization = "declared", operator_scale = "cross_product",
    approximation = approximation, provenance = list(source = "test"),
    csr = list(
      row_ptr = c(0, cumsum(vapply(rows, function(x) length(x$columns),
                                    integer(1)))),
      column_index = as.integer(unlist(lapply(rows, `[[`, "columns"))),
      values = as.numeric(unlist(lapply(rows, `[[`, "values"))),
      diagonal = as.numeric(diag(matrix)),
      complete = identical(approximation, "exact_declared_operator")),
    diagonal = as.numeric(diag(matrix)))
}

phase6b_block_resource <- function(id, matrix, marker_ids,
                                   retained_rank = NULL) {
  split <- split(seq_along(marker_ids),
                 rep(seq_len(ceiling(length(marker_ids) / 2)), each = 2,
                     length.out = length(marker_ids)))
  blocks <- lapply(seq_along(split), function(index) {
    markers <- split[[index]]
    decomposition <- eigen(matrix[markers, markers, drop = FALSE],
                           symmetric = TRUE)
    rank <- if (is.null(retained_rank)) length(markers) else
      min(retained_rank, length(markers))
    list(
      block_id = paste0("block", index), marker_ids = marker_ids[markers],
      eigenvectors = decomposition$vectors[, seq_len(rank), drop = FALSE],
      eigenvalues = decomposition$values[seq_len(rank)])
  })
  retained <- !is.null(retained_rank) &&
    any(vapply(blocks, function(x) ncol(x$eigenvectors) < length(x$marker_ids),
               logical(1)))
  list(
    resource_id = id, marker_ids = marker_ids,
    alleles = phase6b_alleles(marker_ids),
    genotype_coding = "effect_allele_count", centering = "declared",
    standardization = "declared", operator_scale = "cross_product",
    approximation = if (retained) "retained_rank_declared_operator" else
      "exact_declared_block_diagonal_operator",
    provenance = list(source = "test"), blocks = blocks)
}

phase6b_provider <- function(id, trait, resource, score,
                             residual_scale = 1, overlap_group = NULL) {
  list(
    provider_id = id, trait_id = trait,
    operator_resource_id = resource$resource_id,
    score = stats::setNames(as.numeric(score), resource$marker_ids),
    sample_size = 100, residual_scale = residual_scale,
    likelihood_regime = "independent_summary", effect_scale = "standardized",
    population = paste0("population_", id), overlap_group = overlap_group,
    provenance = list(source = "test"))
}

phase6b_fixture <- function(T = 2L, representation = "csr",
                            retained_rank = NULL) {
  marker_ids <- paste0("m", 1:4)
  trait_ids <- paste0("trait", seq_len(T))
  cross_product <- matrix(c(
    4, .25, 0, 0, .25, 3, 0, 0,
    0, 0, 2.5, .125, 0, 0, .125, 2
  ), 4L, byrow = TRUE, dimnames = list(marker_ids, marker_ids))
  resource <- if (identical(representation, "csr")) {
    phase6b_csr_resource("resource", cross_product, marker_ids)
  } else {
    phase6b_block_resource("resource", cross_product, marker_ids,
                           retained_rank)
  }
  providers <- lapply(seq_len(T), function(trait) phase6b_provider(
    paste0("provider", trait), trait_ids[[trait]], resource,
    c(.5, -.25, .35, .15) + .0625 * (trait - 1L),
    residual_scale = .75 + .25 * trait))
  covariance <- diag(seq(.35, .35 + .05 * (T - 1L), length.out = T))
  dimnames(covariance) <- list(trait_ids, trait_ids)
  list(
    providers = providers, resources = list(resource), markers = marker_ids,
    alleles = phase6b_alleles(marker_ids), traits = trait_ids,
    covariance = covariance)
}

phase6b_public_args <- function(fixture, chains = 1L, cores = 1L,
                                seed = 2601L) {
  list(
    providers = fixture$providers,
    operator_resources = fixture$resources,
    global_marker_ids = fixture$markers,
    global_alleles = fixture$alleles,
    trait_ids = fixture$traits,
    initial_marker_covariance = fixture$covariance,
    marker_covariance_prior_degrees_of_freedom =
      length(fixture$traits) + .5,
    marker_covariance_prior_scale = fixture$covariance,
    nit = 6L, nburn = 3L, nthin = 2L, seed = seed,
    nchains = chains, ncores = cores, keep_traces = TRUE,
    memory_limit_bytes = Inf)
}

phase6b_scientific <- function(raw) {
  list(
    posterior = raw$posterior, draws = raw$draws, final = raw$final,
    derived = raw$derived, task_seeds = raw$input$mcmc$task_seeds,
    retained = raw$input$mcmc$retained_transition_indices,
    occupancy = raw$diagnostics$qualification$pattern_occupancy_counts,
    changes = raw$diagnostics$qualification$pattern_change_counts)
}

test_that("public CSR is identical to the qualified Phase 6A kernel", {
  fixture <- phase6b_fixture(2L, "csr")
  args <- phase6b_public_args(fixture)
  fit <- do.call(mtblr_csr, args)
  collection <- sblr:::.mtblr_summary_public_collection(
    fixture$providers, fixture$resources, fixture$markers, fixture$alleles,
    fixture$traits, "csr")
  internal <- sblr:::.blr_cheng_mt_bayesc_summary_qualification(
    collection, fixture$traits, fixture$covariance,
    length(fixture$traits) + .5, fixture$covariance,
    burn_in_iterations = 3L, sampling_iterations = 6L,
    thin_interval = 2L, seed = 2601L, keep_traces = TRUE,
    memory_limit_bytes = Inf)
  raw <- attr(fit, "blr_raw")
  expect_identical(phase6b_scientific(raw), phase6b_scientific(internal))
  expect_silent(sblr:::validate_blr_raw_v2(raw))
  expect_s3_class(fit, "mtblr_fit")
  expect_null(fit$predictions)
  expect_null(fit$residual_covariance_draws)
})

test_that("public block eigen is identical and retained-rank metadata is truthful", {
  fixture <- phase6b_fixture(3L, "block", retained_rank = 1L)
  args <- phase6b_public_args(fixture, seed = 2602L)
  fit <- do.call(mtblr_block_eigen, args)
  collection <- sblr:::.mtblr_summary_public_collection(
    fixture$providers, fixture$resources, fixture$markers, fixture$alleles,
    fixture$traits, "block_eigen")
  internal <- sblr:::.blr_cheng_mt_bayesc_summary_qualification(
    collection, fixture$traits, fixture$covariance,
    length(fixture$traits) + .5, fixture$covariance,
    burn_in_iterations = 3L, sampling_iterations = 6L,
    thin_interval = 2L, seed = 2602L, keep_traces = TRUE,
    memory_limit_bytes = Inf)
  expect_identical(phase6b_scientific(attr(fit, "blr_raw")),
                   phase6b_scientific(internal))
  expect_true(all(vapply(
    attr(fit, "blr_raw")$input$data$operator_resources,
    function(x) identical(x$operator_type, "retained_rank_block_eigen"),
    logical(1))))
})

test_that("public providers support heterogeneous marker maps", {
  fixture <- phase6b_fixture(2L, "csr")
  ids <- rev(fixture$markers[-1L])
  matrix <- diag(c(2, 3, 4))
  dimnames(matrix) <- list(ids, ids)
  second <- phase6b_csr_resource("resource2", matrix, ids)
  second$alleles <- fixture$alleles[match(ids, fixture$markers), , drop = FALSE]
  fixture$resources <- list(fixture$resources[[1L]], second)
  fixture$providers[[2L]] <- phase6b_provider(
    "provider2", fixture$traits[[2L]], second, c(.2, -.1, .3), 1.25)
  fit <- do.call(mtblr_csr, phase6b_public_args(fixture, seed = 2603L))
  raw <- attr(fit, "blr_raw")
  expect_identical(
    unname(raw$input$data$providers$provider2$local_to_global),
    match(ids, fixture$markers))
  expect_true(all(is.finite(fit$bm)))
})

test_that("public summary chains are scheduler neutral", {
  fixture <- phase6b_fixture(4L, "block")
  serial <- do.call(mtblr_block_eigen,
                    phase6b_public_args(fixture, 2L, 1L, 2604L))
  parallel <- do.call(mtblr_block_eigen,
                      phase6b_public_args(fixture, 2L, 2L, 2604L))
  expect_identical(
    phase6b_scientific(attr(serial, "blr_raw")),
    phase6b_scientific(attr(parallel, "blr_raw")))
})

test_that("public summary output follows raw-v2 and formatted contracts", {
  fit <- do.call(mtblr_csr, phase6b_public_args(phase6b_fixture()))
  raw <- attr(fit, "blr_raw")
  expect_identical(fit$bm, raw$posterior$realised_effect_mean)
  expect_identical(fit$dm, raw$posterior$pips)
  expect_identical(fit$cov_b_mean, raw$posterior$marker_covariance_mean)
  expect_identical(raw$derived$predictions, NULL)
  expect_identical(raw$draws$residual_covariance, NULL)
  expect_identical(raw$posterior$residual_covariance_mean, NULL)
  expect_identical(raw$diagnostics$qualification$status,
                   "publicly_supported")
  expect_false(any(c("pi", "pis", "pim") %in% names(fit)))
})

test_that("unsupported targets and the historical hybrid are unreachable", {
  fixture <- phase6b_fixture()
  args <- phase6b_public_args(fixture)
  expect_error(do.call(mtblr_csr, c(args, method = "bayesr")),
               "only method = 'bayesc'", fixed = TRUE)
  fixture$providers[[1L]]$overlap_group <- "shared_samples"
  expect_error(do.call(mtblr_csr, phase6b_public_args(fixture)),
               "overlap-aware likelihood", fixed = TRUE)
  expect_false(exists("mtblr_csr_chains_raw_internal",
                      envir = asNamespace("sblr"), inherits = FALSE))
  expect_false(exists("mtblr_block_eigen_chains_raw_internal",
                      envir = asNamespace("sblr"), inherits = FALSE))
  expect_false(exists(".mtblr_csr_legacy_covariance_hybrid",
                      envir = asNamespace("sblr"), inherits = FALSE))
})
