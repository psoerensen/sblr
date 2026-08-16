phase6a_alleles <- function(marker_ids) {
  data.frame(
    marker_id = marker_ids,
    effect = rep(c("A", "C"), length.out = length(marker_ids)),
    other = rep(c("G", "T"), length.out = length(marker_ids)),
    coding = "effect_allele_count", stringsAsFactors = FALSE)
}

phase6a_provider <- function(id, trait, resource, global_markers, score,
                             residual_scale = 1, sample_size = 100,
                             overlap_group = NULL) {
  sblr:::.blr_new_likelihood_provider(
    provider_id = id, trait_ids = trait,
    operator_resource_id = resource$resource_id,
    local_to_global = stats::setNames(
      match(resource$marker_ids, global_markers), resource$marker_ids),
    sufficient_statistics = list(
      score = matrix(score, ncol = 1L,
                     dimnames = list(
                       marker = resource$marker_ids, trait = trait)),
      residual_scale = residual_scale),
    sample_size = stats::setNames(sample_size, trait),
    likelihood_regime = "independent_summary",
    residual_contract = "fixed_provider_residual_scale",
    population = paste0("population_", id), effect_scale = "standardized",
    overlap_group = overlap_group,
    provenance = list(source = "Phase 6A fixture"))
}

phase6a_collection <- function(T = 2L, representation = "csr",
                               provider_order = seq_len(T)) {
  markers <- paste0("m", 1:4)
  traits <- paste0("trait", seq_len(T))
  alleles <- phase6a_alleles(markers)
  map <- sblr:::.blr_new_global_marker_map(markers, alleles)
  cross_product <- matrix(c(
    4, .25, 0, 0,
    .25, 3, 0, 0,
    0, 0, 2.5, .125,
    0, 0, .125, 2
  ), 4L, 4L, byrow = TRUE, dimnames = list(markers, markers))
  resource <- switch(representation,
    csr = sblr:::.blr_csr_from_dense_resource(
      "resource", cross_product, markers, alleles),
    full = sblr:::.blr_block_eigen_resource(
      "resource", cross_product, markers, alleles,
      blocks = list(1:2, 3:4)),
    retained = sblr:::.blr_block_eigen_resource(
      "resource", cross_product, markers, alleles,
      blocks = list(1:2, 3:4), retained_ranks = c(1L, 1L)))
  scores <- lapply(seq_len(T), function(trait) {
    c(.5, -.25, .35, .15) + .0625 * (trait - 1L)
  })
  providers <- lapply(seq_len(T), function(trait) phase6a_provider(
    paste0("provider", trait), traits[trait], resource, markers,
    scores[[trait]], residual_scale = .75 + .25 * trait,
    sample_size = 80 + 10 * trait))
  names(providers) <- vapply(providers, `[[`, character(1), "provider_id")
  providers <- providers[provider_order]
  collection <- sblr:::.blr_new_provider_collection(
    map, list(resource = resource), providers,
    likelihood_regime = "independent_summary",
    analysis_mode = "joint_multitrait")
  list(collection = collection, traits = traits, markers = markers,
       cross_product = cross_product, resource = resource)
}

phase6a_fit <- function(fixture, chains = 1L, cores = 1L,
                        seed = 1901, update = TRUE,
                        burn = 3L, sampling = 6L, thin = 2L) {
  T <- length(fixture$traits)
  Vb <- diag(seq(.35, .35 + .05 * (T - 1L), length.out = T))
  dimnames(Vb) <- list(fixture$traits, fixture$traits)
  sblr:::.blr_cheng_mt_bayesc_summary_qualification(
    collection = fixture$collection, trait_ids = fixture$traits,
    initial_marker_covariance = Vb,
    marker_covariance_prior_df = T + .5,
    marker_covariance_prior_scale = Vb,
    update_marker_covariance = update,
    update_activity_pattern_probability = update,
    burn_in_iterations = burn, sampling_iterations = sampling,
    thin_interval = thin, chains = chains, cores = cores, seed = seed,
    keep_traces = TRUE, memory_limit_bytes = Inf)
}

phase6a_scientific <- function(raw) {
  list(
    posterior = raw$posterior,
    draws = raw$draws,
    final = raw$final,
    derived = raw$derived,
    task_seeds = raw$input$mcmc$task_seeds,
    retained = raw$input$mcmc$retained_transition_indices,
    occupancy = raw$diagnostics$qualification$pattern_occupancy_counts,
    changes = raw$diagnostics$qualification$pattern_change_counts)
}

test_that("Phase 6A marker conditional matches an independent dense derivation", {
  reference <- new.env(parent = baseenv())
  sys.source(test_path("..", "research", "mtblr_covariance",
                       "mtblr_summary_provider_reference.R"),
             envir = reference)
  traits <- paste0("trait", 1:3)
  patterns <- sblr:::.blr_phase4a_patterns(traits)
  score <- c(.4, -.3, .2)
  diagonal <- c(3.5, 2.25, 1.75)
  Vb <- matrix(c(.7, .18, -.08, .18, .9, .12, -.08, .12, .6), 3L)
  probability <- seq_len(8L)
  probability <- probability / sum(probability)
  expected <- reference$mt_summary_pattern_reference(
    score, diagonal, Vb, probability, patterns)
  actual <- sblr:::mtblr_phase6a_summary_pattern_contract_internal(
    score, diagonal, Vb, probability, patterns)
  expect_equal(actual$probability, expected$probability, tolerance = 1e-13)
  expect_equal(actual$log_weight, expected$log_weight, tolerance = 1e-13)
  for (state in 2:8) {
    expect_equal(drop(actual$active_mean[[state]]),
                 expected$active_mean[[state]], tolerance = 1e-13)
    expect_equal(actual$active_covariance[[state]],
                 expected$active_covariance[[state]], tolerance = 1e-13)
  }
})

test_that("complete CSR and full-rank block eigen share one Cheng target", {
  csr <- phase6a_collection(2L, "csr")
  eigen <- phase6a_collection(2L, "full")
  expect_equal(
    unname(sblr:::.blr_operator_matrix(csr$resource)),
    unname(csr$cross_product),
    tolerance = 1e-14)
  expect_equal(
    unname(sblr:::.blr_operator_matrix(eigen$resource)),
    unname(eigen$cross_product),
    tolerance = 1e-12)
  csr_fit <- phase6a_fit(csr, update = FALSE)
  eigen_fit <- phase6a_fit(eigen, update = FALSE)
  expect_identical(csr_fit$draws$joint_states,
                   eigen_fit$draws$joint_states)
  expect_equal(csr_fit$draws$realised_effects,
               eigen_fit$draws$realised_effects, tolerance = 1e-12)
  expect_equal(csr_fit$posterior$activity_pattern_probabilities,
               eigen_fit$posterior$activity_pattern_probabilities,
               tolerance = 0)
})

test_that("retained block eigen targets its reconstructed operator", {
  retained <- phase6a_collection(2L, "retained")
  reconstructed <- sblr:::.blr_operator_matrix(retained$resource)
  markers <- retained$markers
  csr_resource <- sblr:::.blr_csr_from_dense_resource(
    "resource", reconstructed, markers, phase6a_alleles(markers))
  providers <- lapply(seq_along(retained$traits), function(index) {
    original <- retained$collection$providers[[index]]
    phase6a_provider(
      original$provider_id, retained$traits[index], csr_resource, markers,
      original$sufficient_statistics$score[, 1L],
      original$sufficient_statistics$residual_scale,
      original$sample_size[[1L]])
  })
  names(providers) <- vapply(providers, `[[`, character(1), "provider_id")
  csr <- retained
  csr$collection <- sblr:::.blr_new_provider_collection(
    retained$collection$global_marker_map,
    list(resource = csr_resource), providers,
    "independent_summary", "joint_multitrait")
  retained_fit <- phase6a_fit(retained, update = FALSE)
  csr_fit <- phase6a_fit(csr, update = FALSE)
  expect_identical(retained_fit$draws$joint_states,
                   csr_fit$draws$joint_states)
  expect_equal(retained_fit$draws$realised_effects,
               csr_fit$draws$realised_effects, tolerance = 1e-7)
  expect_identical(
    retained_fit$input$data$operator_resources$resource$approximation,
    "retained_rank_approximation")
})

test_that("heterogeneous maps are provider-order and local-order invariant", {
  markers <- paste0("m", 1:4)
  traits <- c("trait1", "trait2")
  alleles <- phase6a_alleles(markers)
  map <- sblr:::.blr_new_global_marker_map(markers, alleles)
  ids1 <- c("m1", "m3", "m4")
  ids2 <- c("m4", "m2")
  C1 <- matrix(c(3, .25, 0, .25, 2, .125, 0, .125, 1.5), 3L,
               dimnames = list(ids1, ids1))
  C2 <- matrix(c(2.5, .25, .25, 2), 2L,
               dimnames = list(ids2, ids2))
  r1 <- sblr:::.blr_csr_from_dense_resource(
    "r1", C1, ids1, alleles[match(ids1, markers), ])
  r2 <- sblr:::.blr_block_eigen_resource(
    "r2", C2, ids2, alleles[match(ids2, markers), ], blocks = list(1:2))
  p1 <- phase6a_provider("p1", traits[1], r1, markers, c(.4, -.2, .1), .8)
  p2 <- phase6a_provider("p2", traits[2], r2, markers, c(.3, -.1), 1.2)
  make <- function(providers) sblr:::.blr_new_provider_collection(
    map, list(r1 = r1, r2 = r2), providers,
    "independent_summary", "joint_multitrait")
  forward <- list(collection = make(list(p1 = p1, p2 = p2)),
                  traits = traits, markers = markers)
  reverse <- list(collection = make(list(p2 = p2, p1 = p1)),
                  traits = traits, markers = markers)
  expect_equal(phase6a_scientific(phase6a_fit(forward, update = FALSE)),
               phase6a_scientific(phase6a_fit(reverse, update = FALSE)),
               tolerance = 0)

  permutation <- c(3L, 1L, 2L)
  ids1p <- ids1[permutation]
  r1p <- sblr:::.blr_csr_from_dense_resource(
    "r1", C1[permutation, permutation, drop = FALSE], ids1p,
    alleles[match(ids1p, markers), ])
  p1p <- phase6a_provider(
    "p1", traits[1], r1p, markers, c(.4, -.2, .1)[permutation], .8)
  permuted <- list(collection = sblr:::.blr_new_provider_collection(
    map, list(r1 = r1p, r2 = r2), list(p1 = p1p, p2 = p2),
    "independent_summary", "joint_multitrait"),
    traits = traits, markers = markers)
  permuted_fit <- phase6a_fit(permuted, update = FALSE)
  expect_equal(
    phase6a_fit(forward, update = FALSE)$posterior,
    permuted_fit$posterior, tolerance = 1e-12)
})

test_that("provider splitting and BED diagonal residual conditionals reduce", {
  traits <- c("trait1", "trait2")
  patterns <- sblr:::.blr_phase4a_patterns(traits)
  Vb <- matrix(c(.6, .15, .15, .8), 2L)
  score_parts <- c(.3, -.1)
  diagonal_parts <- c(2.5, 1.5)
  phi <- c(.5, 2)
  h <- c(sum(score_parts / phi), .4 / 1.25)
  d <- c(sum(diagonal_parts / phi), 3 / 1.25)
  split <- sblr:::mtblr_phase6a_summary_pattern_contract_internal(
    h, d, Vb, rep(.25, 4L), patterns)
  combined <- sblr:::mtblr_phase6a_summary_pattern_contract_internal(
    c(h[1], h[2]), c(d[1], d[2]), Vb, rep(.25, 4L), patterns)
  expect_identical(split, combined)

  Ve <- diag(c(phi[1], 1.25))
  bed <- sblr:::mtblr_phase4a_pattern_contract_internal(
    c(score_parts[1], .4), marker_sum_squares = diagonal_parts[1],
    marker_covariance = Vb, residual_covariance = Ve,
    activity_pattern_probability = rep(.25, 4L),
    activity_patterns = patterns)
  summary <- sblr:::mtblr_phase6a_summary_pattern_contract_internal(
    c(score_parts[1] / phi[1], .4 / 1.25),
    c(diagonal_parts[1] / phi[1], diagonal_parts[1] / 1.25),
    Vb, rep(.25, 4L), patterns)
  expect_equal(summary$probability, bed$probability, tolerance = 1e-13)
  expect_equal(summary$active_mean, bed$active_mean, tolerance = 1e-13)

  fixture <- phase6a_collection(2L, "csr")
  markers <- fixture$markers
  alleles <- phase6a_alleles(markers)
  C1 <- fixture$cross_product
  C2 <- diag(c(.5, .75, 1, 1.25))
  dimnames(C2) <- list(markers, markers)
  r1 <- sblr:::.blr_csr_from_dense_resource("r1", C1, markers, alleles)
  r2 <- sblr:::.blr_csr_from_dense_resource("r2", C2, markers, alleles)
  rt <- sblr:::.blr_csr_from_dense_resource(
    "rt", 2 * C1 + C2, markers, alleles)
  score1 <- c(.25, -.125, .375, .0625)
  score2 <- c(.125, .25, -.0625, .1875)
  score_t2 <- c(.3, -.2, .1, .25)
  split_providers <- list(
    p1a = phase6a_provider("p1a", "trait1", r1, markers, score1, .5),
    p1b = phase6a_provider("p1b", "trait1", r2, markers, score2, 1),
    p2 = phase6a_provider("p2", "trait2", r1, markers, score_t2, 1))
  combined_providers <- list(
    p1 = phase6a_provider(
      "p1", "trait1", rt, markers, 2 * score1 + score2, 1),
    p2 = phase6a_provider("p2", "trait2", r1, markers, score_t2, 1))
  split_fixture <- list(
    collection = sblr:::.blr_new_provider_collection(
      fixture$collection$global_marker_map,
      list(r1 = r1, r2 = r2), split_providers,
      "independent_summary", "joint_multitrait"),
    traits = traits, markers = markers)
  combined_fixture <- list(
    collection = sblr:::.blr_new_provider_collection(
      fixture$collection$global_marker_map,
      list(r1 = r1, rt = rt), combined_providers,
      "independent_summary", "joint_multitrait"),
    traits = traits, markers = markers)
  split_fit <- phase6a_fit(split_fixture, update = FALSE)
  combined_fit <- phase6a_fit(combined_fixture, update = FALSE)
  expect_identical(split_fit$draws$joint_states,
                   combined_fit$draws$joint_states)
  expect_equal(split_fit$draws$realised_effects,
               combined_fit$draws$realised_effects, tolerance = 1e-12)
})

test_that("BED and independent-summary posterior summaries reduce diagonally", {
  dosage <- matrix(c(0, 1, 2, 0, 1, 2, 1, 0), nrow = 1L)
  bed <- tempfile(fileext = ".bed")
  on.exit(unlink(bed), add = TRUE)
  dosage_to_code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  codes <- unname(dosage_to_code[as.character(dosage[1L, ])])
  packed <- vapply(seq(1L, length(codes), by = 4L), function(index) {
    sum(codes[index:(index + 3L)] * c(1L, 4L, 16L, 64L))
  }, integer(1))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), bed)

  sample_ids <- paste0("id", seq_len(ncol(dosage)))
  trait_ids <- c("trait1", "trait2")
  phenotype <- cbind(
    trait1 = c(-1.2, -.4, .8, .2, 1.1, -.7, .5, -.3),
    trait2 = c(.3, -.8, .4, .9, -.2, .7, -.6, -.1))
  rownames(phenotype) <- sample_ids
  allele_frequency <- rowMeans(dosage) / 2
  genotype <- as.numeric(
    (dosage[1L, ] - 2 * allele_frequency) /
      sqrt(2 * allele_frequency * (1 - allele_frequency)))
  cross_product <- sum(genotype^2)
  score <- as.numeric(crossprod(genotype, phenotype))
  residual_scale <- c(.8, 1.25)
  Vb <- matrix(c(.6, .12, .12, .75), 2L,
               dimnames = list(trait_ids, trait_ids))
  probability <- c(.4, .15, .2, .25)
  Glist <- list(
    n = length(sample_ids), ids = sample_ids, bedfiles = bed,
    rsids = list("m1"), rsidsLD = list("m1"), chr = list(1L),
    pos = list(100), af = list(allele_frequency),
    maf = list(pmin(allele_frequency, 1 - allele_frequency)))
  bed_fit <- sblr:::.blr_phase4a_cheng_mt_bed(
    y = phenotype, Glist = Glist,
    fixed_residual_covariance = diag(residual_scale),
    initial_marker_covariance = Vb,
    marker_covariance_prior_df = 2.5,
    marker_covariance_prior_scale = Vb,
    initial_activity_pattern_probability = probability,
    update_marker_covariance = FALSE,
    update_activity_pattern_probability = FALSE,
    burn_in_iterations = 0L, sampling_iterations = 40L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 520)

  marker_ids <- "m1"
  alleles <- phase6a_alleles(marker_ids)
  resource <- sblr:::.blr_csr_from_dense_resource(
    "resource", matrix(cross_product, 1L, 1L,
                       dimnames = list(marker_ids, marker_ids)),
    marker_ids, alleles)
  providers <- lapply(seq_along(trait_ids), function(trait) {
    phase6a_provider(
      paste0("provider", trait), trait_ids[trait], resource, marker_ids,
      score[trait], residual_scale[trait], length(sample_ids))
  })
  names(providers) <- vapply(providers, `[[`, character(1), "provider_id")
  summary_fixture <- list(
    collection = sblr:::.blr_new_provider_collection(
      sblr:::.blr_new_global_marker_map(marker_ids, alleles),
      list(resource = resource), providers,
      "independent_summary", "joint_multitrait"),
    traits = trait_ids, markers = marker_ids)
  summary_fit <- sblr:::.blr_cheng_mt_bayesc_summary_qualification(
    collection = summary_fixture$collection, trait_ids = trait_ids,
    initial_marker_covariance = Vb,
    marker_covariance_prior_df = 2.5,
    marker_covariance_prior_scale = Vb,
    initial_activity_pattern_probability = probability,
    update_marker_covariance = FALSE,
    update_activity_pattern_probability = FALSE,
    burn_in_iterations = 0L, sampling_iterations = 40L,
    thin_interval = 1L, chains = 1L, cores = 1L, seed = 520,
    memory_limit_bytes = Inf)

  expect_identical(summary_fit$draws$joint_states,
                   bed_fit$draws$joint_states)
  expect_equal(summary_fit$draws$realised_effects,
               bed_fit$draws$realised_effects, tolerance = 1e-12)
  expect_equal(summary_fit$posterior$activity_pattern_probabilities,
               bed_fit$posterior$activity_pattern_probabilities,
               tolerance = 0)
  expect_equal(summary_fit$posterior$realised_effect_mean,
               bed_fit$posterior$realised_effect_mean, tolerance = 1e-12)
})

test_that("Phase 6A T=2 and T=3 short chains return truthful raw v2", {
  for (T in 2:3) {
    raw <- phase6a_fit(phase6a_collection(T, if (T == 2) "csr" else "full"))
    expect_true(sblr:::validate_blr_raw_v2(raw))
    expect_identical(raw$input$data$likelihood_regime,
                     "independent_summary")
    expect_null(raw$derived$predictions)
    expect_null(raw$draws$residual_covariance)
    expect_null(raw$final$residual_covariance)
    expect_equal(dim(raw$draws$realised_effects)[4L], T)
    expect_equal(dim(raw$posterior$activity_pattern_probabilities)[2L], 2^T)
    expect_false(raw$diagnostics$qualification$current_legacy_mt_route_used)
    if (T == 3L) {
      final_effect <- raw$final$realised_effects[1L, , , drop = TRUE]
      provider_residual <- raw$diagnostics$qualification$
        final_provider_residual_scores$chain1
      for (provider_id in names(raw$input$data$providers)) {
        provider <- raw$input$data$providers[[provider_id]]
        resource <- raw$input$data$operator_resources[[
          provider$operator_resource_id]]
        trait <- match(provider$trait_ids, raw$input$data$trait_ids)
        expected <- as.numeric(
          provider$sufficient_statistics$score[, 1L] -
            sblr:::.blr_operator_apply(
              resource,
              final_effect[provider$local_to_global, trait]))
        expect_equal(provider_residual[[provider_id]], expected,
                     tolerance = 1e-7)
      }
    }
  }
})

test_that("Phase 6A chain scheduling is deterministic and task private", {
  fixture <- phase6a_collection(3L, "csr")
  serial <- phase6a_fit(fixture, chains = 2L, cores = 1L, seed = 812)
  parallel <- phase6a_fit(fixture, chains = 2L, cores = 2L, seed = 812)
  expect_identical(phase6a_scientific(serial), phase6a_scientific(parallel))
  expect_identical(serial$input$mcmc$task_seeds,
                   parallel$input$mcmc$task_seeds)
  expect_identical(parallel$diagnostics$workers$logical_task_order,
                   c("chain:0", "chain:1"))
  workers <- parallel$diagnostics$workers
  if (isTRUE(workers$openmp_available) &&
      workers$runtime_maximum_workers >= 2L) {
    expect_identical(workers$actual_team_size, 2L)
    expect_identical(sort(unique(workers$task_worker_ids)), c(0L, 1L))
  } else {
    skip("OpenMP runtime cannot supply two sampler workers")
  }
})

test_that("Phase 6A rejects overlap and preflights summary allocations", {
  fixture <- phase6a_collection(2L, "csr")
  bad <- fixture$collection
  bad$providers[[1L]]$overlap_group <- "overlap1"
  expect_error(
    sblr:::.blr_cheng_mt_bayesc_summary_qualification(
      bad, fixture$traits, diag(2), 2.5, diag(2),
      burn_in_iterations = 0L, sampling_iterations = 2L),
    "overlap-aware")
  estimate1 <- sblr:::.blr_phase6a_memory_estimate(
    fixture$collection, 2L, 1L, 2L, 4L, TRUE,
    memory_limit_bytes = Inf, enforce = FALSE)
  estimate2 <- sblr:::.blr_phase6a_memory_estimate(
    fixture$collection, 2L, 2L, 4L, 8L, TRUE,
    memory_limit_bytes = Inf, enforce = FALSE)
  expect_gt(estimate2$estimated_peak_incremental_bytes,
            estimate1$estimated_peak_incremental_bytes)
  expect_gt(estimate1$components[["summary_operator_resources"]], 0)
  expect_gt(
    estimate1$components[["summary_chain_provider_residual_scores"]], 0)
  expect_error(sblr:::.blr_phase6a_memory_estimate(
    fixture$collection, 2L, 1L, 2L, 4L, TRUE,
    memory_limit_bytes = 1, enforce = TRUE),
    "before native allocation")
})
