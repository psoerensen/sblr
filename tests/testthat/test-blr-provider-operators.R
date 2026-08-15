phase2_marker_fixture <- function() {
  marker_ids <- paste0("m", 1:4)
  alleles <- data.frame(
    marker_id = marker_ids,
    effect = c("A", "C", "G", "T"),
    other = c("G", "T", "A", "C"),
    coding = "effect_allele_count",
    stringsAsFactors = FALSE)
  cross_product <- matrix(c(
    4.0, 0.8, 0.0, 0.0,
    0.8, 3.0, 0.0, 0.0,
    0.0, 0.0, 2.5, 0.4,
    0.0, 0.0, 0.4, 2.0
  ), 4L, 4L, byrow = TRUE,
  dimnames = list(marker_ids, marker_ids))
  list(
    marker_ids = marker_ids,
    alleles = alleles,
    map = sblr:::.blr_new_global_marker_map(marker_ids, alleles),
    cross_product = cross_product)
}

phase2_resource_alleles <- function(fixture, marker_ids) {
  rows <- match(marker_ids, fixture$marker_ids)
  fixture$alleles[rows, c("marker_id", "effect", "other"), drop = FALSE]
}

phase2_summary_provider <- function(id, traits, resource, global_markers,
                                    score = NULL) {
  map <- match(resource$marker_ids, global_markers)
  if (is.null(score)) {
    score <- matrix(
      seq_len(length(map) * length(traits)) / 10,
      nrow = length(map), ncol = length(traits),
      dimnames = list(marker = resource$marker_ids, trait = traits))
  }
  sblr:::.blr_new_likelihood_provider(
    provider_id = id, trait_ids = traits,
    operator_resource_id = resource$resource_id,
    local_to_global = stats::setNames(map, resource$marker_ids),
    sufficient_statistics = list(score = score),
    sample_size = stats::setNames(rep(120, length(traits)), traits),
    likelihood_regime = "independent_summary",
    residual_contract = "marginal_summary_error",
    population = "fixture_population", effect_scale = "standardized",
    overlap_group = NULL,
    provenance = list(source = "Phase 2 deterministic fixture"))
}

test_that("Phase 2 marker maps preserve order, subsets, and absence", {
  fixture <- phase2_marker_fixture()
  expect_true(sblr:::.blr_validate_global_marker_map(fixture$map))
  expect_identical(fixture$map$global_index,
                   stats::setNames(1:4, fixture$marker_ids))

  subset_ids <- c("m4", "m2")
  subset_cross_product <-
    fixture$cross_product[subset_ids, subset_ids, drop = FALSE]
  resource <- sblr:::.blr_dense_resource(
    "subset", subset_cross_product, subset_ids,
    phase2_resource_alleles(fixture, subset_ids))
  provider <- phase2_summary_provider(
    "p_subset", "trait1", resource, fixture$marker_ids)
  collection <- sblr:::.blr_new_provider_collection(
    fixture$map, list(subset = resource), list(p_subset = provider),
    "independent_summary", "single_trait")
  expect_identical(
    sblr:::.blr_provider_presence(collection),
    matrix(c(FALSE, TRUE, FALSE, TRUE), 1L, 4L,
           dimnames = list(provider = "p_subset", marker = fixture$marker_ids)))

  bad <- provider
  bad$local_to_global[[2L]] <- bad$local_to_global[[1L]]
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(subset = resource), list(p_subset = bad),
    "independent_summary", "single_trait"), "unique resource-ordered map")
  bad <- provider
  bad$local_to_global[[1L]] <- 5L
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(subset = resource), list(p_subset = bad),
    "independent_summary", "single_trait"), "unique resource-ordered map")
  bad_resource <- resource
  bad_resource$alleles$effect[[1L]] <- "A"
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(subset = bad_resource), list(p_subset = provider),
    "independent_summary", "single_trait"), "alleles disagree")
  bad_resource <- resource
  bad_resource$genotype_coding <- "other_allele_count"
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(subset = bad_resource), list(p_subset = provider),
    "independent_summary", "single_trait"), "coding disagrees")
  bad <- provider
  dimnames(bad$sufficient_statistics$score)[[1L]] <- rev(subset_ids)
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(subset = resource), list(p_subset = bad),
    "independent_summary", "single_trait"), "score.*marker.*trait order")
})

test_that("dense CSR and block-eigen adapters represent declared operators", {
  fixture <- phase2_marker_fixture()
  alleles <- phase2_resource_alleles(fixture, fixture$marker_ids)
  dense <- sblr:::.blr_dense_resource(
    "dense", fixture$cross_product, fixture$marker_ids, alleles)
  csr <- sblr:::.blr_csr_from_dense_resource(
    "csr", fixture$cross_product, fixture$marker_ids, alleles)
  full <- sblr:::.blr_block_eigen_resource(
    "full", fixture$cross_product, fixture$marker_ids, alleles,
    blocks = list(1:2, 3:4))
  full_reverse <- sblr:::.blr_block_eigen_resource(
    "full_reverse", fixture$cross_product, fixture$marker_ids, alleles,
    blocks = list(3:4, 1:2))
  retained <- sblr:::.blr_block_eigen_resource(
    "retained", fixture$cross_product, fixture$marker_ids, alleles,
    blocks = list(1:2, 3:4), retained_ranks = c(1L, 1L))
  vector <- c(0.2, -0.3, 0.5, 0.7)

  expect_equal(sblr:::.blr_operator_apply(csr, vector),
               sblr:::.blr_operator_apply(dense, vector), tolerance = 1e-14)
  expect_equal(sblr:::.blr_operator_apply(full, vector),
               sblr:::.blr_operator_apply(dense, vector), tolerance = 1e-12)
  expect_equal(sblr:::.blr_operator_apply(full_reverse, vector),
               sblr:::.blr_operator_apply(full, vector), tolerance = 1e-14)
  expect_equal(sblr:::.blr_operator_apply(retained, vector),
               drop(sblr:::.blr_operator_matrix(retained) %*% vector),
               tolerance = 1e-14)
  expect_gt(max(abs(sblr:::.blr_operator_apply(retained, vector) -
                      fixture$cross_product %*% vector)), 1e-3)
  expect_identical(full$approximation,
                   "exact_declared_block_diagonal_operator")
  expect_identical(retained$approximation, "retained_rank_approximation")

  thresholded <- sblr:::.blr_csr_from_dense_resource(
    "thresholded", fixture$cross_product, fixture$marker_ids, alleles,
    zero_tolerance = 0.5)
  expect_false(sblr:::.blr_operator_storage_payload(
    thresholded$storage)$complete)
  expect_false(isTRUE(all.equal(
    sblr:::.blr_operator_apply(thresholded, vector),
    sblr:::.blr_operator_apply(dense, vector))))

  asymmetric <- fixture$cross_product
  asymmetric[1L, 2L] <- asymmetric[1L, 2L] + 0.1
  expect_error(sblr:::.blr_dense_resource(
    "asymmetric", asymmetric, fixture$marker_ids, alleles), "symmetric")
})

test_that("providers share one immutable resource without task state", {
  fixture <- phase2_marker_fixture()
  resource <- sblr:::.blr_dense_resource(
    "shared", fixture$cross_product, fixture$marker_ids,
    phase2_resource_alleles(fixture, fixture$marker_ids))
  providers <- list(
    p1 = phase2_summary_provider(
      "p1", "trait1", resource, fixture$marker_ids),
    p2 = phase2_summary_provider(
      "p2", "trait2", resource, fixture$marker_ids))
  collection <- sblr:::.blr_new_provider_collection(
    fixture$map, list(shared = resource), providers,
    "independent_summary", "independent_traits")

  first <- sblr:::.blr_provider_resource(collection, "p1")
  second <- sblr:::.blr_provider_resource(collection, "p2")
  expect_identical(first$storage, second$storage)
  expect_length(collection$operator_resources, 1L)
  expect_error(first$storage$payload <- diag(4), "locked binding")

  residual_1 <- numeric(4)
  residual_2 <- numeric(4)
  residual_1[[1L]] <- 1
  expect_identical(residual_2, numeric(4))
  expect_false(any(c("residual", "effect", "rng", "diagnostic") %in%
                     names(first)))
})

test_that("provider-local permutations preserve mapped likelihood operations", {
  fixture <- phase2_marker_fixture()
  score <- matrix(c(.2, -.1, .4, .3), 4L, 1L,
                  dimnames = list(marker = fixture$marker_ids,
                                  trait = "trait1"))
  permutation <- c(3L, 1L, 4L, 2L)
  permuted_ids <- fixture$marker_ids[permutation]
  effects <- matrix(c(.1, -.2, .05, .3), 4L, 1L,
                    dimnames = list(marker = fixture$marker_ids,
                                    trait = "trait1"))
  make_resource <- function(kind, id, cross_product, marker_ids) {
    alleles <- phase2_resource_alleles(fixture, marker_ids)
    switch(kind,
      dense = sblr:::.blr_dense_resource(
        id, cross_product, marker_ids, alleles),
      csr = sblr:::.blr_csr_from_dense_resource(
        id, cross_product, marker_ids, alleles),
      block_eigen = sblr:::.blr_block_eigen_resource(
        id, cross_product, marker_ids, alleles,
        blocks = list(seq_along(marker_ids))))
  }
  for (kind in c("dense", "csr", "block_eigen")) {
    original <- make_resource(
      kind, "original", fixture$cross_product, fixture$marker_ids)
    provider <- phase2_summary_provider(
      "p", "trait1", original, fixture$marker_ids, score)
    collection <- sblr:::.blr_new_provider_collection(
      fixture$map, list(original = original), list(p = provider),
      "independent_summary", "single_trait")
    permuted <- make_resource(
      kind, "permuted",
      fixture$cross_product[permutation, permutation, drop = FALSE],
      permuted_ids)
    permuted_provider <- phase2_summary_provider(
      "p_permuted", "trait1", permuted, fixture$marker_ids,
      score[permutation, , drop = FALSE])
    permuted_collection <- sblr:::.blr_new_provider_collection(
      fixture$map, list(permuted = permuted),
      list(p_permuted = permuted_provider),
      "independent_summary", "single_trait")
    expected <- sblr:::.blr_provider_residual_score(
      collection, "p", effects)
    actual <- sblr:::.blr_provider_residual_score(
      permuted_collection, "p_permuted", effects)
    expect_equal(
      actual[match(fixture$marker_ids, permuted_ids), , drop = FALSE],
      expected, tolerance = if (kind == "block_eigen") 1e-12 else 1e-14,
      info = kind)
  }
})

test_that("provider order is computational and absent markers add no score", {
  fixture <- phase2_marker_fixture()
  ids1 <- c("m1", "m3")
  ids2 <- c("m4", "m2")
  resource1 <- sblr:::.blr_dense_resource(
    "r1", fixture$cross_product[ids1, ids1, drop = FALSE], ids1,
    phase2_resource_alleles(fixture, ids1))
  resource2 <- sblr:::.blr_dense_resource(
    "r2", fixture$cross_product[ids2, ids2, drop = FALSE], ids2,
    phase2_resource_alleles(fixture, ids2))
  provider1 <- phase2_summary_provider(
    "p1", "trait1", resource1, fixture$marker_ids)
  provider2 <- phase2_summary_provider(
    "p2", "trait2", resource2, fixture$marker_ids)
  resources <- list(r1 = resource1, r2 = resource2)
  forward <- sblr:::.blr_new_provider_collection(
    fixture$map, resources, list(p1 = provider1, p2 = provider2),
    "independent_summary", "independent_traits")
  reverse <- sblr:::.blr_new_provider_collection(
    fixture$map, resources, list(p2 = provider2, p1 = provider1),
    "independent_summary", "independent_traits")
  forward_scores <- lapply(names(forward$providers), function(id) {
    provider <- forward$providers[[id]]
    sblr:::.blr_provider_global_score(
      provider, sblr:::.blr_provider_resource(forward, id), fixture$map)
  })
  reverse_scores <- lapply(names(reverse$providers), function(id) {
    provider <- reverse$providers[[id]]
    sblr:::.blr_provider_global_score(
      provider, sblr:::.blr_provider_resource(reverse, id), fixture$map)
  })
  names(forward_scores) <- names(forward$providers)
  names(reverse_scores) <- names(reverse$providers)
  expect_identical(forward_scores, reverse_scores[names(forward_scores)])
  expect_true(all(forward_scores$p1[c("m2", "m4"), , drop = FALSE] == 0))
  expect_true(all(forward_scores$p2[c("m1", "m3"), , drop = FALSE] == 0))
})

test_that("independent providers for one effect require compatible scales", {
  fixture <- phase2_marker_fixture()
  resource <- sblr:::.blr_dense_resource(
    "shared", fixture$cross_product, fixture$marker_ids,
    phase2_resource_alleles(fixture, fixture$marker_ids))
  first <- phase2_summary_provider(
    "p1", "trait1", resource, fixture$marker_ids)
  second <- phase2_summary_provider(
    "p2", "trait1", resource, fixture$marker_ids)
  expect_s3_class(sblr:::.blr_new_provider_collection(
    fixture$map, list(shared = resource), list(p1 = first, p2 = second),
    "independent_summary", "single_trait"), "blr_provider_collection_v1")
  second$effect_scale <- "incompatible_scale"
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(shared = resource), list(p1 = first, p2 = second),
    "independent_summary", "single_trait"), "compatible effect scales")
})

test_that("common-sample multi-trait BED remains one structural provider", {
  fixture <- phase2_marker_fixture()
  bed <- sblr:::.blr_new_operator_resource(
    "bed_shared", "bed", fixture$marker_ids,
    phase2_resource_alleles(fixture, fixture$marker_ids),
    "effect_allele_count", "allele_frequency_centered",
    "variance_standardized", "individual_genotypes",
    sblr:::.blr_new_operator_storage_ref(
      "packed_bed", list(
        bed_files = "fixture.bed", source_sample_count = 7L,
        sample_ids = paste0("i", 1:5), selected_rows = c(7L, 2L, 5L, 1L, 6L),
        selected_columns = list(seq_along(fixture$marker_ids)))),
    block_eigen = NULL, approximation = "exact_selected_genotypes",
    provenance = list(source = "Phase 2 structural fixture"))
  traits <- c("trait1", "trait2")
  phenotype <- matrix(seq_len(10) / 10, 5L, 2L,
                      dimnames = list(sample = paste0("i", 1:5), trait = traits))
  provider <- sblr:::.blr_new_likelihood_provider(
    "joint_bed", traits, "bed_shared",
    stats::setNames(seq_along(fixture$marker_ids), fixture$marker_ids),
    sufficient_statistics = list(phenotype = phenotype),
    sample_size = stats::setNames(c(5, 5), traits),
    likelihood_regime = "common_sample",
    residual_contract = "fixed_full_residual_covariance",
    population = "fixture_population", effect_scale = "standardized",
    overlap_group = NULL,
    provenance = list(source = "Phase 2 structural fixture"))
  collection <- sblr:::.blr_new_provider_collection(
    fixture$map, list(bed_shared = bed), list(joint_bed = provider),
    "common_sample", "joint_multitrait")
  expect_length(collection$providers, 1L)
  expect_identical(collection$providers$joint_bed$trait_ids, traits)
  expect_identical(collection$providers$joint_bed$operator_resource_id,
                   "bed_shared")
  bad <- provider
  bad$sufficient_statistics$phenotype <-
    bad$sufficient_statistics$phenotype[rev(seq_len(5L)), , drop = FALSE]
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(bed_shared = bed), list(joint_bed = bad),
    "common_sample", "joint_multitrait"), "sample order")
  singleton <- sblr:::.blr_new_likelihood_provider(
    "p1", "trait1", "bed_shared", provider$local_to_global,
    sufficient_statistics = list(
      phenotype = phenotype[, "trait1", drop = FALSE]),
    sample_size = c(trait1 = 5), likelihood_regime = "common_sample",
    residual_contract = "fixed_full_residual_covariance",
    population = "fixture_population", effect_scale = "standardized",
    overlap_group = NULL,
    provenance = list(source = "Phase 2 structural fixture"))
  expect_error(sblr:::.blr_new_provider_collection(
    fixture$map, list(bed_shared = bed),
    list(p1 = singleton), "common_sample", "joint_multitrait"),
    "non-factorized")
})

test_that("maintained resolved specifications use Phase 2 constructors", {
  chain <- sblr:::.blr_chain_controls(
    2L, 0L, 1L, 11L, 1L, 1L, NULL, FALSE)
  resolved <- sblr:::resolve_blr_spec_from_wrapper(
    "sbayesc", "csr", "trait1", c("m1", "m2"), chain,
    sample_sizes = 100)
  expect_s3_class(resolved$data$operator_resources[[1L]],
                  "blr_operator_resource_v1")
  expect_s3_class(resolved$data$providers[[1L]],
                  "blr_likelihood_provider_v1")
  expect_true(sblr:::validate_blr_resolved_spec(resolved))
  expect_identical(resolved$schema$seed_contract_version, 0L)
  expect_identical(resolved$schema$retention_contract_version, 0L)
})
