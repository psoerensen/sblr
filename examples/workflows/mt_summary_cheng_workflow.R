# Public Phase 6B independent-summary Cheng MT-BayesC-Pi workflow.
# The tiny example is wrapped in a function so sourcing this file is cheap.
run_mtblr_summary_workflow <- function(operator = c("csr", "block_eigen")) {
  operator <- match.arg(operator)
  markers <- paste0("m", 1:4)
  traits <- c("trait1", "trait2")
  alleles <- data.frame(
    marker_id = markers, effect = c("A", "C", "A", "C"),
    other = c("G", "T", "G", "T"),
    coding = "effect_allele_count", stringsAsFactors = FALSE)
  cross_product <- matrix(c(
    4, .25, 0, 0, .25, 3, 0, 0,
    0, 0, 2.5, .125, 0, 0, .125, 2
  ), 4, byrow = TRUE, dimnames = list(markers, markers))

  common <- list(
    resource_id = "ld", marker_ids = markers, alleles = alleles,
    genotype_coding = "effect_allele_count", centering = "declared",
    standardization = "declared", operator_scale = "cross_product",
    provenance = list(source = "small public workflow"))
  if (operator == "csr") {
    resource <- c(common, list(
      approximation = "exact_declared_operator",
      csr = list(
        row_ptr = c(0, 2, 4, 6, 8),
        column_index = c(1L, 2L, 1L, 2L, 3L, 4L, 3L, 4L),
        values = c(4, .25, .25, 3, 2.5, .125, .125, 2),
        diagonal = diag(cross_product), complete = TRUE),
      diagonal = diag(cross_product)))
    fitter <- mtblr_csr
  } else {
    blocks <- lapply(list(1:2, 3:4), function(index) {
      decomposition <- eigen(cross_product[index, index], symmetric = TRUE)
      list(marker_ids = markers[index],
           eigenvectors = decomposition$vectors,
           eigenvalues = decomposition$values)
    })
    resource <- c(common, list(
      approximation = "exact_declared_block_diagonal_operator",
      blocks = blocks))
    fitter <- mtblr_block_eigen
  }
  providers <- lapply(seq_along(traits), function(index) list(
    provider_id = paste0("provider", index), trait_id = traits[index],
    operator_resource_id = "ld",
    score = setNames(c(.5, -.25, .35, .15) + .05 * (index - 1), markers),
    sample_size = 100 + 10 * index, residual_scale = .75 + .25 * index,
    likelihood_regime = "independent_summary", effect_scale = "standardized",
    population = paste0("population", index), provenance = list()))
  marker_covariance <- diag(c(.35, .4))
  dimnames(marker_covariance) <- list(traits, traits)

  fitter(
    providers = providers, operator_resources = list(resource),
    global_marker_ids = markers, global_alleles = alleles,
    trait_ids = traits, initial_marker_covariance = marker_covariance,
    marker_covariance_prior_degrees_of_freedom = 2.5,
    marker_covariance_prior_scale = marker_covariance,
    nburn = 5L, nit = 10L, nthin = 2L, seed = 2606L)
}

# The returned mtblr_fit exposes PIPs in dm, marker-pattern probabilities in
# activity_pattern_probabilities, Dirichlet summaries in pi_mean, and the
# authoritative marker covariance in cov_b_mean. predictions is NULL because
# independent summary providers do not contain individual-level observations.
