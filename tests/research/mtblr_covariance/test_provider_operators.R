source(file.path("tests", "research", "mtblr_covariance",
                 "mtblr_provider_operators.R"))

expect_close <- function(actual, expected, tolerance = 1e-10, label) {
  error <- max(abs(actual - expected))
  if (!is.finite(error) || error > tolerance) {
    stop(label, ": maximum error ", format(error, digits = 8),
         " exceeds ", tolerance, call. = FALSE)
  }
  invisible(error)
}

expect_error_pattern <- function(expression, pattern, label) {
  message <- tryCatch({
    force(expression)
    NA_character_
  }, error = conditionMessage)
  if (is.na(message) || !grepl(pattern, message, perl = TRUE)) {
    stop(label, ": expected error matching `", pattern, "`, got `",
         message, "`.", call. = FALSE)
  }
  invisible(message)
}

permuted_provider <- function(id, original, cross_product, permutation,
                              operator_kind = c("dense", "csr", "eigen")) {
  operator_kind <- match.arg(operator_kind)
  permuted_cross_product <- cross_product[permutation, permutation, drop = FALSE]
  operator <- switch(
    operator_kind,
    dense = mt_dense_operator(permuted_cross_product),
    csr = mt_csr_operator(permuted_cross_product),
    eigen = mt_eigen_operator(permuted_cross_product,
                              blocks = list(seq_along(permutation)),
                              retained_rank = length(permutation))
  )
  mt_summary_provider(
    id = id, trait = original$trait, score = original$score[permutation],
    sample_size = original$sample_size, markers = original$markers[permutation],
    global_markers = original$global_markers, operator = operator,
    effect_allele = original$effect_allele[permutation],
    global_effect_allele = original$global_effect_allele,
    genotype_coding = original$genotype_coding,
    genotype_scale = original$genotype_scale,
    phenotype_scale = original$phenotype_scale,
    population = original$population, error_group = original$error_group
  )
}

global_markers <- c("m1", "m2", "m3")
global_alleles <- c("A", "C", "G")

C1 <- matrix(c(2.2, 0.45,
               0.45, 1.4), 2, 2, byrow = TRUE)
C2 <- matrix(c(1.1, -0.2,
               -0.2, 1.8), 2, 2, byrow = TRUE)

provider_1 <- mt_summary_provider(
  id = "cohort-a-trait-1", trait = 1L, score = c(0.8, -0.3),
  sample_size = 120, markers = c("m1", "m2"),
  global_markers = global_markers, operator = mt_dense_operator(C1),
  effect_allele = global_alleles[c(1, 2)],
  global_effect_allele = global_alleles, population = "population-a"
)
provider_2 <- mt_summary_provider(
  id = "cohort-b-trait-2", trait = 2L, score = c(-0.1, 0.65),
  sample_size = 75, markers = c("m1", "m3"),
  global_markers = global_markers,
  operator = mt_eigen_operator(C2, retained_rank = 2L),
  effect_allele = global_alleles[c(1, 3)],
  global_effect_allele = global_alleles, population = "population-b"
)

effects <- matrix(c(0.20, -0.10,
                    0.05,  0.30,
                   -0.25,  0.15), nrow = 3, byrow = TRUE)

heterogeneous_value <- mt_provider_loglik(list(provider_1, provider_2), effects)
heterogeneous_reversed <- mt_provider_loglik(list(provider_2, provider_1), effects)
expect_close(heterogeneous_value, heterogeneous_reversed, 1e-14,
             "provider-order invariance")

# Provider-local marker order is arbitrary when identifiers, score, alleles,
# operator coordinates, and the induced local-to-global map move together.
permutation <- c(2L, 1L)
for (kind in c("dense", "csr", "eigen")) {
  original_operator <- switch(
    kind,
    dense = mt_dense_operator(C1),
    csr = mt_csr_operator(C1),
    eigen = mt_eigen_operator(C1, blocks = list(1:2), retained_rank = 2L)
  )
  original <- provider_1
  original$id <- paste0("local-order-original-", kind)
  original$operator <- original_operator
  permuted <- permuted_provider(paste0("local-order-permuted-", kind),
                                original, C1, permutation, kind)
  stopifnot(identical(permuted$map, original$map[permutation]),
            identical(permuted$markers, original$markers[permutation]))
  expect_close(mt_provider_loglik(list(original), effects),
               mt_provider_loglik(list(permuted), effects),
               if (kind == "eigen") 1e-12 else 1e-14,
               paste(kind, "provider-local marker-order invariance"))
}

manual_value <-
  sum(provider_1$score * effects[provider_1$map, 1]) -
  0.5 * drop(crossprod(effects[provider_1$map, 1],
                       C1 %*% effects[provider_1$map, 1])) +
  sum(provider_2$score * effects[provider_2$map, 2]) -
  0.5 * drop(crossprod(effects[provider_2$map, 2],
                       C2 %*% effects[provider_2$map, 2]))
expect_close(heterogeneous_value, manual_value, 1e-12,
             "heterogeneous-provider likelihood")

# A marker absent from a provider supplies no likelihood information.  It is
# not evidence that the corresponding effect is zero.
effects_missing_changed <- effects
effects_missing_changed[3, 1] <- 100
expect_close(mt_provider_loglik(list(provider_1), effects),
             mt_provider_loglik(list(provider_1), effects_missing_changed),
             0, "missing-marker likelihood neutrality")

# Two independent providers for the same shared trait can be accumulated
# separately or, when their maps and scales coincide, as one summed operator.
C3 <- matrix(c(0.9, 0.1,
               0.1, 1.3), 2, 2, byrow = TRUE)
provider_3 <- mt_summary_provider(
  id = "cohort-c-trait-1", trait = 1L, score = c(0.25, 0.4),
  sample_size = 63, markers = c("m1", "m2"),
  global_markers = global_markers, operator = mt_csr_operator(C3),
  effect_allele = global_alleles[c(1, 2)],
  global_effect_allele = global_alleles, population = "population-c"
)
combined <- mt_combine_same_trait_providers(list(provider_1, provider_3))
expect_close(mt_provider_loglik(list(provider_1, provider_3), effects),
             mt_provider_loglik(list(combined), effects), 1e-12,
             "separate versus combined independent providers")
stopifnot(combined$sample_size == 183)

# Dense, CSR and full-rank eigen operators represent the same matrix exactly
# up to floating-point eigendecomposition error.
dense <- mt_dense_operator(C1)
csr <- mt_csr_operator(C1)
eigen_full <- mt_eigen_operator(C1, retained_rank = 2L)
probe <- c(0.37, -0.22)
expect_close(mt_operator_apply(csr, probe),
             mt_operator_apply(dense, probe), 1e-14,
             "CSR versus dense operator")
expect_close(mt_operator_apply(eigen_full, probe),
             mt_operator_apply(dense, probe), 1e-12,
             "full-rank eigen versus dense operator")
stopifnot(isTRUE(eigen_full$exact_for_supplied_matrix))

# A retained-rank operator equals its own reconstructed approximation, not the
# original dense operator.  The difference is deliberately material here.
eigen_rank_1 <- mt_eigen_operator(C1, retained_rank = 1L)
retained_matrix <- mt_operator_matrix(eigen_rank_1)
expect_close(mt_operator_apply(eigen_rank_1, probe),
             as.numeric(retained_matrix %*% probe), 1e-14,
             "retained eigen self-consistency")
if (max(abs(mt_operator_apply(eigen_rank_1, probe) - C1 %*% probe)) < 1e-3) {
  stop("The retained-rank fixture does not expose its approximation.",
       call. = FALSE)
}
stopifnot(!isTRUE(eigen_rank_1$exact_for_supplied_matrix))

# Block traversal order is computational only and cannot change the operator.
C_block <- matrix(c(1.8, 0.35, 0,
                    0.35, 1.2, 0,
                    0, 0, 0.7), 3, 3, byrow = TRUE)
blocks_forward <- mt_eigen_operator(C_block, list(1:2, 3L))
blocks_reverse <- mt_eigen_operator(C_block, list(3L, 1:2))
block_probe <- c(-0.2, 0.4, 0.9)
expect_close(mt_operator_apply(blocks_forward, block_probe),
             mt_operator_apply(blocks_reverse, block_probe), 1e-14,
             "eigenblock-order invariance")
expect_close(mt_operator_apply(blocks_forward, block_probe),
             as.numeric(C_block %*% block_probe), 1e-12,
             "full-rank eigenblock representation")

# Residual scores are local: r_dt = s_dt - C_dt b_t on the provider map.
expect_close(mt_provider_residual_score(provider_1, effects[, 1]),
             provider_1$score - C1 %*% effects[provider_1$map, 1], 1e-14,
             "provider residual score")

# Incompatible scales/coding and undeclared sample overlap fail rather than
# being silently pooled.
provider_bad_coding <- provider_3
provider_bad_coding$genotype_coding <- "standardized_genotype"
expect_error_pattern(
  mt_provider_loglik(list(provider_1, provider_bad_coding), effects),
  "incompatible genotype_coding", "incompatible coding"
)
provider_overlap_a <- provider_1
provider_overlap_b <- provider_3
provider_overlap_a$error_group <- "shared-samples"
provider_overlap_b$error_group <- "shared-samples"
expect_error_pattern(
  mt_provider_loglik(list(provider_overlap_a, provider_overlap_b), effects),
  "overlap-aware likelihood", "unmodeled provider overlap"
)
expect_error_pattern(
  mt_summary_provider(
    id = "bad-allele", trait = 1L, score = c(0.8, -0.3),
    sample_size = 120, markers = c("m1", "m2"),
    global_markers = global_markers, operator = dense,
    effect_allele = c("T", "C"), global_effect_allele = global_alleles
  ),
  "aligned to the global effect allele", "allele alignment"
)
provider_bad_global_allele <- provider_2
provider_bad_global_allele$global_effect_allele[2] <- "T"
expect_error_pattern(
  mt_provider_loglik(list(provider_1, provider_bad_global_allele), effects),
  "same ordered global effect alleles", "global allele-reference agreement"
)

# Fixed-seed reproducibility is included even though these operator checks are
# deterministic; it protects future randomized provider fixtures.
set.seed(20260814)
random_effects_1 <- matrix(rnorm(6), 3, 2)
value_1 <- mt_provider_loglik(list(provider_1, provider_2), random_effects_1)
set.seed(20260814)
random_effects_2 <- matrix(rnorm(6), 3, 2)
value_2 <- mt_provider_loglik(list(provider_1, provider_2), random_effects_2)
stopifnot(identical(random_effects_1, random_effects_2), identical(value_1, value_2))

cat("Provider/operator reference checks passed.\n")
cat(sprintf("heterogeneous_loglik=%.12f\n", heterogeneous_value))
cat(sprintf("retained_rank_max_difference=%.12f\n",
            max(abs(mt_operator_apply(eigen_rank_1, probe) - C1 %*% probe))))
