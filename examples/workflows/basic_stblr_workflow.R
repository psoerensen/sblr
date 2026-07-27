## basic_stblr_workflow.R
##
## Canonical end-to-end ST-BLR workflow.
##
## This script shows the intended high-level workflow:
##
##   Glist + phenotype
##     -> make_summary_stats()
##     -> make_sparse_ld()
##     -> stblr_csr()
##     -> summarise_posterior()
##     -> summarise_components()
##     -> summarise_architecture()
##     -> make_credible_sets()
##
## The script is written as a compact template. Before running it, create or load
## a valid `Glist` object and phenotype matrix/vector `y`. The `Glist` object
## should contain genotype/BED paths, marker metadata, marker order, MAF/allele
## frequencies, and after `make_sparse_ld()` also sparse-LD metadata in
## `Glist$sparseLD`.

## ------------------------------------------------------------
## Setup
## ------------------------------------------------------------

## Load sblr if needed.
## library(sblr)

set.seed(1)

nthreads <- 1
niter <- 1000
nburn <- 500
seed <- 1

data_dir <- file.path(tempdir(), "sblr_basic_workflow")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## 1. Prepare Glist
## ------------------------------------------------------------

## This section should be replaced by your own data-loading/preparation code.
##
## Required objects:
##
##   Glist : package-specific genotype/LD carrier object
##   y     : phenotype vector or phenotype matrix
##
## Glist should carry, or be able to point to:
##
##   - genotype/BED paths
##   - marker IDs and marker order
##   - allele frequencies or MAF
##   - chromosome-specific marker lists, for example Glist$rsids
##   - LD marker order, for example Glist$rsidsLD
##   - after make_sparse_ld(), sparse-LD metadata in Glist$sparseLD
##
## Example placeholder:
##
##   Glist <- ...
##   y <- ...
##
## For a package example, prefer a small test fixture rather than private
## absolute paths.

if (!exists("Glist", inherits = FALSE)) {
  stop(
    "Please create or load `Glist` before running this workflow. ",
    "See package tests/examples for small data fixtures."
  )
}

if (!exists("y", inherits = FALSE)) {
  stop(
    "Please create or load phenotype vector/matrix `y` before running this workflow."
  )
}

## ------------------------------------------------------------
## 2. Make summary statistics
## ------------------------------------------------------------

stats <- make_summary_stats(
  Glist,
  y,
  nthreads = nthreads
)

## ------------------------------------------------------------
## 3. Make sparse LD
## ------------------------------------------------------------

Glist <- make_sparse_ld(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_test"),
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

## The updated Glist should now contain sparse-LD metadata, typically in:
##
##   Glist$sparseLD
##
## This updated object is passed to stblr_csr().

## ------------------------------------------------------------
## 4. Fit CSR BayesC
## ------------------------------------------------------------

fit_csr_bayesc <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesc",
  nit = niter,
  nburn = nburn,
  seed = seed
)

names(fit_csr_bayesc)
is.null(fit_csr_bayesc$chains)
is.null(fit_csr_bayesc$ld_swap_chains)


check_stblr_consistency(fit_csr_bayesc)

## ------------------------------------------------------------
## 5. Fit CSR BayesR
## ------------------------------------------------------------

fit_csr_bayesr <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesr",
  nit = niter,
  nburn = nburn,
  seed = seed
)

check_stblr_consistency(fit_csr_bayesr)

names(fit_csr_bayesr)
is.null(fit_csr_bayesr$chains)
is.null(fit_csr_bayesr$ld_swap_chains)



## ------------------------------------------------------------
## 6. Fit sampled selection_s
## ------------------------------------------------------------

## Sampled selection_s is currently supported for CSR BayesC, CSR BayesR,
## and CSR SBayesRC-style models. This example uses CSR BayesC.
##
## selection_s uses the standardized-genotype scale:
##
##   q_j(S) = h_j^(S + 1)
##   h_j    = 2 p_j (1 - p_j)

fit_csr_bayesc_s <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesc",
  estimate_selection_s = TRUE,
  selection_s_init = 0,
  selection_s_prior = c(-3, 2),
  selection_s_proposal_sd = 0.35,
  nit = niter,
  nburn = nburn,
  seed = seed
)

check_stblr_consistency(fit_csr_bayesc_s)

names(fit_csr_bayesc_s)
is.null(fit_csr_bayesc_s$chains)
is.null(fit_csr_bayesc_s$ld_swap_chains)

## Useful sampled-S outputs include:
##
##   fit_csr_bayesc_s$selection_s
##   fit_csr_bayesc_s$selection_s_trace
##   fit_csr_bayesc_s$selection_s_acceptance
##   fit_csr_bayesc_s$selection_s_sd
##   fit_csr_bayesc_s$selection_s_min
##   fit_csr_bayesc_s$selection_s_max

## ------------------------------------------------------------
## 7. Summarise posterior
## ------------------------------------------------------------

posterior_bayesc <- summarise_posterior(
  fit_csr_bayesc,
  nburn = nburn
)

posterior_bayesr <- summarise_posterior(
  fit_csr_bayesr,
  nburn = nburn
)

posterior_selection_s <- summarise_posterior(
  fit_csr_bayesc_s,
  nburn = nburn
)

## Optional plotting. Depending on implementation, this may return a ggplot
## object, a list of plots, or draw directly.
posterior_plot_bayesc <- plot_posterior(
  posterior_bayesc
)

## ------------------------------------------------------------
## 8. Summarise components
## ------------------------------------------------------------

## Component summaries are most relevant for BayesR/SBayesRC-style fits where
## posterior component probabilities are available.
##
## Common BayesR/SBayesRC fields include:
##
##   fit$component_probabilities
##   fit$dm_component_mean
##
## For BayesR, the null component is usually named component_0.
## For SBayesRC, the null component is usually named gamma_0.00.

components_bayesr <- summarise_components(
  fit_csr_bayesr
)

## ------------------------------------------------------------
## 9. Summarise architecture
## ------------------------------------------------------------

## Architecture summaries are descriptive post-processing summaries of the
## fitted posterior. They can be used to inspect how signal relates to MAF,
## prior scaling, PIP, or sampled selection_s behavior.
##
## If the function can infer MAF from the fit/Glist, maf = NULL may be enough.
## Otherwise, provide a marker-aligned MAF vector.
##
## Recommended marker alignment pattern when extracting MAF from Glist:
##
##   idx <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
##   maf_ld <- Glist$maf[[chr]][idx]
##
## Avoid `%in%` for this alignment because it does not preserve order.

architecture_selection_s <- summarise_architecture(
  fit_csr_bayesc_s,
  maf = NULL,
  min_pip = 0.01,
  top_n = 20
)

## ------------------------------------------------------------
## 10. Make credible sets
## ------------------------------------------------------------

## Credible sets are R-level helper outputs constructed from fitted model
## outputs such as fit$dm and LD information. They are not core sampler return
## objects.
##
## Depending on the current implementation, make_credible_sets() may require
## Glist, LD metadata, marker map information, or a PIP threshold/coverage.
## The call below shows the intended high-level pattern.

credible_sets_bayesr <- make_credible_sets(
  fit_csr_bayesr,
  Glist = Glist
)

## ------------------------------------------------------------
## End of workflow
## ------------------------------------------------------------

workflow_objects <- list(
  stats = stats,
  Glist = Glist,
  fit_csr_bayesc = fit_csr_bayesc,
  fit_csr_bayesr = fit_csr_bayesr,
  fit_csr_bayesc_s = fit_csr_bayesc_s,
  posterior_bayesc = posterior_bayesc,
  posterior_bayesr = posterior_bayesr,
  posterior_selection_s = posterior_selection_s,
  posterior_plot_bayesc = posterior_plot_bayesc,
  components_bayesr = components_bayesr,
  architecture_selection_s = architecture_selection_s,
  credible_sets_bayesr = credible_sets_bayesr
)

invisible(workflow_objects)
