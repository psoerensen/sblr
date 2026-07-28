## st_summary_workflow.R
##
## Canonical summary-statistics ST-BLR workflow.
##
## Demonstrates:
##   make_summary_stats() -> make_sparse_ld() -> stblr_csr()
##   summarise_posterior() -> summarise_components() -> summarise_architecture()
##   make_credible_sets() -> check_stblr_consistency()
##
## Requirements:
##   - a valid Glist object
##   - a phenotype vector/matrix y
##
## The script can either use existing objects in the workspace or read a Glist from
## SBLR_EXAMPLE_DATA_DIR/Glist_sparseLD_1k.RDS and simulate phenotypes with mtsim().

library(sblr)
source("./examples/workflows/workflow_helpers.R")

set.seed(1)
chr <- 1L
nthreads <- 1L
niter <- 1000L
nburn <- 500L
seed <- 1L

data_dir <- file.path(tempdir(), "sblr_basic_workflow")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------------
## 1. Prepare Glist and phenotype data
## -------------------------------------------------------------------------

if (!exists("Glist", inherits = FALSE)) {
  Glist <- read_example_glist()
}

if (!exists("y", inherits = FALSE)) {
  sim <- mtsim(
    Glist = Glist,
    chr = chr,
    rsids = Glist$rsidsLD[[chr]],
    nt = 3,
    n_shared = 30,
    n_specific = 10,
    h2 = c(0.4, 0.5, 0.3),
    rg = matrix(c(1.0, 0.7, 0.3, 0.7, 1.0, 0.5, 0.3, 0.5, 1.0), nrow = 3),
    re = 0,
    seed = seed
  )
  y <- as.matrix(scale(sim$y))
}

## -------------------------------------------------------------------------
## 2. Summary statistics and sparse LD
## -------------------------------------------------------------------------

stats <- make_summary_stats(
  Glist = Glist,
  y = y,
  chr = chr,
  nthreads = nthreads
)

Glist <- make_sparse_ld(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_basic"),
  chr = chr,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

## -------------------------------------------------------------------------
## 3. Fit CSR BayesC and BayesR
## -------------------------------------------------------------------------

fit_csr_bayesc <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "sbayesc",
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed
)

fit_csr_bayesr <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "sbayesr",
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed
)

check_fit(fit_csr_bayesc)
check_fit(fit_csr_bayesr)

fits <- list(
  csr_bayesc = fit_csr_bayesc,
  csr_bayesr = fit_csr_bayesr
)

fit_metadata_table(fits)
fit_field_inventory(fits)
summarise_fit_list(fits)

## -------------------------------------------------------------------------
## 4. Posterior, component, and architecture summaries
## -------------------------------------------------------------------------

posterior_bayesc <- summarise_posterior(fit_csr_bayesc, nburn = nburn)
posterior_bayesr <- summarise_posterior(fit_csr_bayesr, nburn = nburn)

head(posterior_bayesc)
head(posterior_bayesr)

## Optional interactive plot.
## plot_posterior(posterior_bayesc)

## Component summaries are most relevant for BayesR/SBayesRC-style models.
components_bayesr <- summarise_components(fit_csr_bayesr)
components_bayesr

## Architecture summaries can use marker-aligned MAF when available.
maf <- marker_aligned_maf(Glist, chr = chr)
architecture_bayesc <- summarise_architecture(fit_csr_bayesc, maf = maf)
architecture_bayesr <- summarise_architecture(fit_csr_bayesr, maf = maf)

architecture_bayesc
architecture_bayesr

## -------------------------------------------------------------------------
## 5. Credible sets from fitted PIPs and LD
## -------------------------------------------------------------------------

trait <- colnames(fit_csr_bayesr$dm)[1] %||% 1L

credible_sets_bayesr <- run_credible_sets(
  fit = fit_csr_bayesr,
  Glist = Glist,
  trait = trait,
  coverage = 0.95
)

credible_sets_bayesr$summary

## -------------------------------------------------------------------------
## 6. Objects returned by the workflow
## -------------------------------------------------------------------------

workflow_objects <- list(
  stats = stats,
  Glist = Glist,
  fits = fits,
  posterior_bayesc = posterior_bayesc,
  posterior_bayesr = posterior_bayesr,
  components_bayesr = components_bayesr,
  architecture_bayesc = architecture_bayesc,
  architecture_bayesr = architecture_bayesr,
  credible_sets_bayesr = credible_sets_bayesr
)

invisible(workflow_objects)
