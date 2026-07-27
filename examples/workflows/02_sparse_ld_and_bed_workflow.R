## 02_sparse_ld_and_bed_workflow.R
##
## Compares summary-statistics/sparse-LD CSR models with individual-level BED
## models using the public wrappers stblr_csr() and stblr_bed().
##
## Demonstrates:
##   make_summary_stats(), make_sparse_ld(), stblr_csr(), stblr_bed(),
##   check_stblr_consistency(), summarise_components(), and basic posterior
##   comparisons.

library(sblr)
source("./examples/workflows/00_workflow_helpers.R")

chr <- 1L
nthreads <- 4L
niter <- 1000L
nburn <- 100L
seed <- 10L

data_dir <- workflow_data_dir()
Glist <- read_example_glist()

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
  seed = 1
)
y <- as.matrix(scale(sim$y))

## -------------------------------------------------------------------------
## 1. Summary statistics and sparse LD
## -------------------------------------------------------------------------

stats <- make_summary_stats(
  Glist = Glist,
  y = y,
  chr = chr,
  nthreads = nthreads
)

Glist <- make_sparse_ld(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_sparse_bed"),
  chr = chr,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

## -------------------------------------------------------------------------
## 2. Fit summary-statistics CSR models
## -------------------------------------------------------------------------

fit_csr_bayesc <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesc",
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed
)

fit_csr_bayesr <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesr",
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed
)

## -------------------------------------------------------------------------
## 3. Fit individual-level BED models
## -------------------------------------------------------------------------
## stblr_bed() uses the BED/genotype information carried by Glist. For subset
## phenotypes, use rownames(y) or names(y) that match Glist$ids.

fit_bed_bayesc <- stblr_bed(
  y = y,
  Glist = Glist,
  method = "bayesc",
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed,
  nchains = 1
)

fit_bed_bayesr <- stblr_bed(
  y = y,
  Glist = Glist,
  method = "bayesr",
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed,
  nchains = 1
)

fits <- list(
  csr_bayesc = fit_csr_bayesc,
  csr_bayesr = fit_csr_bayesr,
  bed_bayesc = fit_bed_bayesc,
  bed_bayesr = fit_bed_bayesr
)

invisible(lapply(fits, check_fit))

fit_metadata_table(fits)
fit_field_inventory(fits)
summarise_fit_list(fits)

## -------------------------------------------------------------------------
## 4. Compare PIP summaries across data levels
## -------------------------------------------------------------------------

mean_pip <- do.call(rbind, lapply(fits, function(fit) colMeans(as.matrix(fit$dm))))
max_pip <- do.call(rbind, lapply(fits, function(fit) apply(as.matrix(fit$dm), 2, max)))

mean_pip
max_pip

## Component summaries for BayesR models.
summarise_components(fit_csr_bayesr)
summarise_components(fit_bed_bayesr)

## -------------------------------------------------------------------------
## 5. Optional heavy benchmark-scale run
## -------------------------------------------------------------------------

if (workflow_should_run_heavy()) {
  fit_bed_bayesr_2c <- stblr_bed(
    y = y,
    Glist = Glist,
    method = "bayesr",
    nit = 2000,
    nburn = 500,
    ncores = nthreads,
    nchains = 2,
    chain_seeds = c(101, 202),
    keep_chains = TRUE,
    seed = seed
  )
  check_fit(fit_bed_bayesr_2c, require_chain_summaries = TRUE)
}

invisible(fits)
