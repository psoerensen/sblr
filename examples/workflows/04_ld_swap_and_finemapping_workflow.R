## 04_ld_swap_and_finemapping_workflow.R
##
## LD-swap diagnostics and credible-set/fine-mapping workflow.
##
## Demonstrates:
##   - stblr_csr(..., updateLDswap = TRUE)
##   - fit$ld_swap and fit$ld_swap_chains
##   - make_credible_sets()
##   - optional extract_stblr_finemap_loci() if available
##
## Credible sets are post-processing outputs built from fit$dm and LD. They are
## not separate sampler return objects.

library(sblr)
source("./examples/workflows/00_workflow_helpers.R")

chr <- 1L
trait <- "D1"
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

stats <- make_summary_stats(Glist = Glist, y = y, chr = chr, nthreads = nthreads)
Glist <- make_sparse_ld(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_finemap"),
  chr = chr,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

## -------------------------------------------------------------------------
## 1. Fit baseline and LD-swap CSR models
## -------------------------------------------------------------------------

mh_conservative <- list(
  updateLDswap = TRUE,
  ld_swap_prob = 0.10,
  ld_swap_r2 = 0.05,
  ld_swap_moves = 5
)

mh_permissive <- list(
  updateLDswap = TRUE,
  ld_swap_prob = 0.50,
  ld_swap_r2 = 0.001,
  ld_swap_moves = 20
)

base_args <- list(
  stats = stats,
  Glist = Glist,
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed
)

fit_bayesc <- do.call(stblr_csr, c(base_args, list(method = "bayesc")))
fit_bayesc_mh <- do.call(stblr_csr, c(base_args, list(method = "bayesc"), mh_conservative))
fit_bayesr <- do.call(stblr_csr, c(base_args, list(method = "bayesr")))
fit_bayesr_mh <- do.call(stblr_csr, c(base_args, list(method = "bayesr"), mh_permissive))

fits <- list(
  bayesC = fit_bayesc,
  bayesC_LDswap = fit_bayesc_mh,
  bayesR = fit_bayesr,
  bayesR_LDswap = fit_bayesr_mh
)

check_fit(fit_bayesc)
check_fit(fit_bayesc_mh, require_ld_swap = TRUE)
check_fit(fit_bayesr)
check_fit(fit_bayesr_mh, require_ld_swap = TRUE)

## -------------------------------------------------------------------------
## 2. LD-swap diagnostics
## -------------------------------------------------------------------------

ld_swap_summary_table(fits)

## -------------------------------------------------------------------------
## 3. Credible sets
## -------------------------------------------------------------------------

credible_sets <- lapply(fits, function(fit) {
  run_credible_sets(
    fit = fit,
    Glist = Glist,
    trait = trait,
    coverage = 0.95
  )
})

lapply(credible_sets, function(cs) cs$summary)

## A single locus can yield multiple CS1, CS2, ... sets because
## make_credible_sets() delegates to the per-locus builder and does not expose
## max_sets. Under-coverage attempts are dropped by default when
## allow_incomplete = FALSE.

## -------------------------------------------------------------------------
## 4. Optional fine-mapping extraction helper
## -------------------------------------------------------------------------

if (exists("extract_stblr_finemap_loci")) {
  finemap <- lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    cs <- credible_sets[[model_name]]
    out <- extract_stblr_finemap_loci(
      fit = fit,
      Glist = Glist,
      locus_sets = cs$locus_sets,
      trait = trait,
      credible_sets = TRUE,
      coverage = 0.95,
      min_r2 = 0.5,
      pip_cutoff = 0.001,
      cs_mode = "multi",
      min_signal_pip = 0.05
    )
    out$model <- model_name
    out
  })
  names(finemap) <- names(fits)
  lapply(finemap, function(x) x$summary)
}

## -------------------------------------------------------------------------
## 5. Top markers for comparison
## -------------------------------------------------------------------------

lapply(names(fits), function(model_name) {
  top_markers(fits[[model_name]], model_name = model_name, trait = trait, n = 20)
})

invisible(list(fits = fits, credible_sets = credible_sets))
