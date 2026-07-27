## 05_selection_s_workflow.R
##
## Fixed and sampled selection_s examples.
##
## selection_s is currently supported for CSR BayesC, CSR BayesR, and
## SBayesRC-style CSR models. It is not supported by the prior, learned, or
## group annotation-aware BayesC backends.
##
## On the standardized-genotype scale:
##   q_j(S) = h_j^(S + 1),  h_j = 2 p_j (1 - p_j)

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

stats <- make_summary_stats(Glist = Glist, y = y, chr = chr, nthreads = nthreads)
Glist <- make_sparse_ld(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_selection_s"),
  chr = chr,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

base_args <- list(
  stats = stats,
  Glist = Glist,
  nit = niter,
  nburn = nburn,
  ncores = nthreads,
  seed = seed
)

## -------------------------------------------------------------------------
## 1. BayesC fixed selection_s sanity checks
## -------------------------------------------------------------------------

fit_c_null <- do.call(stblr_csr, c(base_args, list(method = "bayesc", selection_s = NULL)))
fit_c_minus1 <- do.call(stblr_csr, c(base_args, list(method = "bayesc", selection_s = -1)))

## selection_s = -1 should reproduce the ordinary model because exponent = 0.
max(abs(fit_c_null$dm - fit_c_minus1$dm))
max(abs(fit_c_null$bm - fit_c_minus1$bm))
max(abs(fit_c_null$vle - fit_c_minus1$vle))
max(abs(fit_c_null$vld - fit_c_minus1$vld))

fit_c_s0 <- do.call(stblr_csr, c(base_args, list(method = "bayesc", selection_s = 0)))
fit_c_sneg05 <- do.call(stblr_csr, c(base_args, list(method = "bayesc", selection_s = -0.5)))
fit_c_spos1 <- do.call(stblr_csr, c(base_args, list(method = "bayesc", selection_s = 1)))

## -------------------------------------------------------------------------
## 2. BayesR fixed selection_s sanity checks
## -------------------------------------------------------------------------

fit_r_null <- do.call(stblr_csr, c(base_args, list(method = "bayesr", selection_s = NULL)))
fit_r_minus1 <- do.call(stblr_csr, c(base_args, list(method = "bayesr", selection_s = -1)))

max(abs(fit_r_null$dm - fit_r_minus1$dm))
max(abs(fit_r_null$bm - fit_r_minus1$bm))
max(abs(fit_r_null$vle - fit_r_minus1$vle))
max(abs(fit_r_null$vld - fit_r_minus1$vld))
max(abs(fit_r_null$dm_component_mean - fit_r_minus1$dm_component_mean))

fit_r_s0 <- do.call(stblr_csr, c(base_args, list(method = "bayesr", selection_s = 0)))
fit_r_sneg05 <- do.call(stblr_csr, c(base_args, list(method = "bayesr", selection_s = -0.5)))

## -------------------------------------------------------------------------
## 3. Sampled selection_s for CSR BayesC
## -------------------------------------------------------------------------

fit_c_sampled_s <- do.call(
  stblr_csr,
  c(
    base_args,
    list(
      method = "bayesc",
      estimate_selection_s = TRUE,
      selection_s_init = 0,
      selection_s_prior = c(-3, 2),
      selection_s_proposal_sd = 0.35
    )
  )
)

fits <- list(
  bayesC = fit_c_null,
  bayesC_s0 = fit_c_s0,
  bayesC_sneg05 = fit_c_sneg05,
  bayesC_spos1 = fit_c_spos1,
  bayesC_sampled_s = fit_c_sampled_s,
  bayesR = fit_r_null,
  bayesR_s0 = fit_r_s0,
  bayesR_sneg05 = fit_r_sneg05
)

invisible(lapply(fits, check_fit))

## -------------------------------------------------------------------------
## 4. Inspect selection_s metadata and traces
## -------------------------------------------------------------------------

lapply(fits, function(fit) {
  fit$input[c(
    "selection_s", "selection_s_fixed", "selection_s_exponent",
    "selection_s_scale", "estimate_selection_s"
  )]
})

fit_c_sampled_s[c(
  "selection_s", "selection_s_trace", "selection_s_acceptance",
  "selection_s_sd", "selection_s_min", "selection_s_max"
)]

## -------------------------------------------------------------------------
## 5. Architecture summaries using marker-aligned MAF
## -------------------------------------------------------------------------

maf <- marker_aligned_maf(Glist, chr = chr)

architecture <- do.call(rbind, lapply(names(fits), function(model_name) {
  out <- summarise_architecture(fit = fits[[model_name]], maf = maf)
  out$model <- model_name
  out[, c("model", setdiff(names(out), "model")), drop = FALSE]
}))

architecture

## -------------------------------------------------------------------------
## 6. Optional causal-marker architecture summaries for simulated data
## -------------------------------------------------------------------------

if (!is.null(sim$causal)) {
  causal_shared <- sim$causal$shared
  causal_specific <- sim$causal$specific
  causal_by_trait <- lapply(names(causal_specific), function(trait) {
    unique(c(causal_shared, causal_specific[[trait]]))
  })
  names(causal_by_trait) <- names(causal_specific)

  architecture_causal <- do.call(rbind, lapply(names(fits), function(model_name) {
    out <- summarise_architecture(
      fit = fits[[model_name]],
      maf = maf,
      markers = causal_by_trait
    )
    out$model <- model_name
    out[, c("model", setdiff(names(out), "model")), drop = FALSE]
  }))
  architecture_causal
}

invisible(fits)
