## st_annotation_workflow.R
##
## Annotation-informed ST-BLR workflow.
##
## Demonstrates:
##   - annotation-unaware CSR BayesC/BayesR
##   - fixed annotation-informed BayesC priors
##   - learned annotation effects
##   - grouped BayesC priors
##   - SBayesRC-style annotation-dependent BayesR component probabilities
##   - posterior/component/annotation summaries
##
## This file is intentionally focused on model fitting and summaries. LD-swap,
## fine-mapping, multi-chain diagnostics, and maf_effect_s have separate workflows.

library(sblr)
source("./examples/workflows/workflow_helpers.R")

chr <- 1L
nthreads <- 4L
seed <- 10L
niter <- 1000L
nburn <- 100L
nthin <- 1L

data_dir <- workflow_data_dir()
Glist <- read_example_glist()

## -------------------------------------------------------------------------
## 1. Simulate annotated multi-trait data
## -------------------------------------------------------------------------
## mtsim_annotation() simulates overlapping annotations, enriched causal-marker
## sampling, and marker-specific prior values. Replace this block with study
## phenotypes and biological annotations for real analyses.

sim <- mtsim_annotation(
  Glist = Glist,
  chr = chr,
  rsids = Glist$rsidsLD[[chr]],
  nt = 3,
  n_shared = 30,
  n_specific = 10,
  n_annotations = 5,
  annotation_prob = 0.1,
  enriched_annotations = c(1, 2),
  annotation_enrichment = 5,
  base_pi = 0.001,
  enriched_pi_multiplier = 3,
  enriched_vb_multiplier = 1.5,
  h2 = c(0.4, 0.5, 0.3),
  rg = matrix(c(1.0, 0.7, 0.3, 0.7, 1.0, 0.5, 0.3, 0.5, 1.0), nrow = 3),
  re = 0,
  seed = 1
)

y <- as.matrix(scale(sim$y))

if (exists("summarize_annotation_signal")) {
  summarize_annotation_signal(sim)
}

## -------------------------------------------------------------------------
## 2. Summary statistics, sparse LD, and annotation alignment
## -------------------------------------------------------------------------

stats <- make_summary_stats(
  Glist = Glist,
  y = y,
  chr = chr,
  nthreads = nthreads
)

Glist <- make_sparse_ld(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_annotation"),
  chr = chr,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

marker_id <- Glist$rsidsLD[[chr]]
A <- sim$annot
rownames(A) <- marker_id
stopifnot(nrow(A) == length(marker_id), identical(rownames(A), marker_id))

group <- ifelse(A[, 1] != 0, "annotated", "background")
names(group) <- rownames(A)

## -------------------------------------------------------------------------
## 3. Shared model settings
## -------------------------------------------------------------------------

base_args <- list(
  stats = stats,
  Glist = Glist,
  nit = niter,
  nburn = nburn,
  nthin = nthin,
  ncores = nthreads,
  seed = seed
)

prior_annotations <- list(
  A = A,
  fixed_pi_marker = sim$pi_marker,
  fixed_vb_multiplier = sim$vb_multiplier,
  use_pi_marker = TRUE,
  use_vb_multiplier = TRUE
)

learned_args <- list(
  annotations = A,
  annotation_model = "learned_logistic",
  learn_pi_annot = TRUE,
  learn_vb_annot = TRUE,
  rw_sd_eta_pi = 0.02,
  rw_sd_eta_vb = 0.02,
  annot_update_every = 10
)

group_args <- list(
  annotations = group,
  annotation_model = "group",
  group_names = c("annotated", "background"),
  group_pi_init = c(0.002, 0.001),
  group_vb_multiplier_init = c(1.1, 1.0),
  updatePi = TRUE,
  updateGroupVb = TRUE
)

sbayesrc_args <- list(
  annotations = A,
  annotation_model = "annotation_probit_stick",
  mixture_var = c(0, 0.01, 0.1, 1),
  annotation_intercept_prior = list(
    distribution = "normal", mean = "initial_mixture", sd = 1)
)

## -------------------------------------------------------------------------
## 4. Fit annotation-unaware and annotation-aware models
## -------------------------------------------------------------------------

fit_bayesc <- do.call(stblr_csr, c(base_args, list(method = "sbayesc")))
fit_bayesr <- do.call(stblr_csr, c(base_args, list(method = "sbayesr")))

fit_prior <- do.call(
  stblr_csr_annot,
  c(base_args, list(annotations = prior_annotations,
                    annotation_model = "fixed_marker"))
)

fit_learned <- do.call(stblr_csr_annot, c(base_args, learned_args))
fit_group <- do.call(stblr_csr_annot, c(base_args, group_args))
fit_sbayesrc <- do.call(stblr_csr_annot, c(base_args, sbayesrc_args))

fits <- list(
  bayesC = fit_bayesc,
  bayesR = fit_bayesr,
  prior = fit_prior,
  learned = fit_learned,
  group = fit_group,
  sbayesrc = fit_sbayesrc
)

invisible(lapply(fits, check_fit))

## -------------------------------------------------------------------------
## 5. Model inventory and posterior summaries
## -------------------------------------------------------------------------

fit_metadata_table(fits)
fit_field_inventory(fits)
summarise_fit_list(fits)

posterior_summaries <- lapply(fits, summarise_posterior)
lapply(posterior_summaries, head)

component_summaries <- lapply(
  c("bayesr", "sbayesrc"),
  function(model_name) {
    out <- summarise_components(fits[[model_name]])
    out$model <- model_name
    out[, c("model", setdiff(names(out), "model")), drop = FALSE]
  }
)
component_summaries <- do.call(rbind, component_summaries)
component_summaries

## -------------------------------------------------------------------------
## 6. Annotation-specific outputs
## -------------------------------------------------------------------------

annotation_summaries <- lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  if (is.null(fit$annotation_summary)) {
    return(NULL)
  }
  out <- as.data.frame(fit$annotation_summary)
  out$model <- model_name
  out[, c("model", setdiff(names(out), "model")), drop = FALSE]
})
annotation_summaries <- Filter(Negate(is.null), annotation_summaries)
annotation_summaries

lapply(fits, compact_input)

## -------------------------------------------------------------------------
## 7. Causal recovery summaries for simulated examples
## -------------------------------------------------------------------------

if (!is.null(sim$causal)) {
  causal_shared <- sim$causal$shared
  causal_specific <- sim$causal$specific
  causal_by_trait <- lapply(names(causal_specific), function(trait) {
    unique(c(causal_shared, causal_specific[[trait]]))
  })
  names(causal_by_trait) <- names(causal_specific)

  causal_topn_summary <- function(fit, model_name, causal_by_trait, top_n = c(20, 50, 100)) {
    dm <- extract_posterior(fit, "pips")
    marker <- rownames(dm)
    trait_names <- colnames(dm) %||% names(causal_by_trait)
    out <- list()
    for (trait in trait_names) {
      causal <- intersect(causal_by_trait[[trait]], marker)
      if (!length(causal)) next
      pip <- dm[, trait]
      names(pip) <- marker
      ranked_marker <- names(sort(pip, decreasing = TRUE))
      out[[trait]] <- do.call(rbind, lapply(top_n, function(n) {
        top_marker <- ranked_marker[seq_len(min(n, length(ranked_marker)))]
        n_detected <- length(intersect(causal, top_marker))
        data.frame(
          model = model_name,
          trait = trait,
          top_n = n,
          n_causal = length(causal),
          n_detected = n_detected,
          power = n_detected / length(causal),
          precision = n_detected / length(top_marker),
          stringsAsFactors = FALSE
        )
      }))
    }
    do.call(rbind, out)
  }

  topn_power <- do.call(rbind, lapply(names(fits), function(model_name) {
    causal_topn_summary(fits[[model_name]], model_name, causal_by_trait)
  }))
  topn_power
}

invisible(fits)
