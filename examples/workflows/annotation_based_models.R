# Annotation-based ST-BLR workflow
#
# Fixed-prior, learned annotation, group annotation, and continuous
# overlapping-annotation SBayesRC models are runnable through package wrappers.
#
# The MCMC settings shown here are demonstration settings. Real analyses need
# longer chains and appropriate convergence and posterior predictive checks.

# Packages -----------------------------------------------------------------

library(sblr)

# Data directory setup -----------------------------------------------------

# Data setup ---------------------------------------------------------------
data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR")
if (!nzchar(data_dir)) {
  data_dir <- "C:/Users/au223366/Documents/GitHub/examples/human"
}

nthreads <- 4

Glist <- readRDS(file.path(data_dir, "Glist_sparseLD_1k.RDS"))

# Simulate or load multi-trait phenotype data ------------------------------

# mtsim_annotation() simulates overlapping annotations, enriches causal-marker
# sampling in selected annotations, and creates marker-specific prior values.
# Replace this simulation with phenotype and biological annotation data from
# the study if desired.
sim <- mtsim_annotation(
 Glist = Glist,
 chr = 1,
 rsids = Glist$rsidsLD[[1]],
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
 rg = matrix(
  c(1.0, 0.7, 0.3, 0.7, 1.0, 0.5, 0.3, 0.5, 1.0),
  nrow = 3,
  byrow = TRUE
 ),
 re = 0,
 seed = 1
)
y <- as.matrix(scale(sim$y))

# Inspect the simulated annotation enrichment before model fitting.
summarize_annotation_signal(sim)

# Compute BED sufficient statistics with bed_xtx_xty() --------------------

# Compute summary statistics
stats <- make_stats(
  Glist,
  y,
  nthreads = nthreads
)

# Compute sparse LD
Glist <- make_sparseLD(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_test"),
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = 4
)


# Use the simulated annotation matrix A -----------------------------------

# A is an m x K marker annotation matrix. Rows of A must correspond exactly
# to the marker ordering in stats and the sparse-LD files. Real analyses should
# replace sim$annot with biological annotations in the same marker order.
# A is an m x K marker annotation matrix. Rows of A must correspond exactly
# to the marker ordering in stats and the sparse-LD files. Real analyses should
# replace sim$annot with biological annotations in the same marker order.

marker_id <- Glist$rsidsLD[[1]]
m <- length(marker_id)

A <- sim$annot
rownames(A) <- marker_id

stopifnot(
  nrow(A) == m,
  identical(rownames(A), marker_id)
)

# Annotation-based ST-BLR workflow
#
# Compares annotation-unaware CSR BayesC/BayesR models with annotation-aware
# CSR models: fixed-prior, learned annotation, group annotation, and SBayesRC.
#
# Includes optional LD-swap/MH variants for models where implemented.
#
# Demonstration settings only. Real analyses need longer chains,
# convergence checks, and posterior predictive diagnostics.

# Packages -----------------------------------------------------------------

library(sblr)

# Data setup ---------------------------------------------------------------

data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR")
if (!nzchar(data_dir)) {
  data_dir <- "C:/Users/au223366/Documents/GitHub/examples/human"
}

chr <- 1L
nthreads <- 4L
seed <- 10L

nit <- 1000L
nburn <- 100L
nthin <- 1L

Glist <- readRDS(file.path(data_dir, "Glist_sparseLD_1k.RDS"))

# Simulate or load multi-trait phenotype data ------------------------------

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
  rg = matrix(
    c(
      1.0, 0.7, 0.3,
      0.7, 1.0, 0.5,
      0.3, 0.5, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  ),
  re = 0,
  seed = 1
)

y <- as.matrix(scale(sim$y))

# Inspect simulated annotation enrichment before fitting.
summarize_annotation_signal(sim)

# Compute summary statistics -----------------------------------------------

stats <- make_stats(
  Glist = Glist,
  y = y,
  nthreads = nthreads
)

# Compute sparse LD --------------------------------------------------------

ld_prefix <- file.path(data_dir, "ld_test")

Glist <- make_sparseLD(
  Glist = Glist,
  out_prefix = ld_prefix,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

stopifnot(!is.null(Glist$sparseLD$prefix))

# Prepare annotation matrix ------------------------------------------------

marker_id <- Glist$rsidsLD[[chr]]
m <- length(marker_id)

A <- sim$annot
rownames(A) <- marker_id

stopifnot(
  nrow(A) == m,
  identical(rownames(A), marker_id)
)

# Group annotation from the first simulated annotation ---------------------

group <- ifelse(A[, 1] != 0, "annotated", "background")
names(group) <- rownames(A)

# Shared model settings ----------------------------------------------------

base_args <- list(
  stats = stats,
  Glist = Glist,
  nit = nit,
  nburn = nburn,
  nthin = nthin,
  seed = seed
)

# -------------------------------------------------------------------------
# Annotation-unaware CSR models
# -------------------------------------------------------------------------

fitC <- do.call(
  stblr_csr,
  c(
    base_args,
    list(method = "bayesC")
  )
)

fitC_MH <- do.call(
  stblr_csr,
  c(
    base_args,
    list(
      method = "bayesC",
      updateLDswap = TRUE,
      ld_swap_prob = 0.10,
      ld_swap_r2 = 0.05,
      ld_swap_moves = 5
    )
  )
)

fitR <- do.call(
  stblr_csr,
  c(
    base_args,
    list(method = "bayesR")
  )
)

fitR_MH <- do.call(
  stblr_csr,
  c(
    base_args,
    list(
      method = "bayesR",
      updateLDswap = TRUE,
      ld_swap_prob = 0.50,
      ld_swap_r2 = 0.001,
      ld_swap_moves = 20
    )
  )
)

# -------------------------------------------------------------------------
# Annotation-aware CSR models
# -------------------------------------------------------------------------

fit_prior <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = list(
        A = A,
        fixed_pi_marker = sim$pi_marker,
        fixed_vb_multiplier = sim$vb_multiplier,
        use_pi_marker = TRUE,
        use_vb_multiplier = TRUE
      ),
      annotation_model = "prior"
    )
  )
)

fit_learned <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = A,
      annotation_model = "learned",
      learn_pi_annot = TRUE,
      learn_vb_annot = TRUE,
      rw_sd_eta_pi = 0.02,
      rw_sd_eta_vb = 0.02,
      annot_update_every = 10
    )
  )
)

fit_group <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = group,
      annotation_model = "group",
      group_names = c("annotated", "background"),
      group_pi_init = c(0.002, 0.001),
      group_vb_multiplier_init = c(1.1, 1.0),
      updatePi = TRUE,
      updateGroupVb = TRUE
    )
  )
)

fit_sbayesrc <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = A,
      annotation_model = "sbayesrc",
      gamma = c(0, 0.01, 0.1, 1)
    )
  )
)

# -------------------------------------------------------------------------
# Annotation-aware CSR models with LD-swap/MH
# -------------------------------------------------------------------------
# These calls require LD-swap/MH support to be implemented for each
# annotation-aware backend.

fit_prior_MH <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = list(
        A = A,
        fixed_pi_marker = sim$pi_marker,
        fixed_vb_multiplier = sim$vb_multiplier,
        use_pi_marker = TRUE,
        use_vb_multiplier = TRUE
      ),
      annotation_model = "prior",
      updateLDswap = TRUE,
      ld_swap_prob = 0.10,
      ld_swap_r2 = 0.05,
      ld_swap_moves = 5
    )
  )
)

fit_group_MH <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = group,
      annotation_model = "group",
      group_names = c("annotated", "background"),
      group_pi_init = c(0.002, 0.001),
      group_vb_multiplier_init = c(1.1, 1.0),
      updatePi = TRUE,
      updateGroupVb = TRUE,
      updateLDswap = TRUE,
      ld_swap_prob = 0.10,
      ld_swap_r2 = 0.05,
      ld_swap_moves = 5
    )
  )
)

fit_learned_MH <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = A,
      annotation_model = "learned",
      learn_pi_annot = TRUE,
      learn_vb_annot = TRUE,
      rw_sd_eta_pi = 0.02,
      rw_sd_eta_vb = 0.02,
      annot_update_every = 10,
      updateLDswap = TRUE,
      ld_swap_prob = 0.10,
      ld_swap_r2 = 0.05,
      ld_swap_moves = 5
    )
  )
)

fit_sbayesrc_MH <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = A,
      annotation_model = "sbayesrc",
      gamma = c(0, 0.01, 0.1, 1),
      updateLDswap = TRUE,
      ld_swap_prob = 0.50,
      ld_swap_r2 = 0.001,
      ld_swap_moves = 20
    )
  )
)

# -------------------------------------------------------------------------
# Optional multi-chain examples
# -------------------------------------------------------------------------
# These are shorter demonstration runs. Increase nit/nburn for real analyses.

fit_sbayesrc_2c <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = A,
  annotation_model = "sbayesrc",
  gamma = c(0, 0.01, 0.1, 1),
  nit = 200,
  nburn = 50,
  nthin = 1,
  nchains = 2,
  chain_seeds = c(10, 20),
  keep_chains = TRUE
)

fit_prior_MH_2c <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = list(
    A = A,
    fixed_pi_marker = sim$pi_marker,
    fixed_vb_multiplier = sim$vb_multiplier,
    use_pi_marker = TRUE,
    use_vb_multiplier = TRUE
  ),
  annotation_model = "prior",
  updateLDswap = TRUE,
  ld_swap_prob = 0.10,
  ld_swap_r2 = 0.05,
  ld_swap_moves = 5,
  nit = 200,
  nburn = 50,
  nthin = 1,
  nchains = 2,
  chain_seeds = c(10, 20),
  keep_chains = TRUE
)

# -------------------------------------------------------------------------
# Collect fits
# -------------------------------------------------------------------------

fits <- list(
  bayesC        = fitC,
  bayesC_MH     = fitC_MH,
  bayesR        = fitR,
  bayesR_MH     = fitR_MH,
  prior         = fit_prior,
  prior_MH      = fit_prior_MH,
  learned       = fit_learned,
  learned_MH    = fit_learned_MH,
  group         = fit_group,
  group_MH      = fit_group_MH,
  sbayesrc      = fit_sbayesrc,
  sbayesrc_MH   = fit_sbayesrc_MH
)

# -------------------------------------------------------------------------
# Compact metadata
# -------------------------------------------------------------------------

model_metadata <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    x <- fits[[model_name]]$input
    
    data.frame(
      model_name = model_name,
      method = ifelse(is.null(x$method), NA, x$method),
      model = ifelse(is.null(x$model), NA, x$model),
      backend = ifelse(is.null(x$backend), NA, x$backend),
      data_level = ifelse(is.null(x$data_level), NA, x$data_level),
      annotation_model = ifelse(is.null(x$annotation_model), NA, x$annotation_model),
      annotations = ifelse(is.null(x$annotations), FALSE, x$annotations),
      updateLDswap = ifelse(is.null(x$updateLDswap), FALSE, x$updateLDswap),
      n_markers = nrow(fits[[model_name]]$dm),
      n_traits = ncol(fits[[model_name]]$dm),
      stringsAsFactors = FALSE
    )
  })
)

model_metadata

# -------------------------------------------------------------------------
# Compact fit summaries
# -------------------------------------------------------------------------

summarise_fit <- function(fit, model_name) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)
  
  trait_names <- colnames(dm)
  if (is.null(trait_names)) {
    trait_names <- paste0("T", seq_len(ncol(dm)))
  }
  
  data.frame(
    model = model_name,
    trait = trait_names,
    n_markers = nrow(dm),
    mean_pip = colMeans(dm, na.rm = TRUE),
    sum_pip = colSums(dm, na.rm = TRUE),
    max_pip = apply(dm, 2, max, na.rm = TRUE),
    n_pip_gt_0_001 = colSums(dm > 0.001, na.rm = TRUE),
    n_pip_gt_0_01 = colSums(dm > 0.01, na.rm = TRUE),
    n_pip_gt_0_05 = colSums(dm > 0.05, na.rm = TRUE),
    n_pip_gt_0_5 = colSums(dm > 0.5, na.rm = TRUE),
    n_pip_gt_0_95 = colSums(dm > 0.95, na.rm = TRUE),
    mean_abs_bm = colMeans(abs(bm), na.rm = TRUE),
    max_abs_bm = apply(abs(bm), 2, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

fit_summary <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    summarise_fit(fits[[model_name]], model_name)
  })
)

fit_summary

sum_pip_matrix <- do.call(
  rbind,
  lapply(fits, function(fit) colSums(as.matrix(fit$dm), na.rm = TRUE))
)

sum_pip_matrix

# -------------------------------------------------------------------------
# LD-swap diagnostics
# -------------------------------------------------------------------------

ld_swap_summary <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    x <- fits[[model_name]]$ld_swap
    if (is.null(x)) return(NULL)
    
    x <- as.data.frame(x)
    x$model <- model_name
    x$trait <- rownames(x)
    
    x[, c("model", "trait", setdiff(names(x), c("model", "trait"))), drop = FALSE]
  })
)

ld_swap_summary

# -------------------------------------------------------------------------
# Top marker summaries
# -------------------------------------------------------------------------

top_markers <- function(fit, model_name, trait = 1, top_n = 20) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)
  
  trait_index <- if (is.character(trait)) match(trait, colnames(dm)) else trait
  
  if (length(trait_index) != 1 || is.na(trait_index)) {
    stop("trait must be a valid column name or column index.")
  }
  
  marker <- rownames(dm)
  if (is.null(marker)) marker <- paste0("V", seq_len(nrow(dm)))
  
  trait_name <- colnames(dm)[trait_index]
  if (is.null(trait_name)) trait_name <- paste0("T", trait_index)
  
  out <- data.frame(
    model = model_name,
    trait = trait_name,
    marker = marker,
    pip = dm[, trait_index],
    bm = bm[, trait_index],
    abs_bm = abs(bm[, trait_index]),
    stringsAsFactors = FALSE
  )
  
  out <- out[order(-out$pip, -out$abs_bm), , drop = FALSE]
  utils::head(out, top_n)
}

top_trait1 <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    top_markers(fits[[model_name]], model_name, trait = 1, top_n = 20)
  })
)

top_trait1

# -------------------------------------------------------------------------
# Top-marker overlap
# -------------------------------------------------------------------------

top_overlap_matrix <- function(fits, trait = 1, top_n = 100) {
  top_sets <- lapply(fits, function(fit) {
    dm <- as.matrix(fit$dm)
    marker <- rownames(dm)
    if (is.null(marker)) marker <- paste0("V", seq_len(nrow(dm)))
    
    marker[order(-dm[, trait])[seq_len(top_n)]]
  })
  
  out <- outer(
    names(top_sets),
    names(top_sets),
    Vectorize(function(a, b) {
      length(intersect(top_sets[[a]], top_sets[[b]]))
    })
  )
  
  dimnames(out) <- list(names(top_sets), names(top_sets))
  out
}

top20_overlap_trait1 <- top_overlap_matrix(fits, trait = 1, top_n = 20)
top100_overlap_trait1 <- top_overlap_matrix(fits, trait = 1, top_n = 100)
top500_overlap_trait1 <- top_overlap_matrix(fits, trait = 1, top_n = 500)

top20_overlap_trait1
top100_overlap_trait1
top500_overlap_trait1

# -------------------------------------------------------------------------
# PIP correlations
# -------------------------------------------------------------------------

pip_correlation_signal_markers <- function(fits, trait = 1,
                                           pip_threshold = 0.001,
                                           method = "spearman") {
  pip_mat <- do.call(
    cbind,
    lapply(fits, function(fit) {
      as.matrix(fit$dm)[, trait]
    })
  )
  
  colnames(pip_mat) <- names(fits)
  
  keep <- apply(pip_mat, 1, max, na.rm = TRUE) > pip_threshold
  
  stats::cor(
    pip_mat[keep, , drop = FALSE],
    use = "pairwise.complete.obs",
    method = method
  )
}

pip_cor_signal_trait1 <- pip_correlation_signal_markers(
  fits,
  trait = 1,
  pip_threshold = 0.001
)

pip_cor_signal_trait1

# -------------------------------------------------------------------------
# Annotation summaries where available
# -------------------------------------------------------------------------

annotation_summaries <- lapply(
  names(fits),
  function(model_name) {
    x <- fits[[model_name]]
    
    if (!is.null(x$annotation_summary)) {
      out <- as.data.frame(x$annotation_summary)
      out$model <- model_name
      return(out)
    }
    
    NULL
  }
)

annotation_summaries <- Filter(Negate(is.null), annotation_summaries)
annotation_summaries

# -------------------------------------------------------------------------
# SBayesR / SBayesRC component summaries
# -------------------------------------------------------------------------

bayesr_component_summaries <- lapply(
  c("bayesR", "bayesR_MH", "sbayesrc", "sbayesrc_MH"),
  function(model_name) {
    x <- fits[[model_name]]
    out <- summarise_stblr_bayesr_components(x)
    out$model <- model_name
    out[, c("model", setdiff(names(out), "model")), drop = FALSE]
  }
)

bayesr_component_summaries <- do.call(rbind, bayesr_component_summaries)
bayesr_component_summaries

# -------------------------------------------------------------------------
# Compact input helper
# -------------------------------------------------------------------------

compact_input <- function(fit) {
  x <- fit$input
  x[c(
    "method",
    "model",
    "backend",
    "data_level",
    "annotation_model",
    "annotations",
    "updateLDswap",
    "ld_swap_prob",
    "ld_swap_r2",
    "ld_swap_moves",
    "nchains",
    "keep_chains",
    "n",
    "m",
    "nt",
    "nit",
    "nburn",
    "nthin",
    "seed"
  )]
}

lapply(fits, compact_input)


fit <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = A,
  annotation_model = "sbayesrc",
  nit = 1000,
  nburn = 100
)


fit_prior <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = list(
    A = A,
    fixed_pi_marker = sim$pi_marker,
    fixed_vb_multiplier = sim$vb_multiplier,
    use_pi_marker = TRUE,
    use_vb_multiplier = TRUE
  ),
  annotation_model = "prior",
  nit = 1000,
  nburn = 100
)

fit_learned <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = A,
  annotation_model = "learned",
  learn_pi_annot = TRUE,
  learn_vb_annot = TRUE,
  nit = 1000,
  nburn = 100
)

group <- ifelse(A[, 1] != 0, "annotated", "background")
names(group) <- rownames(A)

fit_group <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = group,
  annotation_model = "group",
  group_names = c("annotated", "background"),
  group_pi_init = c(0.002, 0.001),
  group_vb_multiplier_init = c(1.1, 1.0),
  nit = 1000,
  nburn = 100
)

fitC <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  nit = 1000,
  nburn = 100,
  seed = 10
)

fitC_MH <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  nit = 1000,
  nburn = 100,
  seed = 10,
  updateLDswap = TRUE,
  ld_swap_prob = 0.10,
  ld_swap_r2 = 0.05,
  ld_swap_moves = 5
)

fitR <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  nit = 1000,
  nburn = 100,
  seed = 10
)

fitR_MH <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  nit = 1000,
  nburn = 100,
  seed = 10,
  updateLDswap = TRUE,
  ld_swap_prob = 0.50,
  ld_swap_r2 = 0.001,
  ld_swap_moves = 20
)



## ------------------------------------------------------------
## Collect fits
## ------------------------------------------------------------

fits <- list(
  bayesC        = fitC,
  bayesC_MH     = fitC_MH,
  bayesR        = fitR,
  bayesR_MH     = fitR_MH,
  prior_annot   = fit_prior,
  learned_annot = fit_learned,
  group_annot   = fit_group,
  sbayesrc      = fit
)

## ------------------------------------------------------------
## Compact model metadata
## ------------------------------------------------------------

model_metadata <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    x <- fits[[model_name]]$input
    
    data.frame(
      model_name = model_name,
      method = ifelse(is.null(x$method), NA, x$method),
      model = ifelse(is.null(x$model), NA, x$model),
      backend = ifelse(is.null(x$backend), NA, x$backend),
      data_level = ifelse(is.null(x$data_level), NA, x$data_level),
      annotation_model = ifelse(is.null(x$annotation_model), NA, x$annotation_model),
      annotations = ifelse(is.null(x$annotations), FALSE, x$annotations),
      n_markers = nrow(fits[[model_name]]$dm),
      n_traits = ncol(fits[[model_name]]$dm),
      stringsAsFactors = FALSE
    )
  })
)


model_metadata


## ------------------------------------------------------------
## PIP and effect summaries by model and trait
## ------------------------------------------------------------

summarise_fit <- function(fit, model_name) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)
  
  trait_names <- colnames(dm)
  if (is.null(trait_names)) {
    trait_names <- paste0("T", seq_len(ncol(dm)))
  }
  
  data.frame(
    model = model_name,
    trait = trait_names,
    n_markers = nrow(dm),
    mean_pip = colMeans(dm, na.rm = TRUE),
    sum_pip = colSums(dm, na.rm = TRUE),
    max_pip = apply(dm, 2, max, na.rm = TRUE),
    n_pip_gt_0_001 = colSums(dm > 0.001, na.rm = TRUE),
    n_pip_gt_0_01 = colSums(dm > 0.01, na.rm = TRUE),
    n_pip_gt_0_05 = colSums(dm > 0.05, na.rm = TRUE),
    n_pip_gt_0_5 = colSums(dm > 0.5, na.rm = TRUE),
    n_pip_gt_0_95 = colSums(dm > 0.95, na.rm = TRUE),
    mean_abs_bm = colMeans(abs(bm), na.rm = TRUE),
    max_abs_bm = apply(abs(bm), 2, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

fit_summary <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    summarise_fit(fits[[model_name]], model_name)
  })
)

fit_summary


## ------------------------------------------------------------
## Sum PIP matrix: models x traits
## ------------------------------------------------------------

sum_pip_matrix <- do.call(
  rbind,
  lapply(fits, function(fit) colSums(as.matrix(fit$dm), na.rm = TRUE))
)

sum_pip_matrix

mean_pip_matrix <- do.call(
  rbind,
  lapply(fits, function(fit) colMeans(as.matrix(fit$dm), na.rm = TRUE))
)

mean_pip_matrix

## ------------------------------------------------------------
## Top markers by model and trait
## ------------------------------------------------------------

top_markers <- function(fit, model_name, trait = 1, top_n = 20) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)
  
  trait_index <- if (is.character(trait)) match(trait, colnames(dm)) else trait
  
  if (length(trait_index) != 1 || is.na(trait_index)) {
    stop("trait must be a valid column name or column index.")
  }
  
  marker <- rownames(dm)
  if (is.null(marker)) marker <- paste0("V", seq_len(nrow(dm)))
  
  trait_name <- colnames(dm)[trait_index]
  if (is.null(trait_name)) trait_name <- paste0("T", trait_index)
  
  out <- data.frame(
    model = model_name,
    trait = trait_name,
    marker = marker,
    pip = dm[, trait_index],
    bm = bm[, trait_index],
    abs_bm = abs(bm[, trait_index]),
    stringsAsFactors = FALSE
  )
  
  out <- out[order(-out$pip, -out$abs_bm), , drop = FALSE]
  utils::head(out, top_n)
}

top_trait1 <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    top_markers(fits[[model_name]], model_name, trait = 1, top_n = 20)
  })
)

top_trait1

all_top_markers <- do.call(
  rbind,
  lapply(seq_len(ncol(fitC$dm)), function(trait_index) {
    do.call(
      rbind,
      lapply(names(fits), function(model_name) {
        top_markers(
          fits[[model_name]],
          model_name,
          trait = trait_index,
          top_n = 20
        )
      })
    )
  })
)

all_top_markers


## ------------------------------------------------------------
## How many top markers are shared across models?
## ------------------------------------------------------------

top_marker_sets <- lapply(names(fits), function(model_name) {
  top_markers(fits[[model_name]], model_name, trait = 1, top_n = 20)$marker
})
names(top_marker_sets) <- names(fits)

shared_top_markers <- Reduce(intersect, top_marker_sets)
shared_top_markers

length(shared_top_markers)


pairwise_top_overlap <- expand.grid(
  model1 = names(fits),
  model2 = names(fits),
  stringsAsFactors = FALSE
)

pairwise_top_overlap$n_shared_top20 <- mapply(
  function(a, b) {
    length(intersect(top_marker_sets[[a]], top_marker_sets[[b]]))
  },
  pairwise_top_overlap$model1,
  pairwise_top_overlap$model2
)

pairwise_top_overlap


## ------------------------------------------------------------
## PIP correlation between models
## ------------------------------------------------------------

pip_correlation_by_trait <- function(fits, trait = 1, method = "spearman") {
  pip_mat <- do.call(
    cbind,
    lapply(fits, function(fit) {
      dm <- as.matrix(fit$dm)
      dm[, trait]
    })
  )
  
  colnames(pip_mat) <- names(fits)
  
  stats::cor(
    pip_mat,
    use = "pairwise.complete.obs",
    method = method
  )
}

pip_cor_trait1 <- pip_correlation_by_trait(fits, trait = 1)
pip_cor_trait1

pip_correlations <- lapply(seq_len(ncol(fitC$dm)), function(trait_index) {
  pip_correlation_by_trait(fits, trait = trait_index)
})

names(pip_correlations) <- colnames(fitC$dm)
pip_correlations


## ------------------------------------------------------------
## Annotation summaries where available
## ------------------------------------------------------------

annotation_summaries <- lapply(
  c("prior_annot", "learned_annot", "group_annot", "sbayesrc"),
  function(model_name) {
    x <- fits[[model_name]]
    
    if (!is.null(x$annotation_summary)) {
      out <- x$annotation_summary
      out$model <- model_name
      return(out)
    }
    
    NULL
  }
)

annotation_summaries <- Filter(Negate(is.null), annotation_summaries)
annotation_summaries


# Fixed marker-specific prior model ---------------------------------------


fit_prior_annot <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = list(
    A = A,
    fixed_pi_marker = sim$pi_marker,
    fixed_vb_multiplier = sim$vb_multiplier,
    use_pi_marker = TRUE,
    use_vb_multiplier = TRUE
  ),
  annotation_model = "prior",
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10
)

fit_learn_annot <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = A,
  annotation_model = "learned",
  learn_pi_annot = TRUE,
  learn_vb_annot = TRUE,
  rw_sd_eta_pi = 0.02,
  rw_sd_eta_vb = 0.02,
  annot_update_every = 10,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10
)

group <- ifelse(A[, 1] != 0, "annotated", "background")
names(group) <- rownames(A)

fit_group_annot <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = group,
  annotation_model = "group",
  group_names = c("annotated", "background"),
  group_pi_init = c(0.002, 0.001),
  group_vb_multiplier_init = c(1.1, 1.0),
  updatePi = TRUE,
  updateGroupVb = TRUE,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10
)

fit_sbayesrc <- stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  annotations = A,
  annotation_model = "sbayesrc",
  gamma = c(0, 0.01, 0.1, 1),
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10
)

# stblr_csr_prior_annot() uses the marker-specific inclusion probabilities and
# variance multipliers generated by mtsim_annotation().
fit_prior_annot <- stblr_csr_prior_annot(
 stats = stats,
 ld_prefix = ld_prefix,
 A = A,
 n = Glist$n,
 fixed_pi_marker = sim$pi_marker,
 fixed_vb_multiplier = sim$vb_multiplier,
 use_pi_marker = TRUE,
 use_vb_multiplier = TRUE,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# Learned annotation model ------------------------------------------------

# stblr_csr_learn_annot() learns annotation effects on marker inclusion
# probabilities and, optionally, marker-effect variances. These proposal scales
# are conservative demonstration settings. Real analyses require longer chains,
# convergence checks, and posterior diagnostics.
fit_learn_annot <- stblr_csr_learn_annot(
 stats = stats,
 ld_prefix = ld_prefix,
 A = A,
 n = Glist$n,
 learn_pi_annot = TRUE,
 learn_vb_annot = TRUE,
 rw_sd_eta_pi = 0.02,
 rw_sd_eta_vb = 0.02,
 annot_update_every = 10,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# Group annotation model --------------------------------------------------

# group is a length-m marker group vector in the same order as stats and the
# sparse LD. These are conservative demonstration settings. Real analyses
# require longer chains, convergence checks, and posterior diagnostics.
group <- ifelse(A[, 1] != 0, "annotated", "background")
names(group) <- rownames(A)

fit_group_annot <- stblr_csr_group_annot(
 stats = stats,
 ld_prefix = ld_prefix,
 group = group,
 n = Glist$n,
 group_names = c("annotated", "background"),
 group_pi_init = c(0.002, 0.001),
 group_vb_multiplier_init = c(1.1, 1.0),
 updatePi = TRUE,
 updateGroupVb = TRUE,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# SBayesRC-style annotation model -----------------------------------------

# stblr_csr_sbayesrc_generic() uses continuous overlapping annotations to
# model SBayesRC-style mixture-component probabilities. These are conservative
# demonstration settings. Real analyses require longer chains, convergence
# checks, and posterior diagnostics.
fit_sbayesrc <- stblr_csr_sbayesrc_generic(
 stats = stats,
 ld_prefix = ld_prefix,
 A = A,
 n = Glist$n,
 gamma = c(0, 0.01, 0.1, 1),
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# Inspect posterior summaries ---------------------------------------------

# Printed summaries are deliberately compact because marker-level posterior
# matrices can contain tens of thousands of rows. Full matrices remain
# available in fit$bm, fit$dm, fit$b, and other fit components.
compact_fit_summary <- function(fit, model_name) {
 bm <- as.matrix(fit$bm)
 dm <- as.matrix(fit$dm)
 trait_names <- colnames(dm)
 if (is.null(trait_names)) trait_names <- paste0("T", seq_len(ncol(dm)))

 mean_trace <- function(name) {
  x <- fit[[name]]
  if (is.null(x)) return(rep(NA_real_, ncol(dm)))
  x <- as.matrix(x)
  if (ncol(x) != ncol(dm)) return(rep(NA_real_, ncol(dm)))
  colMeans(x, na.rm = TRUE)
 }

 covariance_summary <- do.call(
  rbind,
  lapply(c("covb", "covg", "cove", "rb", "rg", "re"), function(name) {
   x <- fit[[name]]
   if (is.null(x)) return(NULL)
   x <- as.matrix(x)
   if (nrow(x) != ncol(x) || any(!is.finite(x))) return(NULL)
   off_diagonal <- x[row(x) != col(x)]
   data.frame(
    matrix = name,
    mean_diagonal = mean(diag(x)),
    mean_abs_off_diagonal = if (length(off_diagonal)) {
     mean(abs(off_diagonal))
    } else {
     NA_real_
    },
    stringsAsFactors = FALSE
   )
  })
 )

 list(
  trait_summary = data.frame(
   model = model_name,
   trait = trait_names,
   n_markers = nrow(dm),
   n_traits = ncol(dm),
   pip_sum = colSums(dm, na.rm = TRUE),
   mean_abs_bm = colMeans(abs(bm), na.rm = TRUE),
   mean_vbs = mean_trace("vbs"),
   mean_vgs = mean_trace("vgs"),
   mean_ves = mean_trace("ves"),
   row.names = NULL,
   stringsAsFactors = FALSE
  ),
  covariance_summary = covariance_summary
 )
}

top_marker_summary <- function(fit, trait = 1, top_n = 10) {
 bm <- as.matrix(fit$bm)
 dm <- as.matrix(fit$dm)
 trait_index <- if (is.character(trait)) match(trait, colnames(dm)) else trait
 if (length(trait_index) != 1 || is.na(trait_index) ||
     trait_index < 1 || trait_index > ncol(dm)) {
  stop("trait must identify one column of fit$dm.")
 }
 marker <- rownames(dm)
 if (is.null(marker)) marker <- paste0("V", seq_len(nrow(dm)))
 trait_name <- colnames(dm)[trait_index]
 if (is.null(trait_name)) trait_name <- paste0("T", trait_index)

 out <- data.frame(
  marker = marker,
  trait = trait_name,
  pip = dm[, trait_index],
  abs_bm = abs(bm[, trait_index]),
  bm = bm[, trait_index],
  stringsAsFactors = FALSE
 )
 out <- out[order(-out$pip, -out$abs_bm), , drop = FALSE]
 utils::head(out, top_n)
}

annotation_model_summary <- function(sim, fit, model_name) {
 if (!exists("summarize_annotation_signal", mode = "function")) {
  stop("summarize_annotation_signal() is required for annotation summaries.")
 }
 out <- summarize_annotation_signal(sim, fit)
 keep <- c(
  "annotation", "size", "n_causal", "causal_rate_in_set",
  grep("^(mean_dm_|mean_abs_bm_)", names(out), value = TRUE)
 )
 out <- out[, intersect(keep, names(out)), drop = FALSE]
 out$model <- model_name
 out[, c("model", setdiff(names(out), "model")), drop = FALSE]
}

fits <- list(
 fixed_prior = fit_prior_annot,
 learned_annotation = fit_learn_annot,
 group_annotation = fit_group_annot,
 sbayesrc = fit_sbayesrc
)

for (model_name in names(fits)) {
 cat("\nCompact fit summary:", model_name, "\n")
 print(compact_fit_summary(fits[[model_name]], model_name))
 cat("\nTop 10 markers for trait 1:", model_name, "\n")
 print(top_marker_summary(fits[[model_name]], trait = 1, top_n = 10))
}

for (model_name in c("fixed_prior", "learned_annotation", "group_annotation")) {
 cat("\nAnnotation signal summary:", model_name, "\n")
 print(annotation_model_summary(sim, fits[[model_name]], model_name))
}

# These exported helpers provide compact SBayesRC prior diagnostics.
if (!is.null(fit_sbayesrc$alpha)) {
 sbayesrc_A <- fit_sbayesrc$input$A
 gamma <- fit_sbayesrc$input$gamma
 trait_names <- names(fit_sbayesrc$alpha)
 if (is.null(trait_names)) {
  trait_names <- paste0("T", seq_along(fit_sbayesrc$alpha))
 }

 annotation_gamma <- data.frame(
  annotation = rownames(fit_sbayesrc$alpha[[1]]),
  stringsAsFactors = FALSE
 )
 marker_gamma_summary <- list()
 marker_gamma_pip_correlation <- numeric(length(fit_sbayesrc$alpha))

 for (t in seq_along(fit_sbayesrc$alpha)) {
  marker_gamma <- sbayesrc_marker_gamma_mean(
   sbayesrc_A, fit_sbayesrc$alpha[[t]], gamma
  )
  annotation_gamma[[trait_names[t]]] <- sbayesrc_annotation_gamma_mean(
   fit_sbayesrc$alpha[[t]], gamma
  )
  marker_gamma_summary[[trait_names[t]]] <- summary(marker_gamma)
  marker_gamma_pip_correlation[t] <- stats::cor(
   marker_gamma, fit_sbayesrc$dm[, t], use = "complete.obs"
  )
 }
 names(marker_gamma_pip_correlation) <- trait_names

 cat("\nSBayesRC annotation-level expected gamma:\n")
 print(annotation_gamma)
 cat("\nSBayesRC marker expected gamma summaries:\n")
 print(marker_gamma_summary)
 cat("\nSBayesRC marker expected gamma/PIP correlations:\n")
 print(marker_gamma_pip_correlation)
}

# Real analyses require longer chains, convergence checks, and posterior
# predictive diagnostics before interpreting these compact summaries.
