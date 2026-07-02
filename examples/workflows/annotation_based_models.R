# Annotation-based ST-BLR workflow
#
# Compares annotation-unaware CSR BayesC/BayesR models with annotation-aware
# CSR models: fixed-prior, learned annotation, group annotation, and SBayesRC.
#
# Includes LD-swap/MH variants for comparison with annotation-unaware models.
#
# The MCMC settings shown here are demonstration settings. Real analyses need
# longer chains and appropriate convergence and posterior predictive checks.

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

# mtsim_annotation() simulates overlapping annotations, enriches causal-marker
# sampling in selected annotations, and creates marker-specific prior values.
# Replace this simulation with phenotype and biological annotation data from
# the study if desired.

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

# Inspect the simulated annotation enrichment before model fitting.

summarize_annotation_signal(sim)

# Compute summary statistics -----------------------------------------------

stats <- make_stats(
  Glist = Glist,
  y = y,
  nthreads = nthreads
)

# Compute sparse LD --------------------------------------------------------

Glist <- make_sparseLD(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_test"),
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

# Prepare annotation matrix ------------------------------------------------

# A is an m x K marker annotation matrix. Rows of A must correspond exactly
# to the marker ordering in stats and the sparse-LD files. Real analyses should
# replace sim$annot with biological annotations in the same marker order.

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

# LD-swap/MH settings ------------------------------------------------------

# Conservative BayesC-like comparison.
mh_conservative <- list(
  updateLDswap = TRUE,
  ld_swap_prob = 0.10,
  ld_swap_r2 = 0.05,
  ld_swap_moves = 5
)

# Permissive LD-relocation comparison.
mh_permissive <- list(
  updateLDswap = TRUE,
  ld_swap_prob = 0.50,
  ld_swap_r2 = 0.001,
  ld_swap_moves = 20
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
    list(method = "bayesC"),
    mh_conservative
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
    list(method = "bayesR"),
    mh_permissive
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
# Annotation-aware CSR models with conservative LD-swap/MH
# -------------------------------------------------------------------------

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
      annotation_model = "prior"
    ),
    mh_conservative
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
      updateGroupVb = TRUE
    ),
    mh_conservative
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
      annot_update_every = 10
    ),
    mh_conservative
  )
)

fit_sbayesrc_MH <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(
      annotations = A,
      annotation_model = "sbayesrc",
      gamma = c(0, 0.01, 0.1, 1)
    ),
    mh_permissive
  )
)

# -------------------------------------------------------------------------
# Annotation-aware CSR models with permissive LD-swap/MH
# -------------------------------------------------------------------------
# These settings are useful for checking whether BayesC-like LD-swap moves
# are actively accepted. They are not necessarily the preferred default for
# final analyses.

fit_prior_MH2 <- do.call(
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
    ),
    mh_permissive
  )
)

fit_group_MH2 <- do.call(
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
    ),
    mh_permissive
  )
)

fit_learned_MH2 <- do.call(
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
    ),
    mh_permissive
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
  bayesC = fitC,
  bayesC_MH = fitC_MH,
  bayesR = fitR,
  bayesR_MH = fitR_MH,
  prior = fit_prior,
  prior_MH = fit_prior_MH,
  prior_MH2 = fit_prior_MH2,
  learned = fit_learned,
  learned_MH = fit_learned_MH,
  learned_MH2 = fit_learned_MH2,
  group = fit_group,
  group_MH = fit_group_MH,
  group_MH2 = fit_group_MH2,
  sbayesrc = fit_sbayesrc,
  sbayesrc_MH = fit_sbayesrc_MH
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
  
  trait_index <- if (is.character(trait)) {
    match(trait, colnames(dm))
  } else {
    trait
  }
  
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

pip_correlation_signal_markers <- function(
    fits,
    trait = 1,
    pip_threshold = 0.001,
    method = "spearman"
) {
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
  fits = fits,
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
# BayesR / SBayesRC component summaries
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

# -------------------------------------------------------------------------
# Focused comparison of annotation-aware BayesC-like models
# -------------------------------------------------------------------------

rbind(
  prior = colSums(fit_prior$dm),
  prior_MH = colSums(fit_prior_MH$dm),
  prior_MH2 = colSums(fit_prior_MH2$dm),
  group = colSums(fit_group$dm),
  group_MH = colSums(fit_group_MH$dm),
  group_MH2 = colSums(fit_group_MH2$dm),
  learned = colSums(fit_learned$dm),
  learned_MH = colSums(fit_learned_MH$dm),
  learned_MH2 = colSums(fit_learned_MH2$dm)
)

rbind(
  prior_MH2 = fit_prior_MH2$ld_swap,
  group_MH2 = fit_group_MH2$ld_swap,
  learned_MH2 = fit_learned_MH2$ld_swap
)



# -------------------------------------------------------------------------
# True causal marker sets
# -------------------------------------------------------------------------

causal_shared <- sim$causal$shared
causal_specific <- sim$causal$specific
causal_all <- sim$causal$all

# Trait-specific true causal markers:
# shared markers + trait-specific markers

causal_by_trait <- lapply(names(causal_specific), function(trait) {
  unique(c(causal_shared, causal_specific[[trait]]))
})

names(causal_by_trait) <- names(causal_specific)

causal_by_trait


# -------------------------------------------------------------------------
# Detection power by PIP threshold
# -------------------------------------------------------------------------

causal_detection_summary <- function(fit, model_name, causal_by_trait,
                                     thresholds = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.95)) {
  dm <- as.matrix(fit$dm)
  marker <- rownames(dm)
  
  if (is.null(marker)) {
    stop("fit$dm must have marker rownames.")
  }
  
  trait_names <- colnames(dm)
  if (is.null(trait_names)) {
    trait_names <- names(causal_by_trait)
  }
  
  out <- list()
  
  for (trait in trait_names) {
    causal <- causal_by_trait[[trait]]
    causal <- intersect(causal, marker)
    
    if (length(causal) == 0) {
      next
    }
    
    pip <- dm[, trait]
    names(pip) <- marker
    
    trait_out <- lapply(thresholds, function(thr) {
      detected <- causal[pip[causal] >= thr]
      
      data.frame(
        model = model_name,
        trait = trait,
        threshold = thr,
        n_causal = length(causal),
        n_detected = length(detected),
        power = length(detected) / length(causal),
        stringsAsFactors = FALSE
      )
    })
    
    out[[trait]] <- do.call(rbind, trait_out)
  }
  
  do.call(rbind, out)
}

power_by_threshold <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    causal_detection_summary(
      fit = fits[[model_name]],
      model_name = model_name,
      causal_by_trait = causal_by_trait
    )
  })
)

power_by_threshold


power_matrix <- function(power_by_threshold, threshold = 0.5) {
  x <- subset(power_by_threshold, threshold == !!threshold)
  
  mat <- xtabs(power ~ model + trait, data = x)
  as.matrix(mat)
}

power_matrix(power_by_threshold, threshold = 0.05)
power_matrix(power_by_threshold, threshold = 0.5)
power_matrix(power_by_threshold, threshold = 0.95)

power_matrix <- function(power_by_threshold, threshold_value = 0.5) {
  x <- power_by_threshold[power_by_threshold$threshold == threshold_value, ]
  
  mat <- xtabs(power ~ model + trait, data = x)
  as.matrix(mat)
}

power_matrix(power_by_threshold, threshold_value = 0.05)
power_matrix(power_by_threshold, threshold_value = 0.5)
power_matrix(power_by_threshold, threshold_value = 0.95)


# -------------------------------------------------------------------------
# Detection power among top-N markers
# -------------------------------------------------------------------------

causal_topn_summary <- function(fit, model_name, causal_by_trait,
                                top_n = c(20, 50, 100, 200, 500, 1000)) {
  dm <- as.matrix(fit$dm)
  marker <- rownames(dm)
  
  if (is.null(marker)) {
    stop("fit$dm must have marker rownames.")
  }
  
  trait_names <- colnames(dm)
  if (is.null(trait_names)) {
    trait_names <- names(causal_by_trait)
  }
  
  out <- list()
  
  for (trait in trait_names) {
    causal <- causal_by_trait[[trait]]
    causal <- intersect(causal, marker)
    
    if (length(causal) == 0) {
      next
    }
    
    pip <- dm[, trait]
    names(pip) <- marker
    
    ranked_marker <- names(sort(pip, decreasing = TRUE))
    
    trait_out <- lapply(top_n, function(n) {
      top_marker <- ranked_marker[seq_len(min(n, length(ranked_marker)))]
      detected <- intersect(causal, top_marker)
      
      data.frame(
        model = model_name,
        trait = trait,
        top_n = n,
        n_causal = length(causal),
        n_detected = length(detected),
        power = length(detected) / length(causal),
        precision = length(detected) / length(top_marker),
        stringsAsFactors = FALSE
      )
    })
    
    out[[trait]] <- do.call(rbind, trait_out)
  }
  
  do.call(rbind, out)
}

topn_power <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    causal_topn_summary(
      fit = fits[[model_name]],
      model_name = model_name,
      causal_by_trait = causal_by_trait
    )
  })
)

topn_power


topn_power_matrix <- function(topn_power, top_n_value = 100) {
  x <- topn_power[topn_power$top_n == top_n_value, ]
  
  mat <- xtabs(power ~ model + trait, data = x)
  as.matrix(mat)
}

topn_precision_matrix <- function(topn_power, top_n_value = 100) {
  x <- topn_power[topn_power$top_n == top_n_value, ]
  
  mat <- xtabs(precision ~ model + trait, data = x)
  as.matrix(mat)
}

topn_power_matrix(topn_power, top_n_value = 20)
topn_power_matrix(topn_power, top_n_value = 100)
topn_power_matrix(topn_power, top_n_value = 500)

topn_precision_matrix(topn_power, top_n_value = 20)
topn_precision_matrix(topn_power, top_n_value = 100)

# -------------------------------------------------------------------------
# Rank and PIP of true causal markers
# -------------------------------------------------------------------------

causal_marker_ranks <- function(fit, model_name, causal_by_trait) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)
  marker <- rownames(dm)
  
  if (is.null(marker)) {
    stop("fit$dm must have marker rownames.")
  }
  
  trait_names <- colnames(dm)
  if (is.null(trait_names)) {
    trait_names <- names(causal_by_trait)
  }
  
  out <- list()
  
  for (trait in trait_names) {
    causal <- causal_by_trait[[trait]]
    causal <- intersect(causal, marker)
    
    pip <- dm[, trait]
    names(pip) <- marker
    
    effect <- bm[, trait]
    names(effect) <- marker
    
    rank <- rank(-pip, ties.method = "min")
    names(rank) <- marker
    
    out[[trait]] <- data.frame(
      model = model_name,
      trait = trait,
      marker = causal,
      pip = pip[causal],
      rank = rank[causal],
      bm = effect[causal],
      abs_bm = abs(effect[causal]),
      causal_type = ifelse(causal %in% causal_shared, "shared", "trait_specific"),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

causal_ranks <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    causal_marker_ranks(
      fit = fits[[model_name]],
      model_name = model_name,
      causal_by_trait = causal_by_trait
    )
  })
)

causal_ranks


rank_summary <- aggregate(
  cbind(pip, rank) ~ model + trait + causal_type,
  data = causal_ranks,
  FUN = function(x) c(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
)

rank_summary


rank_summary_simple <- do.call(
  rbind,
  lapply(split(causal_ranks, list(causal_ranks$model, causal_ranks$trait, causal_ranks$causal_type)), function(x) {
    if (nrow(x) == 0) return(NULL)
    
    data.frame(
      model = x$model[1],
      trait = x$trait[1],
      causal_type = x$causal_type[1],
      mean_pip = mean(x$pip, na.rm = TRUE),
      median_pip = median(x$pip, na.rm = TRUE),
      mean_rank = mean(x$rank, na.rm = TRUE),
      median_rank = median(x$rank, na.rm = TRUE),
      n_rank_top20 = sum(x$rank <= 20, na.rm = TRUE),
      n_rank_top100 = sum(x$rank <= 100, na.rm = TRUE),
      n_rank_top500 = sum(x$rank <= 500, na.rm = TRUE),
      n_causal = nrow(x),
      stringsAsFactors = FALSE
    )
  })
)

rank_summary_simple


# -------------------------------------------------------------------------
# Shared vs trait-specific causal detection
# -------------------------------------------------------------------------

causal_type_power <- function(causal_ranks, pip_threshold = 0.5) {
  x <- causal_ranks
  x$detected <- x$pip >= pip_threshold
  
  out <- aggregate(
    detected ~ model + trait + causal_type,
    data = x,
    FUN = mean
  )
  
  names(out)[names(out) == "detected"] <- "power"
  out$pip_threshold <- pip_threshold
  
  out
}

causal_type_power_005 <- causal_type_power(causal_ranks, pip_threshold = 0.05)
causal_type_power_05 <- causal_type_power(causal_ranks, pip_threshold = 0.5)
causal_type_power_095 <- causal_type_power(causal_ranks, pip_threshold = 0.95)

causal_type_power_005
causal_type_power_05
causal_type_power_095


power_matrix(power_by_threshold, threshold_value = 0.5)
topn_power_matrix(topn_power, top_n_value = 100)
topn_power_matrix(topn_power, top_n_value = 500)
rank_summary_simple


topn_power_matrix(topn_power, top_n_value = 100)

mean_rank_matrix <- function(rank_summary_simple, causal_type_value = "shared") {
  x <- rank_summary_simple[rank_summary_simple$causal_type == causal_type_value, ]
  xtabs(mean_rank ~ model + trait, data = x)
}

mean_rank_matrix(rank_summary_simple, "shared")
mean_rank_matrix(rank_summary_simple, "trait_specific")


trait_specific_rank_summary <- rank_summary_simple[
  rank_summary_simple$causal_type == "trait_specific",
]

trait_specific_rank_summary[
  order(trait_specific_rank_summary$trait,
        trait_specific_rank_summary$mean_rank),
  c(
    "model",
    "trait",
    "mean_pip",
    "median_pip",
    "mean_rank",
    "median_rank",
    "n_rank_top20",
    "n_rank_top100",
    "n_rank_top500",
    "n_causal"
  )
]

trait_specific_rank_summary <- rank_summary_simple[
  rank_summary_simple$causal_type == "trait_specific",
]

trait_specific_rank_summary$score <- with(
  trait_specific_rank_summary,
  n_rank_top100 / n_causal
)

trait_specific_rank_summary[
  order(-trait_specific_rank_summary$score,
        trait_specific_rank_summary$mean_rank),
]

trait_specific_top100 <- causal_ranks[
  causal_ranks$causal_type == "trait_specific",
]

trait_specific_top100$detected_top100 <- trait_specific_top100$rank <= 100

aggregate(
  detected_top100 ~ model,
  data = trait_specific_top100,
  FUN = mean
)



