# Annotation-based ST-BLR workflow
#
# Compares annotation-unaware CSR BayesC/BayesR models with annotation-aware
# CSR fixed-prior, learned-annotation, group-prior, and SBayesRC models.
#
# See docs/dev/stblr_backend_computation_inventory.md for a detailed overview
# of backend fields and return-value conventions.
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

# mtsim_annotation() simulates overlapping annotations, enriched causal-marker
# sampling, and marker-specific prior values. Replace this block with study
# phenotypes and biological annotations for real data analyses.

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

summarize_annotation_signal(sim)

# Compute summary statistics -----------------------------------------------

stats <- make_stats(
  Glist = Glist,
  y = y,
  chr = chr,
  nthreads = nthreads
)

# Compute sparse LD --------------------------------------------------------

Glist <- make_sparseLD(
  Glist = Glist,
  out_prefix = file.path(data_dir, "ld_test"),
  chr = chr,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = nthreads
)

# Prepare annotation matrix ------------------------------------------------

# A is an m x K marker annotation matrix. Rows of A must match the marker
# order in stats, sparse LD, and returned marker-level posterior summaries.

marker_id <- Glist$rsidsLD[[chr]]
m <- length(marker_id)

A <- sim$annot
rownames(A) <- marker_id

stopifnot(
  nrow(A) == m,
  identical(rownames(A), marker_id)
)

# Prepare group annotation -------------------------------------------------

group <- ifelse(A[, 1] != 0, "annotated", "background")
names(group) <- rownames(A)

# Shared model settings ----------------------------------------------------

base_args <- list(
  stats = stats,
  Glist = Glist,
  nit = nit,
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
  annotation_model = "learned",
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
  annotation_model = "sbayesrc",
  gamma = c(0, 0.01, 0.1, 1)
)

# LD-swap/MH settings ------------------------------------------------------

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

# Annotation-unaware CSR models -------------------------------------------

fitC <- do.call(
  stblr_csr,
  c(base_args, list(method = "bayesC"))
)

fitC_MH <- do.call(
  stblr_csr,
  c(base_args, list(method = "bayesC"), mh_conservative)
)

fitR <- do.call(
  stblr_csr,
  c(base_args, list(method = "bayesR"))
)

fitR_MH <- do.call(
  stblr_csr,
  c(base_args, list(method = "bayesR"), mh_permissive)
)

# Annotation-aware CSR models ---------------------------------------------

fit_prior <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(annotations = prior_annotations, annotation_model = "prior")
  )
)

fit_learned <- do.call(
  stblr_csr_annot,
  c(base_args, learned_args)
)

fit_group <- do.call(
  stblr_csr_annot,
  c(base_args, group_args)
)

fit_sbayesrc <- do.call(
  stblr_csr_annot,
  c(base_args, sbayesrc_args)
)

# Annotation-aware CSR models with conservative LD-swap/MH -----------------

fit_prior_MH <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(annotations = prior_annotations, annotation_model = "prior"),
    mh_conservative
  )
)

fit_learned_MH <- do.call(
  stblr_csr_annot,
  c(base_args, learned_args, mh_conservative)
)

fit_group_MH <- do.call(
  stblr_csr_annot,
  c(base_args, group_args, mh_conservative)
)

# Annotation-aware CSR models with permissive LD-swap/MH -------------------

# These settings are useful diagnostics for active LD-swap acceptance. They are
# not necessarily preferred defaults for final BayesC-like analyses.

fit_prior_MH2 <- do.call(
  stblr_csr_annot,
  c(
    base_args,
    list(annotations = prior_annotations, annotation_model = "prior"),
    mh_permissive
  )
)

fit_learned_MH2 <- do.call(
  stblr_csr_annot,
  c(base_args, learned_args, mh_permissive)
)

fit_group_MH2 <- do.call(
  stblr_csr_annot,
  c(base_args, group_args, mh_permissive)
)

fit_sbayesrc_MH <- do.call(
  stblr_csr_annot,
  c(base_args, sbayesrc_args, mh_permissive)
)

# Optional multi-chain examples -------------------------------------------

# Short demonstration runs only. Increase nit/nburn and inspect convergence
# diagnostics before using multi-chain output for substantive analysis.

fit_sbayesrc_2c <- do.call(
  stblr_csr_annot,
  c(
    modifyList(
      base_args,
      list(
        nit = 200,
        nburn = 50
      )
    ),
    sbayesrc_args,
    list(
      nchains = 2,
      chain_seeds = c(10, 20),
      keep_chains = TRUE
    )
  )
)

fit_sbayesrc_2c$input[c("backend", "nchains", "keep_chains")]
length(fit_sbayesrc_2c$chains)

fit_prior_MH_2c <- do.call(
  stblr_csr_annot,
  c(
    modifyList(
      base_args,
      list(
        nit = 200,
        nburn = 50
      )
    ),
    list(annotations = prior_annotations, annotation_model = "prior"),
    mh_conservative,
    list(
      nchains = 2,
      chain_seeds = c(10, 20),
      keep_chains = TRUE
    )
  )
)

fit_prior_MH_2c$input[c("backend", "nchains", "keep_chains", "updateLDswap")]
fit_prior_MH_2c$ld_swap

# Collect fits -------------------------------------------------------------

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

# Compact metadata ---------------------------------------------------------

metadata_value <- function(x, name, default = NA) {
  value <- x[[name]]
  if (is.null(value)) default else value
}

model_metadata <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    x <- fits[[model_name]]$input

    data.frame(
      model_name = model_name,
      method = metadata_value(x, "method"),
      model = metadata_value(x, "model"),
      backend = metadata_value(x, "backend"),
      data_level = metadata_value(x, "data_level"),
      annotations = metadata_value(x, "annotations", FALSE),
      annotation_model = metadata_value(x, "annotation_model"),
      updateLDswap = metadata_value(x, "updateLDswap", FALSE),
      n_markers = nrow(fits[[model_name]]$dm),
      n_traits = ncol(fits[[model_name]]$dm),
      stringsAsFactors = FALSE
    )
  })
)

model_metadata

field_inventory <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]

    data.frame(
      model = model_name,
      backend = fit$input$backend,
      data_level = fit$input$data_level,
      annotations = fit$input$annotations,
      annotation_model = ifelse(
        is.null(fit$input$annotation_model),
        NA,
        fit$input$annotation_model
      ),
      has_dm = !is.null(fit$dm),
      has_bm = !is.null(fit$bm),
      has_vle = !is.null(fit$vle),
      has_vld = !is.null(fit$vld),
      has_comp_prob = !is.null(fit$comp_prob),
      has_dm_component_mean = !is.null(fit$dm_component_mean),
      has_ld_swap = !is.null(fit$ld_swap),
      stringsAsFactors = FALSE
    )
  })
)

field_inventory

# Compact fit summaries ----------------------------------------------------

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

# LD-swap diagnostics ------------------------------------------------------

ld_swap_summary <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    x <- fits[[model_name]]$ld_swap
    if (is.null(x)) {
      return(NULL)
    }

    x <- as.data.frame(x)
    x$model <- model_name
    x$trait <- rownames(x)
    x[, c("model", "trait", setdiff(names(x), c("model", "trait")))]
  })
)

ld_swap_summary

# True causal marker recovery ----------------------------------------------

causal_shared <- sim$causal$shared
causal_specific <- sim$causal$specific
causal_all <- sim$causal$all

causal_by_trait <- lapply(names(causal_specific), function(trait) {
  unique(c(causal_shared, causal_specific[[trait]]))
})
names(causal_by_trait) <- names(causal_specific)

causal_detection_summary <- function(
    fit,
    model_name,
    causal_by_trait,
    thresholds = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.95)
) {
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
    causal <- intersect(causal_by_trait[[trait]], marker)
    if (length(causal) == 0L) {
      next
    }

    pip <- dm[, trait]
    names(pip) <- marker

    out[[trait]] <- do.call(
      rbind,
      lapply(thresholds, function(thr) {
        n_detected <- sum(pip[causal] >= thr, na.rm = TRUE)
        data.frame(
          model = model_name,
          trait = trait,
          threshold = thr,
          n_causal = length(causal),
          n_detected = n_detected,
          power = n_detected / length(causal),
          stringsAsFactors = FALSE
        )
      })
    )
  }

  do.call(rbind, out)
}

causal_topn_summary <- function(
    fit,
    model_name,
    causal_by_trait,
    top_n = c(20, 50, 100, 200, 500, 1000)
) {
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
    causal <- intersect(causal_by_trait[[trait]], marker)
    if (length(causal) == 0L) {
      next
    }

    pip <- dm[, trait]
    names(pip) <- marker
    ranked_marker <- names(sort(pip, decreasing = TRUE))

    out[[trait]] <- do.call(
      rbind,
      lapply(top_n, function(n) {
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
      })
    )
  }

  do.call(rbind, out)
}

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
    causal <- intersect(causal_by_trait[[trait]], marker)
    pip <- dm[, trait]
    effect <- bm[, trait]

    names(pip) <- marker
    names(effect) <- marker
    marker_rank <- rank(-pip, ties.method = "min")
    names(marker_rank) <- marker

    out[[trait]] <- data.frame(
      model = model_name,
      trait = trait,
      marker = causal,
      pip = pip[causal],
      rank = marker_rank[causal],
      bm = effect[causal],
      abs_bm = abs(effect[causal]),
      causal_type = ifelse(
        causal %in% causal_shared,
        "shared",
        "trait_specific"
      ),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

topn_power_matrix <- function(topn_power, top_n_value = 100) {
  x <- topn_power[topn_power$top_n == top_n_value, ]
  as.matrix(xtabs(power ~ model + trait, data = x))
}

topn_precision_matrix <- function(topn_power, top_n_value = 100) {
  x <- topn_power[topn_power$top_n == top_n_value, ]
  as.matrix(xtabs(precision ~ model + trait, data = x))
}

mean_rank_matrix <- function(
    rank_summary_simple,
    causal_type_value = "trait_specific"
) {
  x <- rank_summary_simple[
    rank_summary_simple$causal_type == causal_type_value,
  ]
  xtabs(mean_rank ~ model + trait, data = x)
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

rank_summary_simple <- do.call(
  rbind,
  lapply(
    split(causal_ranks, list(
      causal_ranks$model,
      causal_ranks$trait,
      causal_ranks$causal_type
    )),
    function(x) {
      if (nrow(x) == 0L) {
        return(NULL)
      }

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
    }
  )
)

power_by_threshold
topn_power
causal_ranks
rank_summary_simple

topn_power_matrix(topn_power, top_n_value = 20)
topn_power_matrix(topn_power, top_n_value = 50)
topn_power_matrix(topn_power, top_n_value = 100)

topn_precision_matrix(topn_power, top_n_value = 20)
topn_precision_matrix(topn_power, top_n_value = 100)

mean_rank_matrix(rank_summary_simple, "shared")
mean_rank_matrix(rank_summary_simple, "trait_specific")

# Top-marker summaries -----------------------------------------------------

top_markers <- function(fit, model_name, trait = 1, top_n = 20) {
  dm <- as.matrix(fit$dm)
  bm <- as.matrix(fit$bm)

  trait_index <- if (is.character(trait)) {
    match(trait, colnames(dm))
  } else {
    trait
  }

  if (length(trait_index) != 1L || is.na(trait_index)) {
    stop("trait must be a valid column name or column index.")
  }

  marker <- rownames(dm)
  if (is.null(marker)) {
    marker <- paste0("V", seq_len(nrow(dm)))
  }

  trait_name <- colnames(dm)[trait_index]
  if (is.null(trait_name)) {
    trait_name <- paste0("T", trait_index)
  }

  out <- data.frame(
    model = model_name,
    trait = trait_name,
    marker = marker,
    pip = dm[, trait_index],
    bm = bm[, trait_index],
    abs_bm = abs(bm[, trait_index]),
    true_causal = marker %in% causal_all,
    stringsAsFactors = FALSE
  )

  out <- out[order(-out$pip, -out$abs_bm), ]
  utils::head(out, top_n)
}

top_trait1 <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    top_markers(fits[[model_name]], model_name, trait = 1, top_n = 100)
  })
)

head(top_trait1, 40)

# write.csv(top_trait1, file.path(data_dir, "stblr_annotation_top_trait1.csv"),
#           row.names = FALSE)

# Top-marker overlap -------------------------------------------------------

top_overlap_matrix <- function(fits, trait = 1, top_n = 100) {
  top_sets <- lapply(fits, function(fit) {
    dm <- as.matrix(fit$dm)
    marker <- rownames(dm)
    if (is.null(marker)) {
      marker <- paste0("V", seq_len(nrow(dm)))
    }

    ord <- order(-dm[, trait])
    marker[ord[seq_len(min(top_n, length(ord)))]]
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

# PIP correlations ---------------------------------------------------------

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

# Annotation summaries -----------------------------------------------------

annotation_summaries <- lapply(names(fits), function(model_name) {
  x <- fits[[model_name]]
  if (is.null(x$annotation_summary)) {
    return(NULL)
  }

  out <- as.data.frame(x$annotation_summary)
  out$model <- model_name
  out[, c("model", setdiff(names(out), "model"))]
})

annotation_summaries <- Filter(Negate(is.null), annotation_summaries)
annotation_summaries

# BayesR / SBayesRC component summaries -----------------------------------

bayesr_component_summaries <- lapply(
  c("bayesR", "bayesR_MH", "sbayesrc", "sbayesrc_MH"),
  function(model_name) {
    x <- fits[[model_name]]
    out <- summarise_stblr_bayesr_components(x)
    out$model <- model_name
    out[, c("model", setdiff(names(out), "model"))]
  }
)

bayesr_component_summaries

# Posterior component summaries -------------------------------------------

posterior_summaries <- lapply(fits, function(fit) {
  summarise_stblr_posterior(fit)
})

lapply(posterior_summaries, head)

lapply(fits, function(fit) {
  fit$input[c(
    "method",
    "model",
    "backend",
    "data_level",
    "annotations",
    "annotation_model",
    "updateLDswap",
    "nchains",
    "keep_chains"
  )]
})

# Example trace plot call for interactive use:
# plot_stblr_posterior(posterior_summaries$prior_MH)

# Compact input helper -----------------------------------------------------

compact_input <- function(fit) {
  x <- fit$input

  x[c(
    "method",
    "model",
    "backend",
    "data_level",
    "annotations",
    "annotation_model",
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



# Testing S parameter

maf <- Glist$maf[[1]][Glist$rsids[[1]]%in%Glist$rsidsLD[[1]]]
maf_architecture <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    out <- summarise_stblr_maf_architecture(
      fit = fits[[model_name]],
      maf = maf
    )
    out$model <- model_name
    out
  })
)

maf_architecture[
  order(maf_architecture$trait, maf_architecture$selection_s_posthoc),
]

maf_named <- maf
names(maf_named) <- Glist$rsidsLD[[1]]

h_named <- 2 * maf_named * (1 - maf_named)

causal_maf_architecture <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    dm <- as.matrix(fit$dm)
    bm <- as.matrix(fit$bm)
    
    do.call(
      rbind,
      lapply(colnames(dm), function(trait) {
        causal <- causal_by_trait[[trait]]
        causal <- intersect(causal, rownames(dm))
        
        y <- log(dm[causal, trait] * bm[causal, trait]^2 + 1e-12)
        x <- log(h_named[causal])
        
        fit_lm <- lm(y ~ x)
        
        data.frame(
          model = model_name,
          trait = trait,
          selection_s_posthoc_causal = unname(coef(fit_lm)[2]),
          intercept = unname(coef(fit_lm)[1]),
          n_causal = length(causal),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

causal_maf_architecture[
  order(causal_maf_architecture$trait,
        causal_maf_architecture$selection_s_posthoc_causal),
]

causal_maf_architecture <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    dm <- as.matrix(fit$dm)
    bm <- as.matrix(fit$bm)
    
    do.call(
      rbind,
      lapply(colnames(dm), function(trait) {
        causal <- causal_by_trait[[trait]]
        causal <- intersect(causal, rownames(dm))
        
        y <- log(dm[causal, trait] * bm[causal, trait]^2 + 1e-12)
        x <- log(h_named[causal])
        
        fit_lm <- lm(y ~ x)
        sm <- summary(fit_lm)
        
        data.frame(
          model = model_name,
          trait = trait,
          selection_s_posthoc_causal = unname(coef(fit_lm)[2]),
          se = sm$coefficients[2, 2],
          p_value = sm$coefficients[2, 4],
          r2 = sm$r.squared,
          intercept = unname(coef(fit_lm)[1]),
          n_causal = length(causal),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

causal_maf_architecture[
  order(causal_maf_architecture$trait,
        causal_maf_architecture$selection_s_posthoc_causal),
]


aggregate(
  selection_s_posthoc_causal ~ trait,
  data = causal_maf_architecture,
  FUN = function(x) c(
    mean = mean(x),
    sd = sd(x),
    min = min(x),
    max = max(x)
  )
)


maf_architecture_signal <- do.call(
  rbind,
  lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    dm <- as.matrix(fit$dm)
    bm <- as.matrix(fit$bm)
    
    do.call(
      rbind,
      lapply(colnames(dm), function(trait) {
        keep <- dm[, trait] >= 0.01
        keep <- keep & rownames(dm) %in% names(h_named)
        
        y <- log(dm[keep, trait] * bm[keep, trait]^2 + 1e-12)
        x <- log(h_named[rownames(dm)[keep]])
        
        if (length(y) < 5) {
          return(data.frame(
            model = model_name,
            trait = trait,
            selection_s_posthoc_signal = NA_real_,
            intercept = NA_real_,
            n_markers = length(y),
            stringsAsFactors = FALSE
          ))
        }
        
        fit_lm <- lm(y ~ x)
        
        data.frame(
          model = model_name,
          trait = trait,
          selection_s_posthoc_signal = unname(coef(fit_lm)[2]),
          intercept = unname(coef(fit_lm)[1]),
          n_markers = length(y),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

maf_architecture_signal[
  order(maf_architecture_signal$trait,
        maf_architecture_signal$selection_s_posthoc_signal),
]


fit0 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  selection_s = NULL,
  seed = 10
)

fit_minus1 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  selection_s = -1,
  seed = 10
)

max(abs(fit0$dm - fit_minus1$dm))
max(abs(fit0$bm - fit_minus1$bm))
max(abs(fit0$vle - fit_minus1$vle))
max(abs(fit0$vld - fit_minus1$vld))



fit_omit <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  seed = 10
)

fit_null <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  selection_s = NULL,
  seed = 10
)

max(abs(fit_omit$dm - fit_null$dm))
max(abs(fit_omit$bm - fit_null$bm))



fit_s0 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  selection_s = 0,
  seed = 10
)

fit_sneg <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  selection_s = -0.5,
  seed = 10
)

fit_s0$input[c(
  "selection_s",
  "selection_s_fixed",
  "selection_s_exponent",
  "selection_s_scale"
)]

fit_sneg$input[c(
  "selection_s",
  "selection_s_fixed",
  "selection_s_exponent",
  "selection_s_scale"
)]


fits_s <- list(
  bayesC = fit0,
  bayesC_s0 = fit_s0,
  bayesC_sneg05 = fit_sneg
)

topn_power_s <- do.call(
  rbind,
  lapply(names(fits_s), function(model_name) {
    causal_topn_summary(
      fit = fits_s[[model_name]],
      model_name = model_name,
      causal_by_trait = causal_by_trait
    )
  })
)

topn_power_matrix(topn_power_s, top_n_value = 20)
topn_power_matrix(topn_power_s, top_n_value = 50)
topn_power_matrix(topn_power_s, top_n_value = 100)

do.call(
  rbind,
  lapply(names(fits_s), function(model_name) {
    out <- summarise_stblr_maf_architecture(
      fit = fits_s[[model_name]],
      maf = maf,
      markers = causal_by_trait
    )
    out$model <- model_name
    out
  })
)


fit_spos1 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  selection_s = 1,
  seed = 10
)

fit_sneg1_5 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesC",
  selection_s = -1.5,
  seed = 10
)

fits_s_extreme <- list(
  bayesC = fit0,
  bayesC_spos1 = fit_spos1,
  bayesC_sneg15 = fit_sneg1_5
)

do.call(
  rbind,
  lapply(names(fits_s_extreme), function(model_name) {
    out <- summarise_stblr_maf_architecture(
      fit = fits_s_extreme[[model_name]],
      maf = maf,
      markers = causal_by_trait
    )
    out$model <- model_name
    out
  })
)

topn_power_s_extreme <- do.call(
  rbind,
  lapply(names(fits_s_extreme), function(model_name) {
    causal_topn_summary(
      fit = fits_s_extreme[[model_name]],
      model_name = model_name,
      causal_by_trait = causal_by_trait
    )
  })
)

topn_power_matrix(topn_power_s_extreme, top_n_value = 20)
topn_power_matrix(topn_power_s_extreme, top_n_value = 50)
topn_power_matrix(topn_power_s_extreme, top_n_value = 100)



## ------------------------------------------------------------
## Minimal sanity checks for fixed selection_s in CSR BayesR
## ------------------------------------------------------------

fitR0 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = NULL,
  seed = 10
)

fitR_minus1 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = -1,
  seed = 10
)

## selection_s = -1 should reproduce the ordinary BayesR model
## because exponent = selection_s + 1 = 0 and h^0 = 1.

max(abs(fitR0$dm - fitR_minus1$dm))
max(abs(fitR0$bm - fitR_minus1$bm))
max(abs(fitR0$vle - fitR_minus1$vle))
max(abs(fitR0$vld - fitR_minus1$vld))

## BayesR-specific checks
max(abs(fitR0$dm_component_mean - fitR_minus1$dm_component_mean))

max_comp_prob_diff <- max(unlist(
  Map(
    function(a, b) max(abs(a - b)),
    fitR0$comp_prob,
    fitR_minus1$comp_prob
  )
))

max_comp_prob_diff

## ------------------------------------------------------------
## Check omitted selection_s versus explicit NULL
## ------------------------------------------------------------

fitR_omit <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  seed = 10
)

fitR_null <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = NULL,
  seed = 10
)

max(abs(fitR_omit$dm - fitR_null$dm))
max(abs(fitR_omit$bm - fitR_null$bm))
max(abs(fitR_omit$vle - fitR_null$vle))
max(abs(fitR_omit$vld - fitR_null$vld))

max(abs(fitR_omit$dm_component_mean - fitR_null$dm_component_mean))

max_comp_prob_diff_null <- max(unlist(
  Map(
    function(a, b) max(abs(a - b)),
    fitR_omit$comp_prob,
    fitR_null$comp_prob
  )
))

max_comp_prob_diff_null

## ------------------------------------------------------------
## Run fixed-S BayesR models
## ------------------------------------------------------------

fitR_s0 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = 0,
  seed = 10
)

fitR_sneg <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = -0.5,
  seed = 10
)

fitR_s0$input[c(
  "selection_s",
  "selection_s_fixed",
  "selection_s_exponent",
  "selection_s_scale"
)]

fitR_sneg$input[c(
  "selection_s",
  "selection_s_fixed",
  "selection_s_exponent",
  "selection_s_scale"
)]

## ------------------------------------------------------------
## Check BayesR component-probability consistency
## ------------------------------------------------------------

check_bayesr_component_consistency <- function(fit) {
  out <- lapply(names(fit$comp_prob), function(trait) {
    cp <- fit$comp_prob[[trait]]
    
    data.frame(
      trait = trait,
      max_row_sum_error = max(abs(rowSums(cp) - 1)),
      max_dm_error = max(abs(fit$dm[, trait] - (1 - cp[, "component_0"]))),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

check_bayesr_component_consistency(fitR0)
check_bayesr_component_consistency(fitR_minus1)
check_bayesr_component_consistency(fitR_s0)
check_bayesr_component_consistency(fitR_sneg)

## ------------------------------------------------------------
## Compare causal recovery
## ------------------------------------------------------------

fits_R_s <- list(
  bayesR = fitR0,
  bayesR_s0 = fitR_s0,
  bayesR_sneg05 = fitR_sneg
)

topn_power_R_s <- do.call(
  rbind,
  lapply(names(fits_R_s), function(model_name) {
    causal_topn_summary(
      fit = fits_R_s[[model_name]],
      model_name = model_name,
      causal_by_trait = causal_by_trait
    )
  })
)

topn_power_matrix(topn_power_R_s, top_n_value = 20)
topn_power_matrix(topn_power_R_s, top_n_value = 50)
topn_power_matrix(topn_power_R_s, top_n_value = 100)


## ------------------------------------------------------------
## Compare post-hoc MAF architecture on known causal markers
## ------------------------------------------------------------

do.call(
  rbind,
  lapply(names(fits_R_s), function(model_name) {
    out <- summarise_stblr_maf_architecture(
      fit = fits_R_s[[model_name]],
      maf = maf,
      markers = causal_by_trait
    )
    
    out$model <- model_name
    out
  })
)

## ------------------------------------------------------------
## Optional more extreme fixed-S values
## ------------------------------------------------------------

fitR_spos1 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = 1,
  seed = 10
)

fitR_sneg1_5 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = -1.5,
  seed = 10
)

fits_R_s_extreme <- list(
  bayesR = fitR0,
  bayesR_spos1 = fitR_spos1,
  bayesR_sneg15 = fitR_sneg1_5
)

do.call(
  rbind,
  lapply(names(fits_R_s_extreme), function(model_name) {
    out <- summarise_stblr_maf_architecture(
      fit = fits_R_s_extreme[[model_name]],
      maf = maf,
      markers = causal_by_trait
    )
    
    out$model <- model_name
    out
  })
)

topn_power_R_s_extreme <- do.call(
  rbind,
  lapply(names(fits_R_s_extreme), function(model_name) {
    causal_topn_summary(
      fit = fits_R_s_extreme[[model_name]],
      model_name = model_name,
      causal_by_trait = causal_by_trait
    )
  })
)

topn_power_matrix(topn_power_R_s_extreme, top_n_value = 20)
topn_power_matrix(topn_power_R_s_extreme, top_n_value = 50)
topn_power_matrix(topn_power_R_s_extreme, top_n_value = 100)

max(abs(fitR0$dm - fitR_minus1$dm))
max(abs(fitR0$bm - fitR_minus1$bm))
max(abs(fitR0$vle - fitR_minus1$vle))
max(abs(fitR0$vld - fitR_minus1$vld))
max_comp_prob_diff


fitR0 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = NULL,
  seed = 10
)

fitR_minus1 <- stblr_csr(
  stats = stats,
  Glist = Glist,
  method = "bayesR",
  selection_s = -1,
  seed = 10
)

max(abs(fitR0$dm - fitR_minus1$dm))
max(abs(fitR0$bm - fitR_minus1$bm))
max(abs(fitR0$vle - fitR_minus1$vle))
max(abs(fitR0$vld - fitR_minus1$vld))
max(abs(fitR0$dm_component_mean - fitR_minus1$dm_component_mean))

max(unlist(Map(
  function(a, b) max(abs(a - b)),
  fitR0$comp_prob,
  fitR_minus1$comp_prob
)))

check_bayesr_component_consistency <- function(fit) {
  do.call(
    rbind,
    lapply(names(fit$comp_prob), function(trait) {
      cp <- fit$comp_prob[[trait]]
      
      data.frame(
        trait = trait,
        max_row_sum_error = max(abs(rowSums(cp) - 1)),
        max_dm_error = max(abs(fit$dm[, trait] - (1 - cp[, "component_0"]))),
        stringsAsFactors = FALSE
      )
    })
  )
}

check_bayesr_component_consistency(fitR0)
check_bayesr_component_consistency(fitR_minus1)
