#' Simulate Multi-Trait Phenotypes with Annotation Enrichment
#'
#' Simulates multi-trait phenotypes while enriching causal-marker sampling in
#' selected, potentially overlapping marker annotations. It also constructs
#' marker-specific inclusion probabilities and variance multipliers for
#' annotation-aware ST-BLR examples.
#'
#' @param Glist Optional qgg genotype list.
#' @param chr,rsids,ids Optional chromosome, marker, and individual selections.
#' @param causal_rsids Markers eligible to be sampled as causal.
#' @param W Optional genotype matrix.
#' @param n,m Simulated sample and marker counts when `W` is omitted.
#' @param nt Number of traits.
#' @param n_shared,n_specific Numbers of shared and trait-specific causal
#'   markers.
#' @param h2 Trait heritabilities.
#' @param rg,re Genetic-effect and residual correlation specifications.
#' @param effect_sd Standard deviation used to simulate marker effects.
#' @param maf_min,maf_max Minor-allele-frequency range when simulating `W`.
#' @param standardize_W Standardize marker columns before simulation.
#' @param seed Optional simulation seed.
#' @param exact_shared_cor Force the sampled shared effects to have exactly
#'   correlation `rg`.
#' @param annot Optional marker annotation matrix.
#' @param sets Optional list of marker indices per annotation.
#' @param n_annotations Number of annotations to simulate when `annot` is
#'   omitted.
#' @param annotation_prob Annotation membership probability.
#' @param force_at_least_one_annotation Ensure every marker belongs to at least
#'   one simulated annotation.
#' @param enriched_annotations Annotation indices enriched for causal markers.
#' @param annotation_enrichment Causal-marker sampling weight multiplier.
#' @param base_pi Baseline marker inclusion probability.
#' @param enriched_traits Traits receiving annotation-derived prior enrichment.
#' @param enriched_pi_multiplier,enriched_vb_multiplier Annotation-specific
#'   prior multipliers.
#' @param pi_cap,vb_cap Bounds for marker-specific priors.
#' @param center_log_pi,center_log_vb Center annotation-derived log-priors.
#'
#' @return A simulation list containing phenotypes, effects, annotations,
#'   causal-marker summaries, and marker-specific priors.
#' @name mtsim_annotation
#' @export mtsim_annotation
NULL

#' Summarize Annotation Signal
#'
#' Summarizes annotation sizes and causal-marker enrichment. When a fitted
#' model is supplied, posterior inclusion and absolute effect summaries are
#' added for each trait.
#'
#' @param sim Simulation output from [mtsim_annotation()].
#' @param fit Optional fitted annotation-aware ST-BLR model.
#'
#' @return A data frame with annotation-level signal summaries.
#' @name summarize_annotation_signal
#' @export summarize_annotation_signal
NULL
