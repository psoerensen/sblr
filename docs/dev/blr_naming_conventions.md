# BLR naming conventions

## Identifiers and functions

Machine identifiers are lowercase: families `stblr` and `mtblr`; models
`bayesc`, `sbayesc`, `bayesr`, `sbayesr`, `bayesrc`, and `sbayesrc`; operators `csr`, `block_eigen`,
`packed_bed`, and `dense_reference`. `cpg` is an internal historical kernel
term, never a public model. Public fitting functions are exactly
`stblr_csr()`, `stblr_block_eigen()`, `stblr_bed()`, `mtblr_csr()`,
`mtblr_block_eigen()`, and `mtblr_bed()`, plus the single specialized
`stblr_csr_annot()` entry while its annotation models require materially
different inputs.

Methods accept one exact lowercase spelling. Case conversion and mixed-case
aliases are forbidden.

The `s` prefix denotes the summary-statistics data level. Thus `bayesc`,
`bayesr`, and `bayesrc` are individual-level models, while `sbayesc`,
`sbayesr`, and `sbayesrc` are their summary-statistics counterparts. The
paired names reuse the same prior kernel; the prefix never activates MAF
scaling. `selection_s` independently selects `maf_s` for BayesC or
`component_maf_s` for BayesR/BayesRC. The probability-policy identifiers are
`global`, `fixed_marker`, `group`, `learned_logistic`, and
`annotation_probit_stick`.

New fits record `model_semantics_version = 2` and
`model_semantics = "s_prefix_means_summary_statistics"`, plus separate
`prior_kernel`, `data_level`, and `effect_scale_policy` fields. Objects lacking
that marker are not silently reinterpreted.

## Controls

The common terminal control block is `nit`, `nburn`, `nthin`, `seed`,
`nchains`, `ncores`, `chain_seeds`, `keep_chains`, `convergence`,
`convergence_control`, `memory_warning_gb`, and `verbose`, with defaults 1000,
500, 1, 1, 1L, 1L, NULL, FALSE, `c("auto", "none", "core")`, NULL, 8, and
FALSE. Unsupported combinations fail explicitly rather than silently changing
meaning.

`ncores` always means requested concurrent logical MCMC tasks. `nthreads`
always means non-chain preparation or decoding work. `seed` is the fit-local
base; `chain_seeds_requested`, `chain_seeds_resolved`, and
`task_seeds_resolved` distinguish requests from deterministic resolution.

## Results

Final values, posterior means, and traces use `_final`, `_mean`, and `_trace`.
Probability fields are `pi_final`, `pi_mean`, and `pi_trace`; marker mixture
probabilities are `component_probabilities`. MT covariance fields are
`cov_b_mean`, `cov_g_mean`, `cov_e_mean`, `cov_b_final`, `cov_g_final`, and
`cov_e_final`. Across-chain summaries are `bm_chain_mean_sd|min|max` and
`dm_chain_mean_sd|min|max`; they are not posterior standard deviations.

Classes include `blr_fit` and exactly one of `stblr_fit` or `mtblr_fit`.
Backend details live in `diagnostics`; analytical memory is never described as
measured RSS.

Trait-level traces are exactly `vbs`, `vgs`, `ves`, `vle`, and `vld`, with
iteration × trait orientation. Formal convergence quantities use those same
names. Marker effects use `b_final`/`bm` for final/posterior-mean effective
effects and `d_final`/`dm` for final state/posterior non-null probability.
Latent effects, when defined, use `beta_final`/`beta_mean`. Full MT covariance
matrices remain `cov_b_*`, `cov_g_*`, and `cov_e_*` rather than aliases of the
five traces.

Every fit data object owns explicit `genotype_scale`, `effect_scale`,
`phenotype_scale`, `ld_scale`, `n_total`, `n_used`, and `n_by_trait` metadata.
When MAF scaling is requested it also records `selection_maf_source`,
`selection_maf_alignment_status`, and `selection_maf_fallback_used`.
`keep_chains` retains compact logical-chain records;
`convergence_control$keep_traces` independently retains convergence arrays.

## MT BayesR state names

`mixture_var` is the one public component-multiplier control for STBLR and
MTBLR. MT joint states are named `null` and
`<trait-pattern>__component_<positive-index>`. `component_final` is zero-based
with zero reserved for null; `component_probabilities` includes component zero.
`pi_final`, `pi_mean`, and `pi_trace` name the joint-state probability vector.
No `pi`, `pim`, `pis`, `comp_prob`, or case-variant method aliases are added.
