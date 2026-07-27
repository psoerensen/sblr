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

S models are scale policies, not copied kernels: `sbayesc` is BayesC plus
`maf_s`, `sbayesr` is BayesR plus `component_maf_s`, and `sbayesrc` is the
annotation-probit-stick mixture plus `component_maf_s`. The probability-policy
identifiers are `global`, `fixed_marker`, `group`, `learned_logistic`, and
`annotation_probit_stick`.

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
`keep_chains` retains compact logical-chain records;
`convergence_control$keep_traces` independently retains convergence arrays.
