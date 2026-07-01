# ST-BLR Backend Naming

This note records the naming convention for the supported CSR and BED BLR
backends. It is a maintenance guide only; it does not change sampler math,
posterior output semantics, or public user-facing behavior.

## Public Interfaces

Supported user-facing entry points:

- `stblr_csr(..., method = "bayesC")`
- `stblr_csr(..., method = "bayesR")`
- `stblr_csr_annot(..., annotation_model = "prior")`
- `stblr_csr_annot(..., annotation_model = "learned")`
- `stblr_csr_annot(..., annotation_model = "group")`
- `stblr_csr_annot(..., annotation_model = "sbayesrc")`
- `stblr_bed(..., method = "bayesC")`
- `stblr_bed(..., method = "bayesR")`

The explicit public convenience wrapper `stblr_csr_bayesr()` remains supported.
The older public `stblr_bed_marker()` remains available for direct BED marker
workflows.

## Internal R Helpers

Preferred formatter names:

- `.format_stblr_csr_bayesc_fit()`
- `.format_stblr_csr_bayesr_fit()`
- `.format_stblr_bed_bayesc_fit()`
- `.format_stblr_bed_bayesr_fit()`

Preferred fit-helper names:

- `.fit_stblr_csr_bayesc()`
- `.fit_stblr_csr_bayesr()`
- `.fit_stblr_bed_bayesc()`
- `.fit_stblr_bed_bayesr()`

Older internal helper names are retained as compatibility aliases when they may
be used by tests or local scripts. For example,
`.stblr_csr_bayesr_experimental()` delegates to `stblr_csr_bayesr()`, and
`.stblr_bed_marker_bayesr_experimental()` delegates to
`.fit_stblr_bed_bayesr()`.

## Native Backend Naming

The current first pass keeps C++ filenames and Rcpp-exported native symbols
unchanged to avoid churn in generated files and native build artifacts.

Current native files:

- `src/st_cpg_omp_csr.cpp`
- `src/st_cpg_omp_csr_scheduled.cpp`
- `src/st_cpg_omp_csr_bayesr.cpp`
- `src/st_cpg_omp_individual_scheduled_chains.cpp`
- `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp`

Possible future native filenames:

- `src/stblr_csr_bayesc.cpp`
- `src/stblr_csr_scheduled_bayesc.cpp`
- `src/stblr_csr_bayesr.cpp`
- `src/stblr_bed_bayesc.cpp`
- `src/stblr_bed_bayesr.cpp`

Any future native rename should be a separate mechanical change with regenerated
Rcpp exports, compatibility aliases where needed, and focused tests.

## Fit Metadata

Supported fit objects should expose:

- `fit$input$method`
- `fit$input$model`
- `fit$input$backend`
- `fit$input$data_level`
- `fit$input$scheduled`
- `fit$input$nchains`
- `fit$input$keep_chains` when applicable
- `fit$input$updateE`
- `fit$input$updateLDswap` when applicable

Backend values:

- `csr_bayesc`
- `csr_scheduled_bayesc`
- `csr_bayesr`
- `csr_prior_bayesc`
- `csr_annot_bayesc`
- `csr_group_bayesc`
- `csr_sbayesrc`
- `bed_bayesc`
- `bed_bayesr`

Data-level values:

- `summary` for CSR summary-statistics fits
- `individual` for BED individual-level fits

BayesR output convention is unchanged: `fit$dm` is
`P(component > 0)`, `fit$comp_prob` stores marker-by-component probabilities,
and `fit$dm_component_mean` stores posterior mean component index.

BayesC-like annotation-aware CSR backends now follow the standard native
multi-chain output convention. `csr_prior_bayesc`, `csr_annot_bayesc`, and
`csr_group_bayesc` expose `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, and
`bm_max` when multiple chains are requested or compact chains are kept.
`keep_chains = TRUE` returns `fit$chains[[trait]][[chain]]$dm` and `$bm` as
named marker vectors, with compact annotation-specific chain fields where
available.

`csr_prior_bayesc`, `csr_group_bayesc`, and `csr_annot_bayesc` also support
optional active/null LD-swap/MH through the standard `updateLDswap` and
`ld_swap_*` controls. Aggregate diagnostics are returned as `fit$ld_swap`;
when compact chains are kept, chain diagnostics are available as
`fit$ld_swap_chains` and `fit$chains[[trait]][[chain]]$ld_swap`.

## Compatibility Aliases

Compatibility aliases are intentionally kept for local scripts and historical
tests. They should be marked in comments as compatibility aliases and should
delegate to the current helper or public wrapper without changing behavior.

Current intentionally retained older names include:

- `.format_stblr_fit()`
- `.format_stblr_bayesr_fit()`
- `.stblr_csr_bayesr_experimental()`
- `.stblr_bed_marker_bayesr_experimental()`
- native symbols containing `st_cpg`, `stblr_cpg`, `individual_scheduled`, or
  `bed_marker`

## Remaining Unsupported Features

These limitations remain explicit:

- scheduled CSR BayesR is not implemented
- active/active BayesR LD-swap is not implemented
- LD-swap/MH for `csr_sbayesrc` is not implemented
- marker-specific or annotation-specific BayesR LD-swap priors are not
  implemented
- BED `chain_seeds` are not supported
- BED `covar` currently requires pre-adjusted phenotypes

## Annotation-Aware BLR Alignment Plan

Annotation-aware BLR should align with the same public interface, metadata, and
output conventions as CSR and BED BLR. The detailed audit and staged alignment
plan is in `docs/dev/stblr_annotation_backend_design.md`.

The clean `stblr_csr_annot()` interface is the public annotation-aware CSR entry
point. Existing exported annotation wrappers remain compatibility entry points
and should continue returning standardized metadata and annotation aliases.
Direct `stblr_csr(..., annotations = ...)` dispatch can be added later after the
annotation argument contract is stable. Annotation-aware changes should not
rename native C++ files or symbols.
