# ST-BLR Backend Naming

## Experimental packed-BED BayesC names

The canonical public name is `bed_bayesc`, reached through
`stblr_bed(method = "bayesc")`. The lower-level native names
`stblr_cpg_omp_bed_marker_scheduled` and
`stblr_cpg_omp_bed_marker_sparse` are explicitly experimental. Their returned
backend identifiers (`bed_scheduled_bayesc` and `bed_sparse_bayesc`) describe
research implementations and must not be interpreted as canonical aliases.

This note records the naming convention for the supported CSR and BED BLR
backends. It is a maintenance guide only; it does not change sampler math,
posterior output semantics, or public user-facing behavior.

## Public Interfaces

Supported user-facing data preparation entry points:

- `make_summary_stats()`
- `make_sparse_ld()`

Supported user-facing model fitters:

- `stblr_csr(..., method = "bayesC")`
- `stblr_csr(..., method = "bayesR")`
- `stblr_csr_annot(..., annotation_model = "prior")`
- `stblr_csr_annot(..., annotation_model = "learned")`
- `stblr_csr_annot(..., annotation_model = "group")`
- `stblr_csr_annot(..., annotation_model = "sbayesrc")`
- `stblr_bed(..., method = "bayesC")`
- `stblr_bed(..., method = "bayesR")`

Supported user-facing posterior summaries and diagnostics:

- `summarise_posterior()`
- `plot_posterior()`
- `summarise_components()`
- `summarise_architecture()`
- `make_credible_sets()`
- `check_stblr_consistency()`

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
- `fit$input$annotations`
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

BayesR and SBayesRC-like output convention is unchanged: `fit$dm` is
`P(component > 0)`, `fit$comp_prob` stores marker-by-component probabilities,
and `fit$dm_component_mean` stores posterior mean zero-based component index
where the quantity is returned or derivable. CSR BayesR uses `component_0` for
the null component column and `dm = 1 - P(component_0)`.
CSR SBayesRC names components by gamma values; its null component column is `gamma_0.00` and
`dm = 1 - P(gamma_0.00)`.

BayesC-like annotation-aware CSR backends now follow the standard native
multi-chain output convention. `csr_prior_bayesc`, `csr_annot_bayesc`, and
`csr_group_bayesc` expose `dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, and
`bm_max` when multiple chains are requested or compact chains are kept.
`keep_chains = TRUE` returns `fit$chains[[trait]][[chain]]$dm` and `$bm` as
named marker vectors, with compact annotation-specific chain fields where
available.

`csr_prior_bayesc`, `csr_group_bayesc`, `csr_annot_bayesc`, and
`csr_sbayesrc` also support optional active/null LD-swap/MH through the
standard `updateLDswap` and `ld_swap_*` controls. Aggregate diagnostics are
returned as `fit$ld_swap`; when compact chains are kept, chain diagnostics are
available as `fit$ld_swap_chains` and
`fit$chains[[trait]][[chain]]$ld_swap`. The SBayesRC move relocates the full
active `(component, b)` state to a null LD neighbor and uses current
annotation-dependent component probabilities in the MH ratio.

Annotation-aware CSR models return `vle` and `vld` using the same definitions
and conventions as annotation-unaware CSR models. The formatted fields are
iteration-by-trait trace matrices with trait column names.

The cross-backend computation and return inventory is maintained in
`docs/dev/stblr_backend_computation_inventory.md`.

## BayesS-Style Selection-S Terminology

Use `selection_s` for a BayesS-style global MAF-dependent marker-effect
variance scaling parameter. This is distinct from SBayesRC annotation-selection
coefficients such as `alpha`, `eta_pi`, and annotation-dependent component
probabilities.

Fixed-S support is limited to `csr_bayesc`
(`stblr_csr(method = "bayesC", scheduled = FALSE)`), `csr_bayesr`
(`stblr_csr(method = "bayesR")` or `stblr_csr_bayesr()`), and `csr_sbayesrc`
(`stblr_csr_annot(annotation_model = "sbayesrc")`). Sampled trait-specific
`selection_s` support is limited to unscheduled `csr_bayesc`, `csr_bayesr`,
and annotation-aware `csr_sbayesrc`, with default prior `c(-3, 2)` and default
random-walk MH proposal SD 0.35. The `csr_scheduled_bayesc` backend,
BayesC-like annotation-aware CSR backends
(`csr_prior_bayesc`, `csr_annot_bayesc`, and `csr_group_bayesc`), and BED
backends do not support sampled `selection_s`.

Selection-S support summary:

| Backend | Fixed `selection_s` | Sampled `selection_s` |
| --- | --- | --- |
| `csr_bayesc` | yes | yes |
| `csr_bayesr` | yes | yes |
| `csr_sbayesrc` | yes | yes |
| `csr_scheduled_bayesc` | no | no |
| `csr_prior_bayesc` | no | no |
| `csr_annot_bayesc` | no | no |
| `csr_group_bayesc` | no | no |
| `bed_bayesc` | no | no |
| `bed_bayesr` | no | no |

Preferred public argument names are:

- `selection_s = NULL`
- `estimate_selection_s = FALSE`
- `selection_s_init = 0`
- `selection_s_prior = c(-3, 2)`
- `selection_s_proposal_sd = 0.35`

Preferred future fit metadata fields are:

- `fit$input$selection_s`
- `fit$input$estimate_selection_s`
- `fit$input$selection_s_scale`
- `fit$input$selection_s_exponent`
- `fit$selection_s`
- `fit$selection_s_sd`
- `fit$selection_s_min`
- `fit$selection_s_max`
- `fit$selection_s_trace`
- `fit$selection_s_acceptance`

`fit$selection_s_trace` is an iteration x trait matrix. `fit$selection_s` is
the posterior mean by trait, and `fit$selection_s_acceptance` is the MH
acceptance rate by trait. With `keep_chains = TRUE`, compact chain output uses
`fit$chains[[trait]][[chain]]$selection_s` and
`fit$chains[[trait]][[chain]]$selection_s_acceptance`.

For the standard CSR path, fitted `b`/`bm` values are
standardized-genotype-scale effects. A BayesS allele-scale prior with
heterozygosity exponent `S` therefore maps to a standardized-effect prior with
exponent `S + 1`. For CSR BayesC, BayesR, and SBayesRC, fixed sampler-level
`selection_s` scales standardized-genotype effect prior variances by
`h^(selection_s + 1)`, where `h = 2p(1-p)`. Fixed support applies the same
marker-specific variance factor consistently in conditional effect updates,
prior-density or component-probability calculations, active/null LD-swap/MH
prior terms, and marker-effect variance updates.

Sampled CSR BayesC estimates one `S_t` per trait and per chain using the
active-marker log posterior contribution:

```text
log p(S | b, d, vb)
= log p(S)
- 0.5 * sum_{j: d_j = 1} [
    log(q_j(S)) +
    b_j^2 / (vb * q_j(S))
  ]
```

For current CSR standardized effects, `q_j(S) = h_j^(S + 1)`. A random-walk MH
proposal uses `S_new = S_current + Normal(0, selection_s_proposal_sd)` with a
uniform prior over `selection_s_prior`. For sampled `selection_s`, the default
prior is Uniform(-3, 2) and the default proposal SD is 0.35. These tuning
arguments only affect `estimate_selection_s = TRUE`; they do not affect
ordinary BayesC/BayesR/SBayesRC or fixed `selection_s`.

Sampled CSR BayesR uses the active non-null component contribution:

```text
log p(S | b, gamma, vb)
= log p(S)
- 0.5 * sum_{j: gamma_j > 0} [
    log(q_j(S)) +
    b_j^2 / (vb * gamma_j * q_j(S))
  ]
```

where `gamma_j` is the current non-null BayesR component variance multiplier
and `q_j(S) = h_j^(S + 1)`.

Sampled CSR SBayesRC uses the same active non-null component contribution:

```text
log p(S | b, gamma, vb)
= log p(S)
- 0.5 * sum_{j: gamma_j > 0} [
    log(q_j(S)) +
    b_j^2 / (vb * gamma_j * q_j(S))
  ]
```

For SBayesRC, annotations affect component probabilities and alpha updates;
`selection_s` affects marker-specific effect-size prior variance.

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
- annotation-specific BayesR LD-swap priors are not implemented
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

## Phase 17A status vocabulary

Current names are retained unchanged. Scalar CSR names denote canonical public
routes. Names containing `block_eigen` denote internal experimental operator
routes. `sblr()` remains the public legacy multivariate wrapper; `mtblr`,
`mtblr_cpg`, `mtblr_cpg_arma`, `mtblr_cpg_omp`, and `mtblr_eigen` identify its
native algorithm variants. `mtblr_hybrid` and `mtblr_cpg_omp_csr` are
native-only, not additional public APIs. These are support classifications, not
rename or compatibility commitments.

## Phase 17B multivariate contract names

`mtblr` is authoritative for public `sblr(algorithm = "default")`. Phase 17B
“raw” fixtures mean its frozen 20-position native legacy result, not
`stblr_raw_v1`; the public object is the named legacy fit. `Mt*` structures in
`blr_mt_default_audit_types.h` are audit-only vocabulary, not production types.
Alternative names retain their explicit variant dispositions and do not imply
equivalent models or supported schemas.

## Phase 17C corrected multivariate contract names

`mtblr` remains the authoritative supported public legacy, noncanonical route.
“Phase 17C raw reference” means its corrected 20-position native legacy result,
not `stblr_raw_v1`. “Retained count” is accumulator-specific: marker summaries
use the post-burn-relative thinned count, while covariance and updated
probability summaries use their own post-burn contribution counts. Phase 17B
names and fixtures denote historical pre-correction evidence.
