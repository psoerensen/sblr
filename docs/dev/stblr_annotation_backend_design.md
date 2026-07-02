# Annotation-Aware ST-BLR Backend Alignment

## Executive Summary

The annotation-aware implementations are all CSR summary-statistics models that
use disk-backed sparse LD. They are already close to the standard ST-BLR object
shape because the first fields are `bm`, `dm`, and the usual variance traces.
They now follow the current public CSR/BED multi-chain conventions for the
annotation-aware CSR SBayesRC backend and the three BayesC-like annotation
backends: fixed-prior, learned-annotation, and group-prior BayesC. These
BayesC-like backends now also support optional active/null LD-swap/MH for
comparison with annotation-unaware CSR BayesC. The SBayesRC backend supports
optional active/null LD-swap/MH for comparison with annotation-unaware CSR
BayesR.

Annotation-aware CSR models return `vle` and `vld` using the same definitions
and conventions as annotation-unaware CSR models. Native slot 20 stores the
linkage-equilibrium marker-effect component trace, and native slot 21 stores
the LD contribution trace computed as `vgs - vle`.

Do not rename C++ files or native symbols in the alignment phase. The first
R-side alignment step adds a dedicated `stblr_csr_annot()` entry point and
standardized annotation-aware metadata/output aliases while preserving existing
explicit wrappers. Direct `stblr_csr(..., annotations = ...)` dispatch remains a
later design decision.

## C++ Backend Inventory

### `src/st_cpg_omp_csr_prior.cpp`

- Exported Rcpp function: `stblr_cpg_omp_csr_prior()`.
- Model: BayesC-like CSR ST-BLR with fixed marker-specific prior inclusion
  probabilities and/or marker-specific prior variance multipliers.
- Data level: summary statistics with disk-backed sparse LD.
- Annotation/prior input: no annotation matrix is required by C++; R can derive
  marker priors from annotations before calling it. Native inputs are
  `use_pi_marker`, `pi_marker`, `use_vb_multiplier`, and `vb_multiplier`.
- Required core inputs: `wy`, `ww`, `yy`, initial `b`/`d`/`r` state flags,
  `ld_prefix`, `B`, `E`, `ssb_prior`, `sse_prior`, global two-state `pi`,
  `nub`, `nue`, `adjE`, per-trait `n`, `nit`, `nburn`, `nthin`,
  `pi_prior_a`, `pi_prior_b`, `ncores`, and `seed`.
- MCMC controls: `nit`, `nburn`, `nthin`, `updateB`, `updateE`, `updatePi`,
  `ncores`, and `seed`.
- Supports native multi-chain controls: `nchains`, `keep_chains`, and
  `chain_seeds`. Returned chain summaries include `bm_sd`, `bm_min`, `bm_max`,
  `dm_sd`, `dm_min`, and `dm_max`; compact chains expose per-trait/per-chain
  named `dm` and `bm` vectors after R formatting.
- Supports optional LD-swap/MH controls: `updateLDswap`, `ld_swap_prob`,
  `ld_swap_r2`, `ld_swap_max_friends`, and `ld_swap_moves`. The MH ratio uses
  the plain CSR BayesC likelihood/proposal terms plus marker-specific inclusion
  and effect-variance prior ratios when fixed marker priors are active.
- Unsupported compared with current CSR/BED main interfaces: scheduled updates.
- Return layout: 23 slots: `bm`, `dm`, `wy`, `r`, `b`, `d`, marker index,
  `vbs`, `vgs`, `ves`, `covb`, `covg`, `cove`, `vb`, `vg`, `ve`, final `pi`,
  posterior mean `pim`, reserved diagnostics, `nsamples`/`n`, `vle`, `vld`,
  and LD-swap diagnostics. Multi-chain summary/chain slots follow the standard
  CSR convention.
- ST-BLR resemblance: high for base fields; missing current metadata and chain
  fields at the R wrapper level.

### `src/st_cpg_omp_csr_annot.cpp`

- Exported Rcpp function: `stblr_cpg_omp_csr_annot()`.
- Model: BayesC-like CSR ST-BLR with learned annotation effects on
  marker-specific inclusion probabilities and optional variance multipliers.
- Data level: summary statistics with disk-backed sparse LD.
- Annotation input: dense marker x annotation matrix `A`; overlapping
  annotations are allowed.
- Annotation parameters: `eta_pi_init`, `eta_vb_init`, `sigma_eta_pi`,
  `sigma_eta_vb`, `rw_sd_eta_pi`, `rw_sd_eta_vb`, `annot_update_every`,
  `pi_min`, `pi_max`, `vb_multiplier_min`, and `vb_multiplier_max`.
- Required core inputs: same CSR sufficient-statistics, sparse-LD, state,
  variance-prior, global `pi`, `n`, iteration, threading, and seed inputs as
  the fixed-prior backend.
- MCMC controls: `nit`, `nburn`, `nthin`, `updateB`, `updateE`, `updatePi`,
  `learn_pi_annot`, `learn_vb_annot`, `annot_update_every`, `ncores`, and
  `seed`. Annotation effects are updated by random-walk MH.
- Supports native multi-chain controls: `nchains`, `keep_chains`, and
  `chain_seeds`. Returned chain summaries include `bm_sd`, `bm_min`, `bm_max`,
  `dm_sd`, `dm_min`, and `dm_max`; compact chains expose per-trait/per-chain
  named `dm` and `bm` vectors after R formatting, plus compact `eta_pi` and
  `eta_vb` vectors.
- Supports optional LD-swap/MH controls: `updateLDswap`, `ld_swap_prob`,
  `ld_swap_r2`, `ld_swap_max_friends`, and `ld_swap_moves`. The MH ratio uses
  the current marker-specific inclusion probabilities and current
  marker-specific variance multipliers implied by the current learned
  annotation effects at the iteration where the swap is proposed.
- Unsupported compared with current CSR/BED main interfaces: scheduled updates.
- Return layout: 23 slots. Slots 0-17 match the BayesC-like CSR convention;
  slot 18 is posterior mean `eta_pi`, slot 19 is posterior mean `eta_vb`, slot
  20 is `vle`, slot 21 is `vld`, and slot 22 is LD-swap diagnostics.
- ST-BLR resemblance: high for base `bm`/`dm`/variance fields; annotation
  outputs are native but not yet exposed under a standard `annotation_*`
  naming convention.

### `src/st_cpg_omp_csr_group.cpp`

- Exported Rcpp function: `stblr_cpg_omp_csr_group_annot()`.
- Model: group-prior BayesC-like CSR ST-BLR. Each marker belongs to exactly one
  group, with group-specific inclusion probabilities and optional group
  variance multipliers.
- Data level: summary statistics with disk-backed sparse LD.
- Annotation input: zero-based length-m `group_index` and `ngroup`; R wrappers
  translate user group labels to this format.
- Group-prior inputs: `group_pi_init`, `pi_group_prior_a`,
  `pi_group_prior_b`, `group_vb_multiplier_init`, `updateGroupVb`,
  `nub_group`, `ssb_group_prior`, and `normalize_group_vb`.
- Required core inputs: same CSR sufficient-statistics, sparse-LD, state,
  variance-prior, global `pi`, `n`, iteration, threading, and seed inputs as
  the BayesC-like CSR backends.
- MCMC controls: `nit`, `nburn`, `nthin`, `updateB`, `updateE`, `updatePi`,
  `updateGroupVb`, `ncores`, and `seed`.
- Supports native multi-chain controls: `nchains`, `keep_chains`, and
  `chain_seeds`. Returned chain summaries include `bm_sd`, `bm_min`, `bm_max`,
  `dm_sd`, `dm_min`, and `dm_max`; compact chains expose per-trait/per-chain
  named `dm` and `bm` vectors after R formatting, plus compact group-level
  `group_pi`, `group_vb_multiplier`, and `group_nincluded` vectors.
- Supports optional LD-swap/MH controls: `updateLDswap`, `ld_swap_prob`,
  `ld_swap_r2`, `ld_swap_max_friends`, and `ld_swap_moves`. The MH ratio uses
  the current group inclusion probabilities and current group variance
  multipliers at the iteration where the swap is proposed.
- Unsupported compared with current CSR/BED main interfaces: scheduled updates.
- Return layout: 27 slots. Slots 0-21 follow the BayesC-like CSR convention;
  slots 22-25 are `group_pi`, `group_vb_multiplier`, `group_nincluded`, and
  `group_size`; slot 26 is LD-swap diagnostics. Multi-chain summary/chain
  slots follow the standard CSR convention.
- ST-BLR resemblance: high for base fields; group outputs need standard
  annotation naming aliases or summaries.

### `src/st_sbayesrc_omp_csr.cpp`

- Exported Rcpp function: `stblr_cpg_omp_csr_sbayesrc()`.
- Model: SBayesRC-style annotation-dependent BayesR prior over mixture
  variance multipliers `gamma`. It is BayesR-like/SBayesRC-like rather than
  BayesC-like.
- Data level: summary statistics with disk-backed sparse LD.
- Annotation input: dense marker x annotation matrix `A`; overlapping
  annotations are handled through probit stick-breaking regressions.
- Annotation parameters: `gamma`, `alpha_init`, `sigmaSqAlpha_init`,
  `intercept_flat`, `sigmaSqAlpha_a`, `sigmaSqAlpha_b`, `pi_floor`,
  `updateAlpha`, and `alpha_update_every`.
- Required core inputs: CSR sufficient statistics, sparse LD, initial `b`,
  optional mixture component state, optional residual state, variance priors,
  `n`, iteration, threading, and seed inputs.
- MCMC controls: `nit`, `nburn`, `nthin`, `updateAlpha`, `updateB`, `updateE`,
  `alpha_update_every`, `ncores`, and `seed`. There is no `updatePi`; mixture
  probabilities are induced by `A %*% alpha`.
- Supports native multi-chain controls: `nchains`, `keep_chains`, and
  `chain_seeds`, including component-probability chain output.
- Supports optional LD-swap/MH controls: `updateLDswap`, `ld_swap_prob`,
  `ld_swap_r2`, `ld_swap_max_friends`, and `ld_swap_moves`. The move relocates
  the full active `(component, b)` state to a null LD neighbor. The MH ratio
  uses the current annotation-dependent component probabilities implied by the
  current `alpha`; the effect-prior variance term cancels because the component
  and effect value move together and `vb * gamma[component]` is not
  marker-specific.
- Unsupported compared with current CSR/BED main interfaces: scheduled updates.
- Return layout: 25 slots. Slots 0-17 follow the BayesR-like CSR convention
  with `dm = P(component > 0)`; slot 18 is posterior mean `alpha`, slot 19 is
  posterior mean `sigmaSqAlpha`, slots 20-21 are `vle`/`vld`, slot 22 is
  marker component probabilities, slot 23 is posterior component counts, and
  slot 24 is LD-swap diagnostics. Multi-chain summary/chain slots follow the
  standard CSR convention.
- ST-BLR resemblance: high for base fields and BayesR component outputs; the
  annotation coefficient fields need standard names and metadata.

## R Wrapper and Helper Inventory

Public annotation-aware wrappers:

- `stblr_csr_prior_annot()` in `R/stblr-csr-prior-annot.R` calls
  `stblr_cpg_omp_csr_prior()`.
- `stblr_csr_learn_annot()` in `R/stblr-csr-learn-annot.R` calls
  `stblr_cpg_omp_csr_annot()`.
- `stblr_csr_group_annot()` in `R/stblr-csr-group-annot.R` calls
  `stblr_cpg_omp_csr_group_annot()`.
- `stblr_csr_sbayesrc_generic()` in `R/stblr-csr-sbayesrc.R` calls
  `stblr_cpg_omp_csr_sbayesrc()`.

Public helper and diagnostic functions:

- `make_sbayesrc_alpha_init()`
- `sbayesrc_annotation_pi()`
- `sbayesrc_annotation_gamma_mean()`
- `sbayesrc_marker_pi()`
- `sbayesrc_marker_gamma_mean()`
- `mtsim_annotation()`
- `summarize_annotation_signal()`

Internal helpers:

- `.stblr_get_nt_m_names()`, `.stblr_validate_stats()`,
  `.stblr_resolve_architecture()`, `.stblr_make_csr_variance_priors()`,
  `.stblr_init_marker_state()`, `.stblr_init_r_state()`,
  `.stblr_prepare_annotation_matrix()`,
  `.stblr_make_prior_from_annotations()`,
  `.stblr_prepare_group_index()`, `.stblr_expand_group_trait_matrix()`, and
  `.stblr_make_group_priors()` in `R/annotation-helpers.R`.
- `.format_csr_annot_fit()` in `R/stblr-csr-learn-annot.R`.
- `.format_csr_group_annot_fit()` in `R/stblr-csr-group-annot.R`.
- `format_sbayesrc_csr_fit()` in `R/sbayesrc-helpers.R`.

Existing tests and docs:

- `tests/testthat/test-sbayesrc-helpers.R` covers SBayesRC helper validation and
  formatting with synthetic native output.
- No focused tests currently cover `stblr_csr_prior_annot()`,
  `stblr_csr_learn_annot()`, or `stblr_csr_group_annot()` end-to-end.
- Man pages exist for all four public wrappers and exported SBayesRC helpers.
- `docs/annotation_informed_priors_stblr_mtblr.qmd` is broader design/research
  material, not a stable interface specification.

Naming and metadata inconsistencies:

- Wrapper names mix `prior_annot`, `learn_annot`, `group_annot`, and
  `sbayesrc_generic`.
- `fit$input$model` values are `"prior"`, `"annot"`, `"group"`, and
  `"sbayesrc"` rather than the current method/model/backend taxonomy.
- Annotation wrappers do not set `fit$input$method`, `fit$input$backend`,
  `fit$input$data_level`, `fit$input$annotation_model`, `fit$input$nchains`, or
  `fit$input$keep_chains`.
- Annotation matrix input is named `A`, while a future public interface should
  prefer `annotations` with `A` as a compatibility alias.
- Returned annotation fields use native names (`eta_pi`, `eta_vb`, `alpha`,
  `sigmaSqAlpha`, `group_pi`) but not a consistent `annotation_*` convention.

## Workflow Example Summary

`examples/workflows/annotation_based_models.R` expects a local qgg `Glist`, BED
files, and a local data directory. It constructs `cls`, computes sufficient
statistics with `bed_xtx_xty()`, streams sparse LD with
`sparseLD_stream_CSR()`, simulates overlapping annotations with
`mtsim_annotation()`, and verifies that the annotation matrix rows match
`Glist$rsidsLD[[chr]]`.

The workflow demonstrates four model variants:

- fixed marker-specific priors via `stblr_csr_prior_annot()`;
- learned annotation effects via `stblr_csr_learn_annot()`;
- mutually exclusive group priors via `stblr_csr_group_annot()`;
- SBayesRC-style mixture priors via `stblr_csr_sbayesrc_generic()`.

Downstream code uses `fit$bm`, `fit$dm`, variance traces, covariance/correlation
matrices, `summarize_annotation_signal()`, and SBayesRC helper outputs derived
from `fit_sbayesrc$alpha`, `fit_sbayesrc$input$A`, and `fit_sbayesrc$input$gamma`.

The workflow uses current function names and should run where the local data and
native build requirements are available. It is not suitable as a package example
because it has machine-specific paths and expensive BED/LD work. Useful pieces
for future tests are the compact fit-summary functions, marker/annotation row
alignment checks, and tiny synthetic calls to each wrapper.

## Model Taxonomy

| Proposed backend | Data level | Base method | Annotation type | Backend file | Current wrapper | Status |
| --- | --- | --- | --- | --- | --- | --- |
| `csr_prior_bayesc` | summary | bayesc | fixed marker-specific prior probabilities and/or variance multipliers | `src/st_cpg_omp_csr_prior.cpp` | `stblr_csr_prior_annot()` | plausible, documented, untested end-to-end |
| `csr_annot_bayesc` | summary | bayesc | dense marker x annotation matrix with learned inclusion/variance effects | `src/st_cpg_omp_csr_annot.cpp` | `stblr_csr_learn_annot()` | plausible, documented, untested end-to-end |
| `csr_group_bayesc` | summary | bayesc | one group index per marker | `src/st_cpg_omp_csr_group.cpp` | `stblr_csr_group_annot()` | plausible, documented, untested end-to-end |
| `csr_sbayesrc` | summary | sbayesrc / bayesr-like | dense marker x annotation matrix controlling mixture probabilities | `src/st_sbayesrc_omp_csr.cpp` | `stblr_csr_sbayesrc_generic()` | helper-tested, end-to-end wrapper untested |

## Proposed Public Interface

Use Option B first. This is now the active public annotation-aware interface:

```r
stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  ld_prefix = ld_prefix,
  method = "bayesC",
  annotations = annot,
  annotation_model = "group",
  nit = 1000,
  nburn = 100
)
```

and:

```r
stblr_csr_annot(
  stats = stats,
  Glist = Glist,
  ld_prefix = ld_prefix,
  method = "bayesR",
  annotations = annot,
  annotation_model = "sbayesrc",
  nit = 1000,
  nburn = 100
)
```

This is less ambiguous than immediately extending `stblr_csr()` because the
annotation models need different inputs: fixed marker priors, a group vector,
learned annotation effects, or SBayesRC mixture settings. After one stable
release of `stblr_csr_annot()`, `stblr_csr(..., annotations = ...)` can dispatch
to it for common cases.

Keep existing explicit wrappers as compatibility entry points:

- `stblr_csr_prior_annot()`
- `stblr_csr_learn_annot()`
- `stblr_csr_group_annot()`
- `stblr_csr_sbayesrc_generic()`

`stblr_csr_annot()` accepts `annotation_model = "prior"`, `"learned"`,
`"group"`, or `"sbayesrc"` plus common synonyms such as `"fixed_prior"`,
`"annot"`, `"groups"`, and `"SBayesRC"`. It dispatches to the existing wrappers
rather than calling native samplers directly.

## Proposed Internal Helper Naming

Add current-style helpers and make old names delegate where feasible:

- `.fit_stblr_csr_prior_bayesc()`
- `.fit_stblr_csr_annot_bayesc()`
- `.fit_stblr_csr_group_bayesc()`
- `.fit_stblr_csr_sbayesrc()`
- `.format_stblr_csr_prior_bayesc_fit()`
- `.format_stblr_csr_annot_bayesc_fit()`
- `.format_stblr_csr_group_bayesc_fit()`
- `.format_stblr_csr_sbayesrc_fit()`

Do not rename native symbols in this phase.

## Proposed Metadata Convention

All annotation-aware fits set:

- `fit$input$method`: `"bayesc"` or `"bayesr"`/`"sbayesrc"` as appropriate.
- `fit$input$model`: same base model family, such as `"bayesc"` or
  `"sbayesrc"`.
- `fit$input$backend`: one of `csr_prior_bayesc`, `csr_annot_bayesc`,
  `csr_group_bayesc`, or `csr_sbayesrc`.
- `fit$input$data_level`: `"summary"`.
- `fit$input$annotation_model`: `"prior"`, `"learned"`, `"group"`, or
  `"sbayesrc"`.
- `fit$input$annotations`: `TRUE`, indicating that the fit used an
  annotation-aware backend. Prepared matrices or groups remain in existing
  compatibility fields such as `fit$input$A` or `fit$input$group`.
- `fit$input$annotation_names`: annotation column names.
- `fit$input$nchains`: requested native chain count.
- `fit$input$keep_chains`: whether compact per-chain summaries were requested.
- `fit$input$updateE`, `fit$input$updatePi`, and model-specific controls such
  as `updateAlpha`, `learn_pi_annot`, `learn_vb_annot`, and `updateGroupVb`.

Keep `fit$input$A` and older fields for compatibility during the transition.

## Proposed Output Convention

All fits preserve:

- `fit$dm`
- `fit$bm`
- `fit$input`

BayesR/SBayesRC-like fits should also expose:

- `fit$comp_prob`
- `fit$dm_component_mean` when available or derivable

Annotation-specific fields use stable aliases when they can be constructed
reliably:

- `fit$annotation`: prepared annotation input or group metadata.
- `fit$annotation_summary`: compact annotation/group summary table.
- `fit$annotation_effects`: learned `eta_*` or `alpha` coefficients.
- `fit$annotation_pi`: group or marker/component probabilities when available.
- `fit$annotation_variance`: group variance multipliers or `sigmaSqAlpha`.
- `fit$annotation_prior`: fixed marker priors used by `csr_prior_bayesc`.
- `fit$annotation_enrichment`: optional derived summaries from helper code.

Existing native names (`eta_pi`, `eta_vb`, `alpha`, `sigmaSqAlpha`, `group_pi`,
`group_vb_multiplier`) should remain available initially.

## Gap Analysis

| Backend | nchains | chain summaries | keep_chains | updateE | LD-swap/MH | scheduled | standard dm/bm | standard metadata | tests | docs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `csr_prior_bayesc` | yes | yes | yes | yes | yes | no | yes | yes | focused chain and LD-swap tests | man page |
| `csr_annot_bayesc` | yes | yes | yes | yes | yes | no | yes | yes | focused chain and LD-swap tests | man page |
| `csr_group_bayesc` | yes | yes | yes | yes | yes | no | yes | yes | focused chain and LD-swap tests | man page |
| `csr_sbayesrc` | yes | yes | yes | yes | yes | no | yes | yes | focused chain and LD-swap tests | man page |

## Compatibility Strategy

1. Keep existing exported functions working.
2. Add `stblr_csr_annot()` as the clean public annotation entry point.
3. Make old wrappers call new internal helpers where feasible.
4. Do not remove old names or native symbols.
5. Keep old output fields while adding standard aliases.
6. Add tests that old and new wrappers return structurally equivalent fits on
   tiny fixtures.
7. Avoid C++ renaming or sampler math changes in this alignment phase.

## Testing Plan

Use tiny synthetic fixtures rather than expensive real BED workflows.

- Formatter tests for each native return layout.
- Smoke tests for each public wrapper when native compilation is available.
- Metadata tests for `method`, `model`, `backend`, `data_level`,
  `annotation_model`, `nchains`, `updateE`, and `updatePi`.
- Output tests for `dm`, `bm`, row/column names, annotation aliases, and
  SBayesRC `comp_prob`.
- Compatibility tests comparing old wrappers with `stblr_csr_annot()` dispatch.
- Workflow-level tests should only check row alignment and compact summary
  helpers, not stream large sparse LD.

## Recommended Implementation Sequence

1. Standardize formatter output and metadata for `csr_prior_bayesc`, because it
   is the simplest BayesC-like annotation backend.
2. Add `stblr_csr_annot()` dispatch for `annotation_model = "prior"` and keep
   `stblr_csr_prior_annot()` as a compatibility wrapper.
3. Add tests for `fit$dm`, `fit$bm`, metadata, and fixed-prior annotation
   summaries.
4. Align `csr_group_bayesc` and `csr_annot_bayesc` to the same metadata and
   annotation alias conventions.
5. Align `csr_sbayesrc` output to the same conventions, including
   `comp_prob`, `annotation_effects`, `annotation_variance`, and LD-swap
   diagnostics.
6. Add compatibility tests for old and new wrapper structural equivalence.
7. Only later consider scheduled updates for annotation-aware models.
