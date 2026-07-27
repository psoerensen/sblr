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

## Phase 17E production multivariate names

`MtDefaultDataView`, `MtDefaultModelSpec`,
`MtDefaultCovariancePriorView`, `MtDefaultExecutionSpec`,
`MtDefaultInitialState`, and `MtDefaultCoreResult` are production vocabulary for
the authoritative default legacy route. `run_mt_default_core()` is its sole
callable numerical core. These names do not denote a generic MT framework and
do not apply to `mtblr_cpg*`, `mtblr_eigen`, or `mtblr_hybrid`. “Result” means a
binding-neutral numerical result consumed by the still-inline legacy
20-position finalizer, not `stblr_raw_v1`.

## Phase 17F shared ST/MT naming plan

### Phase 17G active multivariate names

`mtblr` alone names the supported public dense MT BayesC sampler.
`mtblr_eigen` and `mtblr_cpg_omp_csr` name unsupported native-only research
evidence and must not appear as `sblr()` algorithm values. The latter's local
`LDCSR` is explicitly noncanonical. Retired names `mtblr_cpg`,
`mtblr_cpg_arma`, `mtblr_cpg_omp`, and `mtblr_hybrid` remain only in historical
documentation and Git history.

Legacy public names remain unchanged. New typed names distinguish final values,
posterior means, and traces and avoid MT-only synonyms.

| Scientific concept | Current scalar name | Current MT legacy name | Phase 17F typed name | Planned canonical shared name | Compatibility decision |
|---|---|---|---|---|---|
| sample size | `sample_size`, `n` | `n` | data-view `n` | `sample_size` | translate only at legacy boundary |
| marker IDs | marker names/input metadata | marker dimnames | none yet | `marker_ids` | shared validation metadata |
| trait IDs | trait names | trait dimnames | none yet | `trait_ids` | shared validation metadata |
| summary cross-products | `wy` | `wy` | `wy` | `wy` | reuse |
| marker diagonal | `diagonal`, `ww` | `ww` | `ww` | `marker_diagonal` with `ww` compatibility | clarify without changing public field |
| sparse LD | CSR row/column/value buffers | experimental CSR prefix | none | `SparseLdCsrView` | reuse canonical scalar representation exactly |
| dense LD | dense/operator-specific LD | nested `XXvalues/XXindices` summaries | same borrowed view | `ld_operator` metadata | do not claim representation equivalence |
| initial marker effects | `initial_effects`, `b` | `b` | initial-state `b` | `initial_effects` | typed canonical name later |
| posterior marker means | `bm` | `bm` | `bm` | `bm` | reuse |
| inclusion probabilities | `dm` | `dm` | `dm` | `dm` | reuse |
| final marker effects | `b` | `b` | `b` | `final_effects` | preserve public `b` |
| final states | `d` | `d` | `d` | `final_states` | preserve public `d` |
| marker variance/covariance | `vb`, `covb` | `vb`, `covb` | `vb`, `covb` | `marker_covariance` | scalar variance is a 1x1 specialization |
| genetic variance/covariance | `vg`, `covg` | `vg`, `covg` | `vg`, `covg` | `genetic_covariance` | same specialization rule |
| residual variance/covariance | `ve`, `cove` | `ve`, `cove` | `ve`, `cove` | `residual_covariance` | do not equate with overlap policy |
| mixture/pattern probabilities | `pi`, `pim`, `pis` | `pi`, `pim` | `pi_final`, `pi_mean` | `pi_final`, `pi_mean`, `pi_trace` | distinguish final/mean/trace |
| iteration count | `iterations`, `nit` | `nit` | `nit` | `iterations` | preserve public argument |
| burn-in | `burnin`, `nburn` | `nburn` | `nburn` | `burnin` | preserve public argument |
| thinning | `thinning`, `nthin` | `nthin` | `nthin` | `thinning` | preserve public argument |
| input metadata | `input` | legacy call metadata | none in numerical types | `input` | binding-owned |
| diagnostics | named optional fields | legacy optional/omitted | none in final result | named diagnostics | binding-owned, present policy model-specific |

### Future MT operator, alignment, and study metadata

Canonical MT CSR must reuse the scalar CSR representation and expose one view per trait/study.
Shared row offsets and column indices are permitted only when identical; LD
values, `ww`, `wy`, and sample size remain trait-specific. Fully independent CSR structures
are required for different ancestries, cohorts,
populations, LD panels, or sparsity thresholds; no artificial union pattern is
allowed. Canonical MT block-eigen must likewise reuse scalar conventions and
allow trait-specific blocks, eigenvectors, eigenvalues, retained ranks,
tolerances, and panel metadata.

Marker row `i` must denote the same canonical marker ID and effect-allele orientation
in every summary and operator. R validation owns marker matching,
trait/study order, allele orientation, duplicates, missing-marker policy, and
panel metadata. The first canonical MT CSR route should use explicit marker
intersection unless a union-with-mask likelihood is specified and tested.
Future metadata must distinguish trait, study, population, ancestry, LD
reference, sample size, marker set, and sample-overlap policy. Independent,
known-overlap, and unknown-overlap studies are distinct; residual covariance
alone must not silently stand for GWAS sample overlap.
## Shared sparse-LD names (Phase 17H)

Use `SparseLdCsrStorage` for the model-neutral owner and `SparseLdCsrView` for
one immutable operator view. `row_ptr`, `column_index`, and `offdiag_xij` are
view vocabulary; the owner's legacy `ptr`, `idx`, and `xij` members remain a
temporary source-compatibility detail. A future MT bundle contains one view per
trait/study and must not introduce an MT-specific CSR representation.
## Phase 17I MT CSR vocabulary

`MtSparseLdBundleView` means an ordered vector of canonical trait/study
operators on one marker domain. `MtCsrDataView` adds aligned `wy`, `yy`, and
`n`. Neither name implies shared LD; sharing is a property of borrowed pointers.
Public MT CSR naming uses `marker_id`, `effect_allele`, `other_allele`,
`allele_frequency`, `trait_id`, `study_id`, `population`, `ancestry`, and
`ld_reference`. Legacy numerical field names remain compatible in `mtblr_fit`.

## Block-filtered names (Phase 17K)

`BlockEigenStorage` is the canonical owner and `BlockEigenView` its immutable
borrower. “Block eigen” describes construction/filtering; runtime data are
`BlockEigenBlockStorage::upper_triangle` float-packed filtered dense blocks.
Do not use names implying retained eigenvectors or low-rank runtime factors.
# Phase 17L naming

“MT block-eigen” names the internal canonical route using `BlockEigenView`; it does not refer to legacy `mtblr_eigen()`. The new native maintenance symbol is `mtblr_block_eigen_internal()` and is not a public API.

# Phase 17M naming

The public function is `mtblr_block_eigen()` and its raw backend identifier is
`mt_block_eigen_bayesc`. `mtblr_eigen()` remains an unsupported research
symbol and is not an implementation component of the public route.

# Phase 17O internal naming

Phase 17O activates internal backend identifier `mt_bed_bayesc` with data
level `individual`. The native execution and inspection symbols are
`mtblr_bed_internal()` and `mtblr_bed_marker_contract_internal()`.
Phase 17P exports `mtblr_bed()` as the public adapter for this backend. The
native symbols remain internal and `mtblr_eigen()` remains unrelated research
code.

`mtblr_bed_chains_internal` is reserved as the proposed Phase 17R internal
name. It is not registered or callable in Phase 17Q; backend naming remains
`mt_bed_bayesc` and raw schema version remains 1.

Phase 17R registers `mtblr_bed_chains_internal` as an internal maintenance
route. Both serial and chains routes retain backend `mt_bed_bayesc`, data level
`individual`, and `mtblr_raw` version 1; no public export is added.

Phase 17S retains `mt_bed_bayesc` for both public single- and multichain calls.
Execution topology is represented by version-1 chain metadata, not a new name.

Core diagnostic quantities use `vbs[trait]`, `vgs[trait]`, `ves[trait]`,
`vle[trait]`, and `vld[trait]`; future extended diagnostics use trait-pair covariance names, probabilities, and
`b/d[marker_id,trait]`, never exposed native zero-based indices.
## Phase 17U convergence names

Internal Tier 1 quantities use the family-neutral `vbs[trait]`, `vgs[trait]`,
`ves[trait]`, `vle[trait]`, and `vld[trait]` names. The binding symbol is
`mtblr_bed_convergence_trace_internal` and remains unexported.

Phase 17V exposes quantity names at `fit$convergence$summary` and reserves
`fit$convergence_traces` for optional post-burn Tier 1 arrays.

Phase 17W reserves future groups `B_cov`, `G_cov`, `E_cov`, `pi_mass`,
`pi_pattern`, `marker_b`, and `marker_d`, preserving fitted trait/model/final
marker order. These names are contractual and not yet public output.
## Phase 18 canonical names

Public models use only `bayesc`, `bayesr`, `bayesrc`, and `sbayesrc`; `cpg`
remains an internal historical kernel term. Public operators use `csr`,
`block_eigen`, and `packed_bed`. Deprecated case and wrapper aliases are not
retained.
