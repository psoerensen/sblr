# ST-BLR Raw Schema v1

## Purpose

The named raw schema replaces fragile positional native return slots with a
stable internal list contract. Different ST-BLR CSR models currently reuse the
same slot positions for different concepts, which makes formatter changes risky
as new diagnostics and model families are added.

Raw schema v1 is an internal development interface. Public `fit` object field
names remain unchanged.

## Top-Level Namespaces

Version 1 uses these top-level names:

```r
schema
meta
marker
trace
variance
pi
diagnostics
chains
prior
group
annotation
component
selection
```

Unused namespaces may be empty named lists or `NULL` for a backend that does
not use them.

## Canonical Dimensions

Marker-level matrices use marker rows and trait columns:

```text
m x nt
```

Trace matrices use MCMC iteration rows and trait columns:

```text
n_trace x nt
```

For ordinary CSR BayesC, `marker$dm` is `P(d = 1)` and `trace$pis` is the
sampled active-marker probability `pi1` for every saved iteration.

## Migrated Backends

Phase 1 migrated ordinary summary-statistic CSR BayesC:

```text
src/st_cpg_omp_csr.cpp
```

This backend returns `schema`, `meta`, `marker`, `trace`, `variance`, `pi`,
`diagnostics`, `chains`, and `selection`. The `prior`, `group`, `annotation`,
and `component` namespaces are present as empty lists.

Phase 2 migrates summary-statistic CSR BayesR:

```text
src/st_cpg_omp_csr_bayesr.cpp
```

CSR BayesR uses the `component` namespace. `raw$component$prob` is a list of
length `nt`, with one `m x K` marker-by-component posterior probability matrix
per trait. The null component is always named `component_0`, and the formatted
`fit$dm` is derived as `1 - P(component_0)`.

For BayesR, `raw$trace$pis` and formatted `fit$pis` are the total active-marker
probability trace, `1 - pi_component_0`. Final and posterior mean mixture
probabilities are stored in `raw$pi$final` and `raw$pi$mean` as `nt x K`
matrices with component names.

The R formatter `.format_stblr_raw_v1()` consumes only this named schema and
maps it back to the existing user-facing fit fields such as `bm`, `dm`, `vbs`,
`vgs`, `ves`, `vle`, `vld`, `pi`, `pim`, `pis`, `chains`, `ld_swap`, and
sampled `selection_s` summaries. For BayesR it also maps `component$prob` to
`fit$comp_prob` and `component$dm_component_mean` to
`fit$dm_component_mean`.

Phase 3 migrates summary-statistic CSR SBayesRC:

```text
src/st_sbayesrc_omp_csr.cpp
```

CSR SBayesRC uses both the `component` and `annotation` namespaces.
`raw$component$prob` is a list of length `nt`, with one `m x K`
marker-by-component posterior probability matrix per trait. The null component
is always named `gamma_0.00`, and formatted `fit$dm` is derived as
`1 - P(gamma_0.00)`.

For SBayesRC, `raw$trace$pis` and formatted `fit$pis` are the total
active-marker probability trace: marker-averaged `1 - P(gamma_0.00)`.
Final and posterior mean component probabilities are stored in `raw$pi$final`
and `raw$pi$mean` as `nt x K` matrices with gamma component names.

`raw$annotation` stores SBayesRC annotation parameters, including `alpha` as a
trait list of `nAnno x (K - 1)` matrices and `sigmaSqAlpha` as
`(K - 1) x nt` matrices. The R formatter maps these back to existing
SBayesRC fit fields such as `alpha`, `sigmaSqAlpha`, `annotation_summary`,
`annotation_pi`, and `annotation_effects`.

Phase 4 migrates the BayesC-like annotation/prior summary-statistic CSR
backends:

```text
src/st_cpg_omp_csr_prior.cpp  -> csr_prior_bayesc
src/st_cpg_omp_csr_group.cpp  -> csr_group_bayesc
src/st_cpg_omp_csr_annot.cpp  -> csr_annot_bayesc
```

All three use BayesC-style state semantics: `raw$marker$dm` is `P(d = 1)`,
and `raw$pi$final` / `raw$pi$mean` are `nt x 2` matrices named `pi0` and
`pi1`.

Marker-prior BayesC uses `raw$prior`. It stores the resolved marker-specific
prior probabilities and marker-effect variance multipliers as `m x nt`
matrices, including fixed inputs when the priors are not sampled.

Group-prior BayesC uses `raw$group`. `raw$group$pi_mean`,
`raw$group$vb_multiplier_mean`, and `raw$group$n_included_mean` use
`ngroup x nt` layout, with `raw$group$group_index` retaining the native
0-based marker group index.

Learned-annotation BayesC uses `raw$annotation`. It stores learned `eta_pi`
and `eta_vb` summaries as `nAnno x nt` matrices. SBayesRC-specific
`alpha`/`sigmaSqAlpha` fields are not invented for this backend.

For these three backends, `raw$trace$pis` is the existing tracked
active-marker probability trace. It is the sampled global `pi1` for
marker-prior and learned-annotation BayesC, and the marker-weighted current
group probability for group-prior BayesC. Unused namespaces remain empty.

Phase 5 completes the migration for the remaining backends: scheduled CSR
BayesC, BED BayesC/BayesR, and the individual-level (non-CSR) backends:

```text
src/st_cpg_omp_csr_scheduled.cpp
src/st_cpg_omp_individual.cpp
src/st_cpg_omp_individual_scheduled.cpp
src/st_cpg_omp_individual_scheduled_chains.cpp
src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp
```

These backends reuse the same BayesC/BayesR namespace conventions described
above for their respective model family (CSR BayesC or BayesR semantics for
`raw$marker`, `raw$pi`, and, for BayesR, `raw$component`).

## All Active Backends Are Migrated

Every active, R-reachable native backend now returns `stblr_raw_v1` output
(schema `class = "stblr_raw"`, `version = 1L`). There are no remaining
non-migrated active backends.

R-side wrappers detect the schema explicitly:

```r
raw$schema$class == "stblr_raw"
as.integer(raw$schema$version) == 1L
```

(via the internal helper `.is_stblr_raw_v1()`) before routing through
`.format_stblr_raw_v1()`.

## Legacy Positional Backend Output Is Unsupported

Because every active backend now emits `stblr_raw_v1`, the positional
backend-output path is no longer a supported fallback at the wrapper level.
Every R-facing model-fitting wrapper (`stblr_csr()`, `stblr_csr_annot()`,
`stblr_bed()`, `stblr_csr_bayesr()`, the BED marker helpers, and the finemap
local re-fit helper) stops with a clear error via
`.stblr_stop_unsupported_raw_output()` if a backend ever returns something
that is not a valid `stblr_raw_v1` object, instead of silently reformatting
positional output.

The underlying positional formatters (`.format_stblr_fit()`,
`.format_stblr_bayesr_fit()`, `.format_stblr_csr_bayesr_fit()`,
`.format_csr_annot_fit()`, `.format_csr_group_annot_fit()`,
`format_sbayesrc_csr_fit()`) are intentionally retained as internal
compatibility code. They are exercised directly by unit tests that construct
legacy-shaped positional lists (see `test-bayesr-csr-backend.R` and
`test-sbayesrc-helpers.R`), but no production code path calls them anymore.

## Stable Formatted Fit Contract

`stblr_raw_v1` is an internal development interface and may change between
backend phases. The formatted `fit` object returned to users is the stable
contract: public field names (`bm`, `dm`, `comp_prob`, `chains`, `ld_swap`,
`dm_sd`/`dm_min`/`dm_max`, `bm_sd`/`bm_min`/`bm_max`, `selection_s*`, etc.),
their dimensions, and `fit$input` metadata do not change when the raw schema
changes underneath them. [`check_stblr_consistency()`](../../R/check-stblr-backend-consistency.R)
is the executable description of this contract and should keep passing for
every backend regardless of raw-schema revisions.

## Present-but-`NULL` Diagnostic Fields

Optional diagnostic fields are always present as named elements on the
formatted `fit` object, set to `NULL` rather than omitted, when a backend or
configuration does not produce them. This lets callers use `fit$field` rather
than `"field" %in% names(fit)` to test for availability. Examples:

- `fit$ld_swap` and `fit$ld_swap_chains` (`.stblr_ensure_ld_swap_fields()`) are
  always present, `NULL` unless `updateLDswap = TRUE` (or chain-level swap
  diagnostics were kept).
- `fit$dm_sd`/`fit$dm_min`/`fit$dm_max`/`fit$bm_sd`/`fit$bm_min`/`fit$bm_max`
  (`.stblr_ensure_chain_summary_fields()`) are always present. They are
  `NULL` when the backend has no chain-summary data for them (e.g.
  `nchains > 1` but the native backend did not return per-chain summaries),
  except for the single-chain convention below.

## Single-Chain Degenerate Chain-Summary Convention

When a fit uses exactly one chain (`meta$nchains == 1`) and the native raw
object does not itself supply `dm_sd`/`dm_min`/`dm_max`/`bm_sd`/`bm_min`/
`bm_max`, `.stblr_ensure_chain_summary_fields()` synthesizes them from `dm`/`bm`
instead of leaving them `NULL`:

- `dm_sd` and `bm_sd` are set to zero (same dimensions as `dm`/`bm`).
- `dm_min`/`dm_max` and `bm_min`/`bm_max` are set equal to `dm`/`bm`.

This is intentional: with a single chain there is no cross-chain variation,
so the degenerate summary (zero spread, min == max == the point estimate) is
the mathematically correct value rather than a missing one. It also means
`check_stblr_consistency()`'s single-chain checks
(`chain_summary.single.dm_sd_zero`, `chain_summary.single.dm_equal`, etc.)
can assert on these fields whenever a backend chooses to populate them for a
single chain.
