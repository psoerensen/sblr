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

## Non-Migrated Backends

The old positional formatter remains active for non-migrated backends,
including scheduled CSR BayesC, BED backends, finemap paths, and
individual-level scheduled backends.

Wrappers should explicitly detect:

```r
raw$schema$class == "stblr_raw"
as.integer(raw$schema$version) == 1L
```

before routing through `.format_stblr_raw_v1()`.

## Planned Migration Order

The remaining follow-up order is scheduled and BED backends after the
summary-statistic CSR schemas have stabilized.
