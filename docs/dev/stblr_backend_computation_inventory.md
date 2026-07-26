# ST-BLR Backend Computation and Return Inventory

## Experimental packed-BED BayesC routes

`stblr_bed_marker()` retains two explicitly experimental implementations outside
the canonical `stblr_bed(method = "bayesc")` path. The scheduled single-chain
export `stblr_cpg_omp_bed_marker_scheduled()` is retained as a deterministic
Phase 11B reference implementation. The sparse export
`stblr_cpg_omp_bed_marker_sparse()` retains a distinct `null_update_prob`
marker-scheduling experiment. Both decode fit-local packed BED storage, are
noncanonical, and are unsupported for general use; neither is a fallback from
the canonical public route.

This document records what the supported BLR backends compute, accumulate,
return from native code, and expose in formatted R fit objects. It is an audit
artifact for backend alignment; it does not change sampler definitions.

## Scope

Audited public backends:

| backend | method/model | data level | C++ file | R interface | R formatter/helper |
| --- | --- | --- | --- | --- | --- |
| `csr_bayesc` | BayesC | summary | `src/st_cpg_omp_csr.cpp` | `stblr_csr(method = "bayesC")` | `.format_stblr_fit()` |
| `csr_scheduled_bayesc` | BayesC | summary | `src/st_cpg_omp_csr_scheduled.cpp` | `stblr_csr(method = "bayesC", scheduled = TRUE)` | `.format_stblr_fit()` |
| `csr_bayesr` | BayesR | summary | `src/st_cpg_omp_csr_bayesr.cpp` | `stblr_csr(method = "bayesR")`, `stblr_csr_bayesr()` | `.format_stblr_csr_bayesr_fit()` |
| `bed_bayesc` | BayesC | individual | `src/st_cpg_omp_individual_scheduled_chains.cpp` | `stblr_bed(method = "bayesC")` | `.format_stblr_fit()` |
| `bed_bayesr` | BayesR | individual | `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp` | `stblr_bed(method = "bayesR")` | `.format_stblr_bayesr_fit()` |
| `csr_prior_bayesc` | BayesC fixed annotation prior | summary | `src/st_cpg_omp_csr_prior.cpp` | `stblr_csr_annot(annotation_model = "prior")` | `.format_stblr_fit()` plus prior aliases |
| `csr_annot_bayesc` | BayesC learned annotation effects | summary | `src/st_cpg_omp_csr_annot.cpp` | `stblr_csr_annot(annotation_model = "learned")` | `.format_csr_annot_fit()` |
| `csr_group_bayesc` | BayesC group priors | summary | `src/st_cpg_omp_csr_group.cpp` | `stblr_csr_annot(annotation_model = "group")` | `.format_csr_group_annot_fit()` |
| `csr_sbayesrc` | SBayesRC-style BayesR | summary | `src/st_sbayesrc_omp_csr.cpp` | `stblr_csr_annot(annotation_model = "sbayesrc")` | `format_sbayesrc_csr_fit()` |

`src/st_sbayesrc_omp_csr_annot.cpp` exists in the source tree but is not the
current audited public SBayesRC path.

The likelihood-independent SBayesRC probit stick-breaking probability and
annotation-coefficient update utilities live in
`src/st_bayesrc_annotation_prior.h`. CSR SBayesRC and individual-level BED
BayesRC share these utilities.

Sequence 2 added the `bed_bayesrc` native backend. Sequence 3 routes it through
public `stblr_bed(method = "bayesrc")`. It combines the
existing packed-BED BayesR likelihood calculations with the shared probit
stick-breaking annotation prior and uses structurally explicit full marker
sweeps. Adaptive scheduling,
`selection_s`, LD-swap, and annotation-dependent effect variances are disabled.

Shared packed-BED storage, decoding, residual, variance, and CPO utilities for
ordinary BED BayesR and internal BED BayesRC are defined inline in
`src/st_bed_bayesr_common.h`. Model-specific marker priors, scheduling, chain
runners, and raw-output aggregation remain in their respective translation
units.

The internal BED BayesRC checkpoint has focused regression coverage for fixed
BayesR-prior reduction, learned annotation priors, component/prior identities,
traits and chains, deterministic seeds and thread counts, variance/CPO
diagnostics, native validation, and schema structure. Public readiness still
depends on successful compiled-platform execution of those native tests.

## Backend Feature Matrix

Values are `yes`, `no`, `partial`, or `n/a`.

| backend | return type | nchains | chain_seeds | keep_chains | updateE | LD-swap/MH | scheduled | dm | bm | dm/bm sd/min/max | vbs/vgs/ves | vle/vld | covb | pi/pim | comp_prob | dm_component_mean | alpha/sigmaSqAlpha | ncomp | annotation aliases | ld_swap | ld_swap_chains | chains |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `csr_bayesc` | legacy slots | yes | yes | yes | yes | yes | no | yes | yes | yes | yes | yes | yes | yes | n/a | n/a | n/a | n/a | n/a | yes | yes | yes |
| `csr_scheduled_bayesc` | legacy slots | yes | yes | no | yes | no | yes | yes | yes | yes | yes | yes | yes | yes | n/a | n/a | n/a | n/a | n/a | no | no | no |
| `csr_bayesr` | `Rcpp::List` | yes | yes | yes | yes | yes | no | yes | yes | yes | yes | yes | yes | yes | yes | yes | n/a | yes | n/a | yes | yes | yes |
| `bed_bayesc` | legacy slots | yes | no | no | yes | no | yes | yes | yes | yes | yes | yes | yes | yes | n/a | n/a | n/a | n/a | n/a | no | no | no |
| `bed_bayesr` | legacy slots | yes | no | no | yes | no | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | n/a | n/a | n/a | no | no | no |
| `csr_prior_bayesc` | legacy slots | yes | yes | yes | yes | yes | no | yes | yes | yes | yes | yes | yes | yes | n/a | n/a | n/a | n/a | yes | yes | yes | yes |
| `csr_annot_bayesc` | legacy slots | yes | yes | yes | yes | yes | no | yes | yes | yes | yes | yes | yes | yes | n/a | n/a | n/a | n/a | yes | yes | yes | yes |
| `csr_group_bayesc` | legacy slots | yes | yes | yes | yes | yes | no | yes | yes | yes | yes | yes | yes | yes | n/a | n/a | n/a | n/a | yes | yes | yes | yes |
| `csr_sbayesrc` | legacy slots | yes | yes | yes | yes | yes | no | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes |

## Common Quantity Definitions

`dm` is the posterior marker inclusion probability. For BayesR-like backends it
is `P(component > 0)` and is formatted from `comp_prob`.

`bm` is the posterior mean marker effect.

`vbs`, `vgs`, and `ves` are post-burn-in/thinned trace matrices by iteration
and trait for marker-effect variance, genetic variance, and residual variance.

`vle` is the linkage-equilibrium marker-effect component trace. In CSR and BED
code this is accumulated from the current sampled marker effects and inclusion
state as the diagonal marker contribution before the sparse-LD cross terms.

`vld` is the LD-induced genetic variance contribution trace and is accumulated
as `vgs - vle`. The formatted R object stores both as iteration-by-trait
matrices with row names `Iter1`, `Iter2`, ... and trait column names.

`covb`, `covg`, and `cove` are posterior mean trait covariance matrices; `vb`,
`vg`, and `ve` are the final covariance state matrices.

`pi` is the final inclusion or component-probability state. `pim` is the
posterior mean inclusion or component-probability state. BayesR-like formatters
also expose `final_pi` and `mean_pi` aliases where supported.

`dm_sd`, `dm_min`, `dm_max`, `bm_sd`, `bm_min`, and `bm_max` are across-chain
summaries of marker-level posterior means. They are returned for multi-chain
native backends that compute chain summaries; some BED and scheduled backends
compute summaries but do not return compact per-chain objects.

For the standard CSR workflow, `make_summary_stats(..., scale = TRUE)` constructs
summary statistics from standardized genotype columns,
`(dosage - 2p) / sqrt(2p(1 - p))`, and `make_sparse_ld()` stores correlations
between the same standardized columns. The native CSR samplers update `b`
against `wy`, `ww`, and sparse LD cross-products from that scale. Therefore
CSR `bm` is a standardized-genotype-scale marker effect in the supported
package path. `ww_j = x_j' x_j`; with standardized columns, `ww_j / n` is
approximately one, not `h_j = 2p_j(1 - p_j)`.

The diagonal LE trace `vle = sum_j ww_j * b_j^2 / n` is consequently
consistent with standardized-genotype-scale effects in the standard CSR path.
If a BayesS-style prior is defined on allele-scale effects with heterozygosity
exponent `S`, then the corresponding standardized-effect prior for current CSR
samplers is `b_j ~ N(0, v_m h_j^(S + 1))`. Fixed global
`selection_s` implements this fixed scaling for `csr_bayesc`, `csr_bayesr`,
and `csr_sbayesrc`. Sampled trait-specific `selection_s` is currently
implemented for annotation-unaware unscheduled `csr_bayesc`
(`stblr_csr(method = "bayesC", scheduled = FALSE)`), `csr_bayesr`
(`stblr_csr(method = "bayesR")`), and annotation-aware `csr_sbayesrc`
(`stblr_csr_annot(annotation_model = "sbayesrc")`). Its default bounded
uniform prior is `c(-3, 2)`, and its default random-walk MH proposal SD is
0.35. For CSR BayesC, included-marker prior variances are
`vb * h_j^(selection_s + 1)`. For CSR BayesR, non-null component prior
variances are `vb * mixture_var_m * h_j^(selection_s + 1)`, with component 0
remaining the point-mass null. For CSR SBayesRC, annotations continue to
control marker-specific component probabilities, while fixed `selection_s`
scales non-null effect prior variances as
`vb * gamma_m * h_j^(selection_s + 1)`. The scheduled CSR BayesC backend,
BayesC-like annotation-aware CSR backends (`csr_prior_bayesc`,
`csr_annot_bayesc`, and `csr_group_bayesc`), and BED backends do not currently
support sampler-level `selection_s`.

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

User-facing model capability map:

| Model/API | Fixed `selection_s` | Sampled `selection_s` | Annotation-aware | `comp_prob` | Annotation outputs | LD-swap |
| --- | --- | --- | --- | --- | --- | --- |
| `stblr_csr(method = "bayesC")` | yes | yes | no | no | no | yes |
| `stblr_csr(method = "bayesR")` | yes | yes | no | yes | no | yes |
| `stblr_csr_annot(annotation_model = "prior")` | no | no | yes | no | yes | yes |
| `stblr_csr_annot(annotation_model = "learned")` | no | no | yes | no | yes | yes |
| `stblr_csr_annot(annotation_model = "group")` | no | no | yes | no | yes | yes |
| `stblr_csr_annot(annotation_model = "sbayesrc")` | yes | yes | yes | yes | yes | yes |

`comp_prob` follows backend-specific null component names. CSR BayesR uses
`component_0`, so `dm = 1 - P(component_0)`.
CSR SBayesRC uses `gamma_0.00`, so `dm = 1 - P(gamma_0.00)`.
Fine-mapping diagnostics are represented by
marker PIPs in `dm` and optional LD-swap summaries in `ld_swap` when
`updateLDswap = TRUE`. Credible sets are produced by R helper functions from
PIPs and LD rather than returned as a separate native sampler object.

Sampled CSR BayesC uses a bounded uniform prior and the Metropolis-Hastings
kernel

```text
log p(S_t | b_t, d_t, vb_t)
  = log p(S_t)
    - 0.5 * sum_{j: d_jt = 1} [
        log(q_j(S_t)) + b_jt^2 / (vb_t * q_j(S_t))
      ]
```

with `q_j(S_t) = h_j^(S_t + 1)`. Native code receives aligned `log(h_j)` and
computes dynamic marker scales as `exp((S_t + 1) * log(h_j))`.

Sampled CSR BayesR uses the analogous trait-chain-specific MH update with
active non-null components:

```text
log p(S_t | b_t, gamma_t, vb_t)
  = log p(S_t)
    - 0.5 * sum_{j: gamma_jt > 0} [
        log(q_j(S_t)) + b_jt^2 / (vb_t * gamma_jt * q_j(S_t))
      ]
```

where `gamma_jt` is the current BayesR non-null component variance multiplier
and `q_j(S_t) = h_j^(S_t + 1)`.

Sampled CSR SBayesRC uses the same active-component kernel with the current
SBayesRC component labels:

```text
log p(S_t | b_t, gamma_t, vb_t)
  = log p(S_t)
    - 0.5 * sum_{j: gamma_jt > 0} [
        log(q_j(S_t)) + b_jt^2 / (vb_t * gamma_jt * q_j(S_t))
      ]
```

Here `gamma_jt` is the current non-null SBayesRC component variance multiplier
and `q_j(S_t) = h_j^(S_t + 1)`. Annotations affect component probabilities and
alpha updates; `selection_s` affects marker-specific effect-size prior
variance.

### Fixed `selection_s` Performance Path

The default `selection_s = NULL` path should be treated as the baseline CSR
BayesC/BayesR/SBayesRC implementation. Native code should not allocate an all-ones prior
scale vector for this path, and hot marker-update loops should not read
`prior_scale[j]`, multiply by a unit marker scale, or branch on
`selection_s`.

Current CSR BayesC, CSR BayesR, and CSR SBayesRC native code uses separate unscaled helper
paths when `selection_s_prior_scale` is absent or empty. The scaled helper paths
are used only when a fixed `selection_s` is supplied. This intentionally
duplicates a small amount of marker-update, variance-update, and LD-swap code
so the default path stays close to the previous unscaled C++ code.
Sampled `selection_s` for CSR BayesC, CSR BayesR, and CSR SBayesRC is a third path: the
sampler recomputes a dynamic marker-scale vector once per trait-chain
iteration and reuses the scaled beta, `vb`, and LD-swap helpers.

Audit conclusion:

- `selection_s = NULL`: not expected to be meaningfully slower than the
  previous C++ implementation from fixed-`selection_s` support, because the
  marker-specific prior-scale lookup and multiply are not on the default hot
  path.
- `selection_s = -1`: expected to give identical fitted values to `NULL`, but
  may be slower because it explicitly passes and validates a unit scale vector.
- `selection_s = 0` or `selection_s = -0.5`: expected to be modestly slower
  than `NULL`, because marker-specific prior variances are used.
- CSR BayesR is the backend most likely to show overhead for fixed
  `selection_s`, because the component-posterior loop evaluates multiple
  mixture components per marker.
- CSR SBayesRC has the same component-loop consideration as CSR BayesR and
  additionally combines fixed MAF scaling with annotation-dependent component
  probabilities.

Manual timing helper for ad hoc development checks, assuming `stats` and
`Glist` are already defined in the current R session:

```r
bench_csr_selection_s <- function(method = c("bayesC", "bayesR"),
                                  selection_s = NULL,
                                  nrep = 3,
                                  seed = 10) {
  method <- match.arg(method)
  out <- numeric(nrep)

  for (i in seq_len(nrep)) {
    gc()
    out[i] <- system.time({
      stblr_csr(
        stats = stats,
        Glist = Glist,
        method = method,
        selection_s = selection_s,
        seed = seed
      )
    })[["elapsed"]]
  }

  data.frame(
    method = method,
    selection_s = if (is.null(selection_s)) NA_real_ else selection_s,
    mean_seconds = mean(out),
    sd_seconds = sd(out),
    min_seconds = min(out),
    max_seconds = max(out),
    nrep = nrep
  )
}

rbind(
  bench_csr_selection_s("bayesC", NULL),
  bench_csr_selection_s("bayesC", -1),
  bench_csr_selection_s("bayesC", -0.5),
  bench_csr_selection_s("bayesR", NULL),
  bench_csr_selection_s("bayesR", -1),
  bench_csr_selection_s("bayesR", -0.5)
)
```

`summarise_architecture()` is a post-hoc descriptive diagnostic for
the relationship between posterior marker signal and marker heterozygosity. It
does not change backend return fields and is not an MCMC estimate of a BayesS
selection parameter.

## Raw C++ Return Slots

Slot numbers here are zero-based native indices. R formatters index them as
one-based list positions.

### `csr_bayesc`

Legacy base slots:

| slot | meaning |
| --- | --- |
| 0 | `bm` |
| 1 | `dm` |
| 2 | `wy` |
| 3 | `r` |
| 4 | final sampled `b` |
| 5 | final sampled `d` |
| 6 | marker index |
| 7 | `vbs` trace |
| 8 | `vgs` trace |
| 9 | `ves` trace |
| 10 | `covb` |
| 11 | `covg` |
| 12 | `cove` |
| 13 | final `vb` |
| 14 | final `vg` |
| 15 | final `ve` |
| 16 | final `pi` |
| 17 | posterior mean `pim` |
| 18 | diagnostics/reserved |
| 19 | sample-count diagnostics |
| 20 | `vle` trace |
| 21 | `vld` trace |
| 22 | `pis` trace or LD-swap diagnostics, detected by formatter |

When chain summaries are returned, slots 23-28 are `bm_sd`, `bm_min`,
`bm_max`, `dm_sd`, `dm_min`, and `dm_max`. When compact chains are kept, slots
29-31 are flattened chain `dm`, chain `bm`, and chain LD-swap diagnostics.
When sampled `selection_s` is enabled, native CSR BayesC appends
`selection_s_trace` and `selection_s_acceptance` after the usual slots. With
`keep_chains = TRUE`, it also appends flattened chain-level `selection_s`
traces and chain-level `selection_s_acceptance`.

### `csr_scheduled_bayesc`

Slots 0-22 match `csr_bayesc`, except LD-swap is unsupported and slot 22 is the
`pis` trace. Slots 23-28 are chain summaries when `nchains > 1`.
`keep_chains = TRUE` is rejected by the R wrapper for scheduled CSR BayesC.

### `csr_bayesr`

This backend returns an `Rcpp::List` with named fields, including:

`bm`, `dm`, `wy`, `r`, `b`, `component`, `vbs`, `vgs`, `ves`, `covb`, `covg`,
`cove`, `vb`, `vg`, `ve`, `pi`, `pim`, `vle`, `vld`, `bm_sd`, `bm_min`,
`bm_max`, `dm_sd`, `dm_min`, `dm_max`, `dm_component_mean`, `comp_prob`,
`ncomp`, `mixture_var`, optional `updateE_diagnostics`, optional `ld_swap`,
optional `chains`, and optional `ld_swap_chains`.

The R formatter resets `dm` to `1 - comp_prob[[trait]][, "component_0"]` and
validates marker-by-component probabilities.

BayesR component-probability columns use `component_0` for the null component;
therefore `dm = 1 - P(component_0)`.

Fixed and sampled `selection_s` are supported for this annotation-unaware CSR BayesR
backend. The R wrapper aligns `Glist$maf` to `Glist$rsidsLD` with
`match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])`, computes
`h = pmax(2p(1-p), 1e-8)`, and passes either `h^(selection_s + 1)` to native
code as `selection_s_prior_scale` for fixed S or aligned `log(h)` for sampled
S. Native BayesR uses the current scale in non-null component posterior
weights, conditional effect draws, the global `vb` variance update through
`sum b_j^2 / (mixture_var_m * scale_j)`, and active/null LD-swap prior density
terms.

### `bed_bayesc`

Legacy slots:

| slot | meaning |
| --- | --- |
| 0-6 | `bm`, `dm`, `wy`, `r`, final `b`, final `d`, marker index |
| 7-9 | `vbs`, `vgs`, `ves` |
| 10-15 | `covb`, `covg`, `cove`, final `vb`, final `vg`, final `ve` |
| 16-17 | final `pi`, posterior mean `pim` |
| 18 | diagnostics: log-CPO, mean log-CPO, mean seconds, max seconds |
| 19 | sample-count diagnostics |
| 20-21 | `vle`, `vld` |
| 22 | `pis` trace |
| 23-28 | `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`, `dm_max` |

BED BayesC does not return compact per-chain objects or LD-swap diagnostics.

### `bed_bayesr`

Slots 0-21 match `bed_bayesc` except slot 5 is the final component state.
Additional slots:

| slot | meaning |
| --- | --- |
| 22 | flattened `comp_prob` by trait |
| 23-28 | `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`, `dm_max` |
| 29 | `dm_component_mean` |

BED BayesR does not return compact per-chain objects or LD-swap diagnostics.

### `csr_prior_bayesc`

Slots 0-17 match `csr_bayesc`. Slot 18 is reserved diagnostics, slot 19 is
sample-count diagnostics, slots 20-21 are `vle` and `vld`, and slot 22 is
LD-swap diagnostics with a formatter flag. R-side multi-chain wrappers append
slots 23-28 for chain summaries and, with `keep_chains = TRUE`, slots 29-31 for
chain `dm`, chain `bm`, and chain LD-swap diagnostics.

### `csr_annot_bayesc`

Slots 0-17 match `csr_bayesc`. Slot 18 is posterior mean `eta_pi`, slot 19 is
posterior mean `eta_vb`, slots 20-21 are `vle` and `vld`, and slot 22 is
LD-swap diagnostics. Multi-chain summary slots 23-28 are standard. Compact
chain slots are 29 `chain_dm_raw`, 30 `chain_bm_raw`, 31 `chain_ld_swap_raw`,
32 `chain_eta_pi_raw`, and 33 `chain_eta_vb_raw`.

### `csr_group_bayesc`

Slots 0-21 follow the BayesC-like CSR convention. Additional slots:

| slot | meaning |
| --- | --- |
| 22 | `group_pi` |
| 23 | `group_vb_multiplier` |
| 24 | `group_nincluded` |
| 25 | `group_size` |
| 26 | LD-swap diagnostics |
| 27-32 | `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`, `dm_max` |
| 33-38 | compact chain `dm`, `bm`, LD-swap, `group_pi`, `group_vb_multiplier`, `group_nincluded` |

### `csr_sbayesrc`

Legacy slots:

| slot | meaning |
| --- | --- |
| 0 | `bm` |
| 1 | `dm = P(component > 0)` |
| 2 | `wy` |
| 3 | `r` |
| 4 | final sampled `b` |
| 5 | final sampled component index |
| 6 | marker index |
| 7-9 | `vbs`, `vgs`, `ves` |
| 10-15 | `covb`, `covg`, `cove`, final `vb`, final `vg`, final `ve` |
| 16 | final `pi0`/active summary |
| 17 | posterior mean `pim` |
| 18 | flattened posterior mean `alpha` |
| 19 | posterior mean `sigmaSqAlpha` |
| 20-21 | `vle`, `vld` |
| 22 | flattened `comp_prob` |
| 23 | `ncomp` |
| 24 | LD-swap diagnostics |
| 25-30 | `bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`, `dm_max` |
| 31-36 | compact chain `dm`, `bm`, LD-swap, `comp_prob`, `alpha`, `sigmaSqAlpha` |
| 37-38 | sampled `selection_s_trace`, sampled `selection_s_acceptance` |
| 39-40 | compact chain sampled `selection_s_trace`, sampled `selection_s_acceptance` |

The R formatter derives `dm_component_mean` from `comp_prob` as posterior mean
zero-based component index, matching BayesR output semantics.

SBayesRC component-probability columns are named by gamma values. The null
component column is `gamma_0.00`; therefore `dm = 1 - P(gamma_0.00)`.

Fixed and sampled `selection_s` are supported for this annotation-aware CSR
SBayesRC backend. The R wrapper reuses the CSR BayesC/BayesR MAF alignment
helper: `match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])`, followed by
`h = pmax(2p(1-p), 1e-8)`. For fixed S it passes `h^(selection_s + 1)`. For
sampled S it passes aligned `log(h)` and native code computes the current
dynamic scale as `exp((S_t + 1) * log(h_j))`. Native SBayesRC uses this scale
in non-null component posterior weights, conditional effect draws, the global
`vb` update through `sum b_j^2 / (gamma_m * scale_j)`, and active/null LD-swap
Gaussian prior-density terms. Annotation/alpha updates operate on component
labels and do not require MAF-scale changes.

## Final R Fit Fields

All audited backends expose `fit$dm`, `fit$bm`, `fit$input`, `fit$vbs`,
`fit$vgs`, `fit$ves`, `fit$vle`, `fit$vld`, and covariance fields where native
code returns them.

CSR BayesC, CSR BayesR, annotation-aware CSR BayesC-like backends, and
CSR SBayesRC support compact `fit$chains` when `keep_chains = TRUE`. Chain
entries use `fit$chains[[trait]][[chain]]$dm` and `$bm` as named marker vectors.

LD-swap/MH-capable CSR backends expose aggregate diagnostics in `fit$ld_swap`;
with compact chains they also expose `fit$ld_swap_chains` and chain-level
`$ld_swap`.

BayesR-like fits expose `fit$comp_prob`; CSR BayesR, BED BayesR, CSR SBayesRC,
and BED BayesRC expose `fit$dm_component_mean`. CSR BayesR, CSR SBayesRC, and
BED BayesRC also expose `fit$ncomp` where native code returns component counts.

Annotation-aware fits expose standardized metadata:

- `fit$input$annotation_model`
- `fit$input$annotations = TRUE`
- `fit$annotation` when an annotation matrix or group input is available
- fixed prior: `fit$annotation_prior`
- learned effects: `fit$annotation_effects`, `fit$eta_pi`, `fit$eta_vb`
- group priors: `fit$annotation_pi`, `fit$annotation_variance`,
  `fit$annotation_summary`, `fit$group_pi`, `fit$group_vb_multiplier`,
  `fit$group_nincluded`, `fit$group_size`
- SBayesRC and BED BayesRC: `fit$annotation_pi`, `fit$annotation_summary`,
  `fit$annotation_effects`, `fit$annotation_variance`, `fit$alpha`,
  `fit$sigmaSqAlpha`

Annotation-unaware CSR and BED fits expose `fit$input$annotations = FALSE`.

## Common Quantity Comparison

| quantity | computed where | exposed where | dimensions and names | chains | LD-swap |
| --- | --- | --- | --- | --- | --- |
| `dm` | all backends | all formatted fits | marker x trait matrix with marker row names and trait column names | compact chain named marker vectors where supported | unchanged |
| `bm` | all backends | all formatted fits | marker x trait matrix | compact chain named marker vectors where supported | unchanged |
| `dm_sd/min/max` | multi-chain backends | all audited backends when native summaries are returned | marker x trait matrices matching `dm` | summary only for scheduled/BED; compact chains for supported CSR | unchanged |
| `bm_sd/min/max` | multi-chain backends | all audited backends when native summaries are returned | marker x trait matrices matching `bm` | summary only for scheduled/BED; compact chains for supported CSR | unchanged |
| `vbs/vgs/ves` | all audited native backends | all formatted fits | iteration x trait traces | CSR BayesR additionally stores chain traces in compact chains | unchanged |
| `vle/vld` | all audited native backends | all formatted fits | iteration x trait traces | CSR BayesR additionally stores chain traces in compact chains | unchanged |
| `covb` | all audited native backends | all formatted fits | trait x trait matrix | no chain-level covariance convention | unchanged |
| `pi/pim` | all BayesC and BayesR-like native backends | all formatted fits | scalar/vector/matrix depending method | BayesR compact chains include final/mean pi | unchanged |
| `input` | R wrappers | all formatted fits | named metadata list | includes `nchains` and `keep_chains` where meaningful | includes LD-swap settings where meaningful |
| `chains` | CSR compact-chain backends | CSR BayesC/BayesR and annotation-aware CSR when `keep_chains = TRUE` | trait list of chain lists | named marker vectors | chain-level LD-swap when supported |
| `ld_swap` | LD-swap-capable CSR backends | capable CSR fits when enabled | trait x attempted/accepted/rate data frame | `ld_swap_chains` with compact chains | n/a |

`summarise_posterior()` and `plot_posterior()` operate on standard
trace fields including `vbs`, `vgs`, `ves`, `vle`, and `vld`.

`summarise_architecture()` operates on fitted `dm` and `bm` marker
summaries plus user-supplied MAF or heterozygosity values. It reports a
post-hoc regression slope named `selection_s_posthoc` and should not be treated
as a sampled or fixed sampler parameter.

## Gaps and Inconsistencies

| gap | classification | disposition |
| --- | --- | --- |
| Raw legacy slot layouts were undocumented and fragile. | documentation gap | Documented here. |
| Annotation-unaware `fit$input` did not consistently include `annotations = FALSE`. | bug | Fixed in R metadata formatting. |
| CSR SBayesRC returned `comp_prob` but did not expose BayesR-like `dm_component_mean`. | bug | Fixed by deriving posterior mean component index from `comp_prob` in the formatter. |
| BED and scheduled CSR backends compute chain summaries but do not return compact `fit$chains`. | intentional difference / future enhancement | Documented; changing requires native return expansion and public behavior decisions. |
| BED backends do not support LD-swap/MH. | intentional difference | Documented. |
| BED backends do not support `chain_seeds`. | intentional difference / future enhancement | Documented. |
| Scheduled CSR BayesC rejects `keep_chains = TRUE` and LD-swap/MH. | intentional difference / future enhancement | Documented. |
| `src/st_sbayesrc_omp_csr_annot.cpp` is present but is not part of the current public backend set. | documentation gap | Documented as out of scope for current public audit. |
| Full cross-backend field inventory was not tested directly. | missing test | Added `test-stblr-backend-field-inventory.R`. |

## Phase 17A remaining-backend computation inventory

The active non-packed-BED numerical surface comprises seven canonical scalar
CSR routes, three internal block-eigen operator routes, five public legacy
multivariate algorithm routes, and two native-only multivariate exports.
Preparation utilities and shared sampler helpers are not independent routes.

Canonical CSR routes own typed contexts/results and named raw conversion.
Block-eigen routes reuse the canonical scalar cores after fit-local dense-block
construction but remain experimental. The multivariate family retains legacy
monolithic/positional boundaries; `mtblr_cpg_omp` additionally derives a seed
from the OpenMP worker identity. See `blr_framework_phase17a_report.md` for the
complete reachability, ownership, RNG, reference, and disposition matrices.

### Phase 17B default multivariate execution

Authoritative `mtblr` owns one fit-local `std::mt19937`, trait-major
marker/effect/state/residual arrays, covariance matrices, and one serial Gibbs
loop. It consumes in-memory dense summaries converted to nested vectors,
produces a 20-position native result, and relies on R for naming/orientation and
BayesC truncation to 18 base fields. Three native/formatted fixture pairs use a
narrow `1e-12` tolerance for Armadillo ulp variation and exact structure to
freeze this behavior. Unconditional set-local and latent `B` updates mean
`updateB=FALSE` is not honored; this is a correction prerequisite, not changed
by the audit.

## Validation Pointers

## Phase 17G current multivariate inventory

| Route | Reachability | Status | Statistical role | Representation | Action |
|---|---|---|---|---|---|
| `mtblr()` | public through `sblr()` | authoritative | corrected MT BayesC reference | dense summary representation | retain |
| `mtblr_eigen()` | native/internal | research | historical eigen evidence | legacy eigen convention | retain temporarily |
| `mtblr_cpg_omp_csr()` | native/internal | research | sparse research evidence | local noncanonical `LDCSR` | retain temporarily |

Phase 17G removed `mtblr_cpg()`, `mtblr_cpg_arma()`, `mtblr_cpg_omp()`, and
`mtblr_hybrid()`. They are historical dispositions, not active routes. The
retained local `LDCSR` assumes one shared operator and is not a candidate
canonical representation.

The lightweight consistency checker validates common marker matrices, optional
`vle`/`vld` traces, chain summaries, compact chain marker names, chain
component probabilities, and LD-swap diagnostics when present. It remains
optional for compatibility objects that lack newer fields.

## Design History

Chain-summary harmonization across CSR and BED backends (shared
`src/st_chain_utils.h` seed/task helpers, canonical exact-CSR and BED
scheduled-chains backends, standard `bm_sd/bm_min/bm_max/dm_sd/dm_min/dm_max`
fields) followed the design work in
`docs/dev/archive/stblr_backend_harmonization_design.md` and
`docs/dev/archive/stblr_csr_multichain_design.md`. Those recommendations are
now implemented; the archived documents are kept for historical rationale
only, not as a current reference.

## Phase 17C corrected default multivariate execution

Only the authoritative public `mtblr` body changed. Set-local `sampleBset` and
`sampleB_latent`, plus the existing global `sampleB`, all respect `updateB`.
Post-burn accumulation begins at `it >= nburn`; marker thinning uses
`(it - nburn) % nthin`; and `bm`/`dm`, `covb`, `covg`, `cove`, and `pim` use
their own positive contribution counts. Disabled covariance/probability updates
retain zero posterior summaries. Native positions 0--19 and R naming,
orientation, and conditional correlation fields remain unchanged. Inputs remain
in-memory with no numerical disk I/O.

## Phase 17E typed default multivariate execution

`mtblr()` now constructs immutable borrowed views for `wy/ww/yy`, sparse nested
`XXvalues/XXindices`, sample counts, models, sets, and covariance priors. Small
execution controls are copied. Existing by-value mutable `b/B/E/pi` inputs move
into `MtDefaultInitialState`; `run_mt_default_core()` owns all workspaces and
its fit-local RNG and returns an owning `MtDefaultCoreResult`. Dense `XX` is not
copied by the typed boundary. The inline positional finalizer reads the result
and borrowed `wy`; no MCMC-time I/O, native aggregation, or binding conversion
was added.

## Phase 17F finalized default multivariate result

`MtDefaultCoreResult` retains raw marker/covariance/probability accumulators,
traces, final state, and contribution counts. `finalize_mt_default_result()`
moves those buffers, performs all posterior divisions, and returns
`MtDefaultFinalResult` with finalized `bm/dm`, covariance means,
`pi_final/pi_mean`, traces, final state, legacy zero fields, and retained-count
metadata. `mtblr()` only allocates/shapes/copies the legacy positions, including
borrowed `wy`. No extra dense LD copy or MCMC-time I/O is introduced.

Future canonical MT LD execution must reuse the scalar CSR ABI (zero-based
`uint64` row offsets, canonical column indices/value buffers, marker diagonal,
summary statistics, sample sizes, immutable borrowing) rather than create an
MT-specific format. An MT bundle must permit one operator per trait: identical
row/column structure may be shared with trait-specific values, while different
cohorts/ancestries/panels/thresholds may use independent structures without a
union sparsity pattern. The same rule applies to the future canonical scalar
block-eigen representation: one compatible decomposition per trait, shareable
only for truly identical LD.
## Phase 17H CSR representation

Canonical scalar CSR backends now share one owning definition,
`sblr::core::SparseLdCsrStorage`, through the temporary `STLDCSR` alias.
Ordinary CSR BayesC actively borrows `SparseLdCsrView`; the other scalar CSR
families retain model-specific adapters. `mtblr_cpg_omp_csr()` retains its
local one-LD `LDCSR` strictly as noncanonical internal research evidence.
## Phase 17I internal MT CSR route

`mtblr_csr_internal()` is a registered maintenance route using canonical
`SparseLdCsrStorage/View`, one operator per trait, the shared corrected MT core,
the existing finalizer, and the shared legacy adapter. It is not routed by
`sblr()` and is distinct from noncanonical `mtblr_cpg_omp_csr()` research.
Phase 17J adds the public serial `mtblr_csr()` adapter. It uses the canonical
Phase 17I trait-specific CSR core and adds no Gibbs loop, OpenMP path, or
research-backend dependency.
Phase 17K terminology: scalar “block eigen” routes construct operators using
eigendecomposition or ridge regularization, then store reconstructed dense
within-block cross-products as float-packed upper triangles. BayesC, BayesR,
and SBayesRC borrow the same canonical view; no eigensystem persists at runtime.
# Internal MT block-eigen execution (Phase 17L)

`mtblr_block_eigen_internal()` builds canonical `BlockEigenStorage` owners from BED descriptors, retains matching transformed `wy`, and delegates through `MtBlockEigenDataView` to the Phase 17I shared MT BayesC core. It does not use the retained research `mtblr_eigen()` implementation.

# Public MT block-eigen execution (Phase 17M)

`mtblr_block_eigen()` validates same-BED provenance and delegates to
`mtblr_block_eigen_raw_internal()`. The raw adapter builds and samples once;
summary `ww` never enters the runtime operator.

# Planned internal MT packed-BED computation (Phase 17N)

Phase 17O adds internal backend `mt_bed_bayesc`. It uses one
`PackedBedMatrix`, one `BedPackedGenotypeView`, a reusable decoded double
marker, `arma::mat` sample residuals, full-`E` masked-latent marker algebra,
and final `X'R` reconstruction into `mtblr_raw` version 1. It is not routed
through any public function.
## Phase 17P public MT BED adapter

`mtblr_bed()` performs Glist/BED/sample/marker alignment, phenotype validation
and centering, public prior/initialization preparation, analytical memory
estimation, and metadata formatting in R. All individual-level MCMC computation
remains solely in the unchanged Phase 17O `mtblr_bed_internal()` backend.

Phase 17Q adds no backend computation. Its future plan prepares owner/view,
phenotype, maps, `wy`, and order once, then invokes the unchanged Phase 17O core
once per independent chain and aggregates only after every chain succeeds.

Phase 17R implements that internal computation: one prepared dataset feeds one
unchanged `run_mt_bed_bayesc_core()` invocation per chain. Static OpenMP is
chain-level only; aggregation invokes the existing finalizer exactly once.

Phase 17S adds no backend computation. The public adapter validates controls,
estimates memory, calls `mtblr_bed_chains_internal()` once, enriches metadata,
and delegates formatting to `.as_mtblr_fit()`.

MT BED formal convergence diagnostics remain unavailable. Typed chains already
hold B/G/E diagonal iteration traces; Phase 17U is reserved for their post-burn
rank/folded R-hat and ESS/MCSE engine.
