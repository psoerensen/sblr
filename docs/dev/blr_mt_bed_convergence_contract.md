# MT BED convergence-diagnostic contract

## 1. Purpose

This document specifies formal convergence diagnostics for individual-level MT
BED multichain fits. Phase 17U activates the internal Tier 1 engine; the public
`mtblr_bed()` adapter does not call it until a later activation phase.

## 2. Current availability

Each typed Phase 17R chain owns trait-major `vbs`, `vgs`, and `ves` matrices of
length `nit + nburn`. They contain the diagonal of B when B is updated, the
diagonal of derived G every iteration, and the diagonal of E when E is updated.
Compact chains expose these traces only when `keep_chains=TRUE`; before
aggregation the typed chain results contain them regardless of that choice.
Phase 17U extracts exactly the post-burn portion into a binding-neutral Tier 1
bundle independently of compact-chain retention.

## 3. Current limitations

There are no iteration traces for marker latent/effective effects, inclusion
states, off-diagonal B/G/E entries, model probabilities, full covariance
matrices, or marker residual scores. Pooled `vbs/vgs/ves` are iterationwise
chain means and are not valid formal diagnostic inputs. `bm_sd` and `dm_sd` are
stability summaries of chain posterior means, not R-hat or posterior SD.

## 4. Statistical authority

The target is Vehtari, Gelman, Simpson, Carpenter, and Bürkner (2021),
*Bayesian Analysis* 16(2), 667–718, DOI 10.1214/20-BA1221: rank-normalized
split R-hat, folding, localization for tail ESS, and modern multichain ESS.
Classical unsplit Gelman–Rubin and `coda::gelman.diag()` are not targets.

## 5. Post-burn policy

For public one-based iterations use exactly `nburn + 1` through
`nburn + nit`. Burn-in is excluded. Every chain must have exactly `nit`
post-burn draws.

## 6. Thinning policy

Use every available post-burn trace draw. Do not apply additional diagnostic
thinning. Current `nthin` governs marker posterior-summary retention, not the
low-dimensional iteration traces. Future selected-marker diagnostic draws are
also recorded every post-burn iteration.

## 7. Chain splitting

For N post-burn draws set `half=floor(N/2)`. Split into draws `1:half` and
`(N-half+1):N`. If N is odd, discard the central draw. The result is `2*nchains`
chains of length `half`; never pad, recycle, or interpolate.

## 8. Rank normalization

Pool one scalar quantity across split chains, assign deterministic average
ranks for ties, and transform rank r among S values by
`qnorm((r - 3/8)/(S + 1/4))`. Ties and discrete values are never jittered.
The plus sign is the audited `posterior` 1.6.1/Blom back-transformation
`S - 2c + 1` with `c=3/8`; Phase 17T's minus sign was corrected in Phase 17U.

## 9. Rank R-hat

For M split chains of length N, let W be the mean within-chain sample variance,
`B=N*var(split-chain means)`, and
`var_plus=((N-1)/N)*W+B/N`. Then
`rhat_rank=sqrt(var_plus/W)` on rank-normalized values.

## 10. Folded R-hat

Take absolute deviations of unsplit post-burn draws from their pooled median,
then split, rank-normalize, and apply the same R-hat calculation to obtain
`rhat_folded`.

## 11. Final R-hat

Report `rhat=max(rhat_rank,rhat_folded)`. Reporting only the smaller component
is forbidden.

## 12. Bulk ESS

`ess_bulk` is the multichain ESS of rank-normalized split draws. Phase 17U
reproduces `posterior` 1.6.1's FFT autocovariance normalization, mean
within-chain variance, between-chain contribution, rho sequence, Geyer initial
positive sequence, and initial monotone sequence. Its numerical-stability rule
sets `tau >= 1/log10(M*N)`, so ESS is bounded above by
`M*N*log10(M*N)` rather than by nominal draws. Whether that bound was applied
is recorded per quantity. Define `ess_bulk_per_chain=ess_bulk/nchains`.

## 13. Tail ESS

Compute pooled raw-scale q05 and q95. Form indicators `draw<=q05` and
`draw<=q95`, calculate split-chain ESS for each, and use their minimum.
Define `ess_tail_per_chain=ess_tail/nchains`. Constant indicators caused by
ties are unavailable with an explicit status; bulk ESS is not substituted.

## 14. Mean ESS

`ess_mean` applies the selected split-chain autocorrelation estimator to raw
post-burn draws for the posterior-mean estimand.

## 15. MCSE

Phase 17U implements only `mcse_mean=posterior_sd/sqrt(ess_mean)` and
`mcse_mean_over_sd=mcse_mean/posterior_sd`. `ess_sd`, `mcse_sd`, quantile MCSE,
and median MCSE are deferred.

## 16. Constant traces

All-constant traces receive status `constant` and NA R-hat, ESS, and MCSE. They
must not be reported as R-hat 1, nominal ESS, or MCSE zero. A constant chain
while others vary receives `constant_chain_mismatch`, retains calculable
diagnostic values for inspection, and forces warning status; it is not an
automatic perfect result.

## 17. Fixed parameters

When `updateB=FALSE`, B rows have status `not_updated`; when `updateE=FALSE`, E
rows have status `not_updated`. Their current zero trace placeholders are never
diagnosed as stochastic. G remains eligible because it is derived every
iteration.

## 18. Insufficient chains

With fewer than two chains use `unavailable_single_chain`. A two- or
three-chain result may be computed but uses
`computed_fewer_than_four_chains`; record the advisory
`fewer_than_recommended_chains=TRUE`.

## 19. Insufficient draws

Fewer than four post-burn draws per original chain receives
`insufficient_draws`. Four or five draws permit split R-hat but leave ESS and
MCSE unavailable, with status `computed_partial`. Six draws permit the selected
ESS kernel because each split chain then has three draws. Metric-level
availability fields preserve valid R-hat rather than discarding it. Nonfinite
values receive `nonfinite`.

## 20. Thresholds

Advisory defaults are R-hat > 1.01, bulk ESS < `100*nchains`, tail ESS <
`100*nchains`, and `mcse_mean_over_sd > 0.05`. These do not prove convergence.

## 21. Warning policy

Automatic diagnostics with one chain record unavailability without warning.
With at least two chains, issue at most one main-thread warning per fit if any
threshold fails. An explicit request that cannot be computed warns once. The
warning reports affected quantity count, maximum R-hat, minimum bulk and tail
ESS per chain, maximum relative MCSE, and the detailed-output location.

## 22. Tier 1 scope

The shared Phase 18 implementation diagnoses `vbs[trait]`, `vgs[trait]`,
`ves[trait]`, `vle[trait]`, and `vld[trait]`, subject to update ownership and
scientific applicability. It consumes only per-chain post-burn traces. For
MTBLR, the first three are the B/G/E diagonals; full covariance matrices remain
separate `cov_*_mean` and `cov_*_final` fields. `vle` is the represented
operator-diagonal contribution of current effective effects and
`vld = vgs - vle` at every iteration.

## 23. Tier 2 scope

Tier 2 requires new lower-triangle B/G/E traces in column-major lower-triangle
order `(trait1,trait1),(trait2,trait1),...`, with public pair names, plus
`pi_null` and `pi_active`. Pattern-specific `pi[model_name]` is opt-in and
memory-controlled. Diagonal E exposes diagonal entries only; disabled updates
use `not_updated`.

## 24. Tier 3 scope

Tier 3 records post-burn `b[marker_id,trait]` and `d[marker_id,trait]` for
explicit marker IDs or one-based indices, preserving requested order. It is
opt-in, validates all selections, and warns against a memory threshold rather
than silently truncating. All-marker default traces are forbidden.

## 25. `keep_chains` independence

Formal diagnostics must work with `keep_chains=FALSE`. Compact chain records
and diagnostic trace capture are separate. The diagnostic engine consumes
typed traces before compact-chain discard and does not retain marker means or
final states merely to diagnose Tier 1.

## 26. Diagnostic trace bundle

Binding-neutral `MtBedConvergenceTraceBundle` metadata contains chain
count, post-burn length, ordered quantity descriptors, and contiguous
quantity-major values with chain then iteration varying within quantity. Tier
1 contains only post-burn B/G/E diagonals. It is Rcpp-free, read-only during
diagnosis, independent of `keep_chains`, and discardable after summarization.

## 27. Output schema

The internal helper additively retains raw version 1 at
`raw$diagnostics$convergence`; future public activation formats it at
`fit$convergence`. Include version, requested/computed flags, scope,
overall status, chain/draw counts, thresholds, availability, summary, overview,
warning messages, and trace-retention policy. Optional retained post-burn
diagnostic traces belong at `fit$convergence_traces`, never `fit$chains`.

## 28. Summary table

One rectangular row per scalar contains `quantity`, `group`, `trait`, `trait2`,
`marker_id`, `model_name`, `updated`, `status`, `rhat`, `rhat_rank`,
`rhat_folded`, `ess_bulk`, `ess_bulk_per_chain`, `ess_tail`,
`ess_q05`, `ess_q95`, `ess_tail_per_chain`, `ess_mean`, `posterior_sd`, `mcse_mean`,
`mcse_mean_over_sd`, `nchains`, `draws_per_chain`,
`split_draws_per_chain`, metric-level availability fields, the ESS-stability
bound indicator, and four threshold flags. Unused identities are NA.

## 29. Overview

Return `overall_status`, `n_computed`, `n_unavailable`, `n_flagged`, maxima and
quantity names for R-hat and relative MCSE, minima/per-chain minima and names
for bulk/tail ESS, and `fewer_than_four_chains`. Status is one of
`not_requested`, `unavailable`, `ok`, `warning`, or `partial`; never a
definitive `converged` status.

## 30. Memory

For `nit` post-burn draws, C chains, T traits, Q=T(T+1)/2, K selected markers,
and P patterns:

```text
tier1_trace_bytes       = 8*C*nit*3*T
covariance_trace_bytes  = 8*C*nit*3*Q
probability_trace_bytes = 8*C*nit*2
full_pi_trace_bytes     = 8*C*nit*P
selected_b_trace_bytes  = 8*C*nit*K*T
selected_d_trace_bytes  = 4*C*nit*K*T
```

Process one quantity at a time so diagnostic workspace is O(C*nit). Captured
storage, workspace, optional retained traces, and summary output are reported
separately and remain analytical estimates, not measured RSS. Diagnostic
memory never charges packed genotype storage.

## 31. Complexity

For S usable values per scalar, rank sorting is O(S log S), R-hat after ranking
is O(S), FFT autocorrelation/ESS is O(S log S), and workspace is O(S).
Thousands of marker/model quantities are not negligible.

## 32. Future public controls

After `keep_chains`, prefer
`convergence=c("auto","none","core")` and `convergence_control=NULL`.
Auto computes Tier 1 for at least two chains and records quiet unavailability
otherwise; none captures nothing; core explicitly requests Tier 1 and warns
once if unavailable. Control defaults are `warn=TRUE`, `rhat_threshold=1.01`,
`ess_per_chain_threshold=100`, `mcse_mean_over_sd_threshold=0.05`, and
`keep_traces=FALSE`. The future engine uses base R/stats without a mandatory
runtime dependency; `posterior` is optional for development comparison only,
and coda is not the target.

## 33. Phase 17U implementation

Phase 17U captures typed post-burn per-chain B/G/E diagonal traces independently
of `keep_chains` and implements rank/folded R-hat, bulk/tail/mean ESS, mean
MCSE, statuses, summary, overview, validation, memory estimates, and frozen plus
optional `posterior` oracles without public formals.

## 34. Phase 17V plan

Activate public controls, aggregated warnings, `fit$convergence`, optional
trace retention, memory metadata, documentation, and exact default reduction.

## 35. Extended diagnostic plan

Only later add lower-triangle covariances, null/active probabilities, optional
pattern probabilities, selected-marker b/d traces, and their memory controls.

## Phase 17V public activation

Public `mtblr_bed()` now supports `convergence=c("auto","none","core")` and a
strict `convergence_control` list. Trace capture occurs only for requested
multichain diagnostics or explicit trace retention. `none` returns a validated
zero-row not-requested object; automatic single-chain use returns quiet
unavailability. Warnings are aggregated, advisory, and do not prove
convergence. Diagnostics remain independent of `keep_chains`, use per-chain
post-burn rather than pooled traces, and expose no Tier 2 or marker scope.

## Phase 17W extended-scope boundary

The implementation-ready Tier 2/Tier 3 ordering, deduplication, selection,
safety, output, and staging rules are owned by the extended convergence
contract. Tier 1 mathematics and public behavior remain unchanged, and no
extended diagnostic is implemented in Phase 17W.
## Phase 18 shared-engine ownership

The Phase 17U mathematics is now owned by the family-neutral internal
`.blr_convergence_*` layer. MT BED supplies the same Tier 1 bundle and remains
numerically identical after identifier/metadata normalization. ST adapters use
the same result schema without changing the statistical definitions.

## Phase 21 additive extension

Phase 21 implements this extension through the shared bundle and scalar engine.
MT BED owns strict-lower covariance, probability/annotation, and explicit
selected-marker buffers inside each joint chain. The same public descriptor,
memory-guard, warning, and retention contracts now consume ST trait-by-chain
component, annotation/group, and selected-marker native buffers without
changing the MT recorder or convergence mathematics.

`extended` adds strict-lower raw B/G/E entries, applicable model probability
states, BayesRC coefficient states, and explicitly selected b/d/component
states. The completed-iteration native checkpoint is unchanged. Diagonal-E
off-diagonals are structural rows without stochastic capture; covariance
diagonals remain exclusively in Tier 1.
