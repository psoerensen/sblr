# Unified BLR Framework Phase 17T report

## 1. Executive summary

The MT BED convergence-diagnostic contract is implementation-ready. No
diagnostic implementation is active.

## 2. Repository baseline

Baseline: clean `master` at `661566f9da2882b205ae50e373e65e651e159054`,
R 4.4.1 UCRT, GCC/G++/GFortran 13.2, OpenMP, R BLAS/R LAPACK. The source suite
passed 5,265 expectations with zero failures/warnings and two opt-in skips.
The built check passed with zero errors/warnings and three established notes.
Hosted CI is not locally visible.

## 3. Phase 17S verification

Public `mtblr_bed()` runs joint chains through the Phase 17R chains route,
pools posterior accumulators, uses primary-chain final state, averages traces,
and optionally retains compact chains. Production is unchanged.

## 4. Existing diagnostic inventory

Repository search found legacy scalar `coda::effectiveSize`, Geweke, lag-one
autocorrelation, and simple MCSE summaries, plus documentation disclaimers.
It found no MT BED rank R-hat, folded R-hat, modern ESS, or MCSE engine.

| Quantity | Per-chain iteration trace | Current other representation | Updated/fixed | Formal now | Additional trace | Memory order | Tier |
|---|---|---|---|---|---|---|---|
| B diagonal | `vbs` | final/posterior covariance | `updateB` | possible before discard | no | O(CNT) | 1 |
| G diagonal | `vgs` | final/posterior covariance | derived every iteration | possible before discard | no | O(CNT) | 1 |
| E diagonal | `ves` | final/posterior covariance | `updateE` | possible before discard | no | O(CNT) | 1 |
| B/G/E off-diagonal | none | posterior mean/final only | mixed | no | lower triangle | O(CN T²) | 2 |
| model probabilities | none | final/mean only | `updatePi` | no | null/active or patterns | O(CN) or O(CNP) | 2 |
| marker b/d | none | posterior mean/final only | stochastic | no | selected markers | O(CNKT) | 3 |
| marker stability | none | chain posterior mean SD/min/max | summary | no | selected-marker draws | O(CNKT) | 3 |

## 5. Current trace availability

Typed chains contain trait-major B/G/E diagonal traces of length `nit+nburn`.
B/E are written only when updated; G is written every iteration. `nthin` does
not thin these traces.

## 6. Current trace limitations

No iteration traces exist for markers, off-diagonal covariances, covariance
matrices, or probabilities. Pooled traces are chain means, not diagnostic
draws. Compact per-chain traces are conditional on `keep_chains`.

## 7. Statistical reference

Vehtari et al. (2021), DOI 10.1214/20-BA1221, is authoritative. Classical
unsplit Gelman–Rubin and coda Gelman diagnostics are rejected as targets.

## 8. Post-burn policy

Use public iterations `nburn+1` through `nburn+nit`; exclude all burn-in.

## 9. Split policy

Split each chain into first/last `floor(nit/2)` draws and discard an odd center.

## 10. Rank normalization

Use average ranks for ties and `qnorm((r-3/8)/(S-1/4))`, without jitter.

## 11. R-hat

Compute split rank R-hat and split folded rank R-hat; report their maximum.

## 12. Bulk ESS

Use rank-normalized split draws with multichain autocovariances and Geyer
initial-positive then initial-monotone paired sequences.

## 13. Tail ESS

Use the minimum ESS of raw pooled 5% and 95% threshold indicators.

## 14. Mean ESS and MCSE

Phase 17U targets raw-scale mean ESS, posterior SD, mean MCSE, and relative mean
MCSE. SD and quantile MCSE are deferred.

## 15. Constant/fixed handling

Constant diagnostics are NA with `constant`; disabled B/E are `not_updated`.
Neither is represented as perfect convergence.

## 16. Chain/draw sufficiency

One chain is unavailable; two/three chains compute with an advisory; fewer
than four draws is insufficient; nonfinite draws are unavailable.

## 17. Thresholds

Defaults: R-hat 1.01, bulk/tail ESS 100 per chain, relative MCSE 0.05, and a
fewer-than-four-chain advisory.

## 18. Warning behavior

Warnings are advisory, aggregated once on the main R thread, and never stop a
fit. Automatic single-chain unavailability is quiet.

## 19. Tier 1 diagnostics

First scope is updated B diagonals, derived G diagonals, and updated E
diagonals, named by public trait.

## 20. Tier 2 diagnostics

Lower-triangle B/G/E and null/active probability diagnostics require new
traces. Full-pattern probabilities are opt-in and memory-controlled.

## 21. Tier 3 diagnostics

Selected-marker b/d traces are explicit opt-in by ID or one-based index; all
markers are never retained by default.

## 22. `keep_chains` independence

Diagnostics consume a separate typed post-burn bundle before optional compact
chain discard and therefore do not require `keep_chains=TRUE`.

## 23. Diagnostic computation location

Selected architecture: native typed trace capture, dependency-free auditable R
diagnostic engine, and optional posterior-package test oracle.

## 24. Dependency decision

Base R/stats will implement the authority contract. `posterior` is optional for
development comparison; coda is not the target. DESCRIPTION remains unchanged.

## 25. Output schema

Future additive raw version-1 output is `raw$diagnostics$convergence`; formatted
output is `fit$convergence`; optional traces are `fit$convergence_traces`.

## 26. Summary and overview

The contract defines a rectangular per-quantity table, explicit availability
statuses/flags, aggregate extrema, counts, and non-definitive overall statuses.

## 27. Future public API

Preferred controls are `convergence=c("auto","none","core")` and a named
`convergence_control`, after `keep_chains`. They are not added in Phase 17T.

## 28. Memory scaling

Tier 1 is `8*C*nit*3*T`; covariance, probability, full-pi, and selected-marker
formulas are formalized separately. Workspace processes one quantity at a time.

## 29. Complexity

Ranking and FFT ESS are O(S log S), R-hat after ranking is O(S), and workspace
is O(S) per processed scalar.

## 30. Deterministic oracles

Pure-R fixtures cover well-mixed, location/scale mismatch, drift, positive and
negative autocorrelation, tails, ties, binary, constant, partially constant,
nonfinite, one/two/four chains, odd lengths, short chains, and fixed status.

## 31. Optional reference comparison

Future development may compare to `posterior::rhat`, `ess_bulk`, `ess_tail`,
`ess_mean`, and `mcse_mean` when installed; package checks never require it. A
local optional comparison on the deterministic well-mixed fixture matched
R-hat exactly (0.9746794), mean ESS within 0.5%, MCSE within 0.3%, bulk ESS
within 3%, and tail ESS within 24%; the larger tail difference is retained as
an explicit oracle-tolerance/staging item rather than hidden by a dependency.

## 32. Existing-route protection

All Phase 15/17, MT public, scalar BED, raw schema, wrapper, and native owners
remain unchanged.

## 33. Mutation sensitivity

All 35 required statistical, scope, dependency, production, schema, and CI
mutations are guarded.

## 34. Synthetic benchmark

The synthetic benchmark reports trace/storage formulas and oracle timing only;
it makes no production-runtime claim. Across the requested grid, maxima were
19.2 MB Tier 1, 201.6 MB covariance, 0.64 MB null/active probability, 1.311 GB
full-pi, 6.4 GB selected-b, 3.2 GB selected-d, and 0.32 MB per-quantity
workspace. One hundred deterministic calls took approximately 0.05 s R-hat,
0.03 s bulk ESS, and 0.02 s mean MCSE on this host.

## 35. Installed-check behavior

Numerical contract oracles are portable. Only repository source protection
assertions skip narrowly without a checkout.

## 36. Tests and CI

Phase 17T is included in the exact fast filter. Focused Phase 17T passed 54
expectations. The exact fast tier passed 2,462 expectations with zero failures
or warnings and one established opt-in skip. The final full source suite passed
5,319 expectations with zero failures or warnings and two established opt-in
skips. Installed package tests passed.

## 37. Package check

The built-tarball check passed with zero errors and zero warnings. The three
unchanged classified notes are the long Phase 17C fixture path, installed size
(5.2 MB, including 4.2 MB in `libs`), and legacy scalar-backend `std::cout`
symbols. No new note was introduced.

## 38. Diff hygiene

Final audit confirms contract/docs/tests/tooling plus one CI-filter edit only.
Production, DESCRIPTION, NAMESPACE, Rd, wrappers, registration, schemas,
fixtures, and native sources are unchanged. Compiled/check artifacts are
removed before handoff, and normal versus ignore-EOL statistics are compared.

## 39. Deviations and blockers

No contract deviation or implementation blocker is known. The optional
`posterior` tail-ESS oracle differs by about 24% on the deliberately synthetic
well-mixed fixture; Phase 17U must resolve the exact reference implementation
details before production activation. This is an explicit staging item, not a
Phase 17T contract blocker.

## 40. Recommended next phase

> implement the internal Tier 1 MT BED convergence engine defined by Phase 17T, using post-burn per-chain `vbs`, `vgs`, and `ves` traces independently of `keep_chains`, calculating rank-normalized split and folded R-hat, bulk ESS, tail ESS, mean ESS, and MCSE mean, returning additive internal diagnostic summaries while leaving the public `mtblr_bed()` signature and numerical sampler unchanged.

## 41. Readiness marker

PHASE 17T COMPLETE — MT BED CONVERGENCE-DIAGNOSTIC CONTRACT FORMALIZED
