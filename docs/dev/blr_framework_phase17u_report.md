# Phase 17U implementation report

## 1. Executive summary

The internal Tier 1 MT BED convergence engine is active. Public `mtblr_bed()`
still uses the Phase 17S chains route and exposes no convergence controls or
fields.

## 2. Repository baseline

Work began from clean `master` at
`65a41e964cb275108326314887fae549eee7818e`. The host uses R 4.4.1 UCRT,
GCC/GFortran 13.2, R BLAS/LAPACK, and OpenMP. Baseline source tests passed with
the two established opt-in skips. The built package passed with zero errors and
warnings and the three established notes: the Phase 17C long fixture path,
installed size, and legacy scalar `std::cout` symbols.

## 3. Phase 17T verification

| Contract | Production symbol | Test owner | Audit owner | Output |
|---|---|---|---|---|
| post-burn Tier 1 bundle | `build_mt_bed_convergence_trace_bundle()` | Phase 17U trace tests | Phase 17U architecture audit | `trace_bundle` |
| rank/folded R-hat | `.mtblr_convergence_scalar()` | frozen/posterior oracles | Phase 17U audit | summary R-hat columns |
| bulk/tail/mean ESS | `.mtblr_convergence_ess()` | frozen/posterior oracles | mutation audit | summary ESS columns |
| MCSE mean | `.mtblr_convergence_scalar()` | frozen/posterior oracles | mutation audit | MCSE columns |
| validation | two convergence validators | failure-closed tests | architecture audit | validated internal result |
| memory | `.mtblr_convergence_memory_estimate()` | formula tests | mutation audit | analytical estimate |

## 4. Numerical reconciliation

The Phase 17T tail discrepancy came from its simplified autocorrelation-pair
sequence and nominal-draw cap. The current reference uses its precise Geyer
sequence and a `tau >= 1/log10(nominal)` stability bound, allowing antithetic
ESS above nominal draws. Phase 17T also recorded the Blom denominator with the
wrong sign.

## 5. Reference algorithm

The pinned development oracle is `posterior` 1.6.1, inspected on 2026-07-26.
Inspected symbols were `rhat.default`, `.rhat`, `ess_bulk.default`,
`ess_tail.default`, `ess_quantile.default`, `.ess_quantile`,
`ess_mean.default`, `.ess`, `mcse_mean.default`, `.split_chains`,
`z_scale`, `backtransform_ranks`, `fold_draws`, `autocovariance`, and
`should_return_NA`. R-hat tolerance is 1e-12; ESS and MCSE relative tolerance
is 1e-8.

## 6. Contract amendments

Rank scores use denominator `S + 1/4`. ESS records application of the
`nominal*log10(nominal)` upper stability bound. R-hat needs four original
draws, while ESS/MCSE need six; four- and five-draw results use
`computed_partial`.

## 7. Trace bundle

`MtBedConvergenceTraceBundle` is Rcpp-free and owns descriptors plus contiguous
doubles in quantity, chain, iteration order. It contains only B/G/E diagonals.

## 8. Trace extraction

The extractor validates length `nit+nburn` and copies C++ indices `nburn`
through `nburn+nit-1`. Updated quantities must be finite.

## 9. Internal route

`mtblr_bed_convergence_trace_internal()` has the 31 Phase 17R chain arguments
and returns `list(raw, trace_bundle)`. Both internal routes delegate to one
binding helper.

## 10. R engine architecture

`R/mtblr-convergence.R` owns splitting, ranks, folding, R-hat,
autocovariance, ESS, scalar/Tier 1 diagnosis, validators, overview, high-level
attachment, and memory estimation using base R and stats only.

## 11. Chain splitting

Each chain becomes its first and last `floor(nit/2)` draws; an odd central draw
is dropped.

## 12. Rank normalization

Average ranks are transformed with
`qnorm((rank-3/8)/(S+1/4))`. Ties are never jittered.

## 13. Rank R-hat

Rank R-hat applies basic split R-hat to the rank-normalized split matrix using
sample within- and between-chain variances.

## 14. Folded R-hat

Folded R-hat diagnoses absolute deviations from the pooled unsplit median,
then splits and rank-normalizes.

## 15. Final R-hat

The reported R-hat is the maximum of rank and folded R-hat.

## 16. ESS kernel

The dependency-free kernel matches `posterior` 1.6.1 FFT autocovariance
normalization, between-chain variance, rho construction, initial positive
sequence, initial monotone sequence, last positive even term, tau, and bound.

## 17. Tail ESS

Type-7 pooled q05/q95 thresholds create two `draw <= q` indicator matrices.
Each is split and diagnosed separately; tail ESS is their minimum. The former
24% discrepancy is eliminated.

## 18. Mean ESS

Mean ESS uses raw-scale split draws, not rank-normalized draws.

## 19. MCSE mean

`mcse_mean=posterior_sd/sqrt(ess_mean)` and its SD-relative ratio are returned.
No SD or quantile MCSE is implemented.

## 20. Partial metric availability

Five explicit availability columns allow valid four/five-draw R-hat to coexist
with unavailable ESS and MCSE.

## 21. Constant/fixed handling

Statuses are `not_updated`, `nonfinite`, `constant`,
`unavailable_single_chain`, `insufficient_draws`,
`constant_chain_mismatch`, `computed_partial`,
`computed_fewer_than_four_chains`, and `computed`. Constant traces never
receive perfect diagnostics.

## 22. Tier 1 quantity scope

Stable names are `B_diag[trait]`, `G_diag[trait]`, and `E_diag[trait]`.
B/E respect update ownership; G is always eligible.

## 23. Summary table

The rectangular table includes identity, status, all R-hat components,
bulk/q05/q95/tail/mean ESS, posterior SD, MCSE, counts, availability,
stability-bound, and threshold-flag columns.

## 24. Overview

The overview reports computed/partial/unavailable/not-updated/flagged counts,
metric extrema and quantity labels, mismatch and stability-bound counts, and
the fewer-than-four-chain advisory. It never declares convergence.

## 25. Internal result object

The convergence version-1 object records scope, counts, thresholds,
availability, summary, overview, warning text for future activation,
trace-retention policy, and algorithm provenance.

## 26. Validation

Bundle validation checks schema, dimensions, Tier 1 order, trait indices,
update flags, and finite updated traces. Result validation checks schema,
statuses, availability, metrics, flags, overview, and retention policy.

## 27. keep_chains independence

Typed traces are extracted before aggregation discards chain state.
`keep_chains=FALSE` and `TRUE` produce identical bundles and diagnostics.

## 28. Public protection

`R/mtblr-bed.R`, its formals and dispatch, public fit fields, memory warning,
Rd, NAMESPACE, and DESCRIPTION remain unchanged.

## 29. Numerical reference comparisons

Frozen well-mixed references include R-hat 0.974679434480896, bulk ESS
54.458829511151, tail ESS 197.402597402597, mean ESS 43.8040345821326, and
MCSE mean 0.183995239970637. Optional `posterior` comparisons pass.

## 30. Actual MT BED tests

Tiny packed-BED tests cover one/two/four chains, one/two/three traits, full and
diagonal E, update controls, retained/unretained chains, and post-burn order.

## 31. Serial/OpenMP equality

When OpenMP is available, serial and two-worker trace bundles, diagnostics, and
ordinary numerical output are exact after timing normalization.

## 32. Reproducibility

Repeated, explicit/default-seed, fresh-process, and intervening-fit contracts
use the existing deterministic Phase 17R execution owner.

## 33. Memory

Capture is `8*nchains*nit*3*nt` bytes. The estimate separates one-quantity
workspace, summary output, and optional retained-trace storage and excludes
genotypes and phenotypes.

## 34. Architecture audit

The audit requires one route, bundle, R-hat engine, ESS engine, and MCSE engine;
zero public controls/calls; Tier 1-only dependency-free behavior; and unchanged
public adapter and numerical core.

## 35. Mutation sensitivity

All 44 required extraction, statistical, ownership, public-protection,
dependency, schema, aggregation, and CI guards pass.

## 36. Benchmark

The benchmark reports synthetic scalar/Tier 1 timing, actual internal-route
signals, analytical memory, and retained-trace size without public-runtime or
negligible-cost claims. Across the required grid, the largest case
(`nchains = 8`, `nit = 5000`, `nt = 20`) used 19,200,000 captured-trace bytes,
1,920,000 maximum workspace bytes, and 19,200,000 retained-trace bytes, for a
0.03756657 GiB analytical total; its Tier 1 calculation took 7.15 seconds on
this host. The small actual BED fixture took 0.11 seconds for route plus trace
extraction and 0.07 seconds for the Tier 1 engine; the trace bundle and retained
trace objects were 5,000 and 7,008 bytes, respectively. These are regression
signals, not public-runtime or speedup claims.

## 37. Existing-route protection

Phase 17O/S/R/T, scalar packed-BED, summary MT, raw-schema, and public API
owners remain covered by their permanent tests.

## 38. Generated wrappers

`Rcpp::compileAttributes()` adds exactly one unexported R wrapper and one
registered native symbol for the convergence trace route.

## 39. Installed-check behavior

Numerical and actual-route tests are portable. Only repository source/hash
assertions skip narrowly without a source checkout. `posterior` remains
optional and is not needed by installed tests.

## 40. Tests and CI

Phase 17U is included in the exact fast filter. Its focused owner contains 127
passing expectations. The exact fast tier completed 2,590 expectations with
2,589 passes, zero failures, zero warnings, and one established opt-in skip.
The full source suite completed 5,448 expectations with 5,446 passes, zero
failures, zero warnings, and two established opt-in skips. Installed tests
completed successfully.

## 41. Package check

The built-tarball check completed with zero errors, zero warnings, and the same
three established notes as baseline: the long Phase 17C fixture path, installed
size (5.2 MB, including 4.2 MB under `libs`), and legacy scalar backend
`std::cout` symbols. Phase 17U introduced no new note.

## 42. Diff hygiene

Final audits check fixtures, generated wrappers, NAMESPACE/DESCRIPTION/Rd,
normal versus ignored-EOL statistics, and compiled/check artifacts.

## 43. Deviations and blockers

The Phase 17T rank denominator, stability rule, and short-chain minimum were
amended to match the authoritative current reference. No implementation
blocker remains.

## 44. Recommended next phase

> activate the validated internal Tier 1 convergence engine in public `mtblr_bed()` by adding `convergence` and `convergence_control`, integrating convergence-trace memory into public analytical warnings, exposing `fit$convergence` and optional `fit$convergence_traces`, issuing at most one aggregated advisory warning, and preserving `convergence="none"` and the existing single-chain numerical output exactly.

## 45. Readiness marker

PHASE 17U COMPLETE — INTERNAL MT BED TIER 1 CONVERGENCE ENGINE ACTIVE
