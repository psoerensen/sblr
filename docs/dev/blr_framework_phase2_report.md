# Unified BLR Framework: Phase 2 Report

## 1. Executive summary

Phase 2 migrated the existing unscheduled CSR BayesC sampler behind the Phase
1 typed specification, borrowed CSR view, and typed result vocabulary. The
optimized sampler logic was moved without redesigning its mathematics or hot
loop. The existing native entry still constructs the exact `stblr_raw_v1`
object, and the public R interface and formatter are unchanged.

## 2. Repository baseline

- Branch: `master`.
- Starting commit: `7df2d1ddec1c5f45be03e399dc7afb419afd56d8`.
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1, Rtools44, GNU C++17, OpenMP. The initial compiled load
  completed successfully in 275.1 seconds.
- Full baseline: 3,283 passed, zero failures, warnings, or skips.
- Focused baseline: 1,222 passed, zero failures, warnings, or skips.
- Phase 1 moderate baseline used 2,000 markers, two independent traits, 80
  sampling iterations, 20 burn-in iterations, and one or two chains/cores.
  Its three-run medians were approximately 0.21, 0.34, and 0.25 seconds for
  minimal output (1/1, 2/1, and 2/2 chains/cores), and 0.17, 0.34, and 0.27
  seconds for ordinary output. Sampled peak RSS ranged from 144.1 to 149.5 MB.

The exact baseline commands were:

```text
Rscript -e "Rcpp::compileAttributes('.')"
Rscript -e "pkgload::load_all('.', compile = TRUE)"
Rscript -e "devtools::test('.')"
Rscript tools/benchmarks/blr_phase1_csr_bayesc.R --fixture=tiny --nrep=3
Rscript tools/benchmarks/blr_phase1_csr_bayesc.R --fixture=moderate --nrep=3
```

## 3. Extraction seam

The extraction seam is inside the unchanged `stblr_cpg_omp_csr()` native
entry after existing Rcpp decoding, validation, native conversion, shared CSR
construction, scaling checks, and marker-order construction, but before the
trait-chain OpenMP loop.

At this seam the binding constructs a Phase 1 `ResolvedSpec`, borrowed
immutable `CsrBayesCDataView`, typed priors, initial state, controls, output
request, LD-friend view, and marker-order view. `run_csr_bayesc()` performs
trait-chain execution and aggregation and returns a typed result. The adapter
then performs only R-specific matrix orientation, naming, actual-`NULL`
handling, class/version metadata, and raw-list construction.

The block-eigen instantiation remains on its original generic path. It is not
routed through the CSR data view and its mathematics and wrapper are unchanged.

## 4. Files changed

- `src/st_cpg_omp_csr.cpp`: adds the typed adapter and routes only ordinary
  unscheduled CSR BayesC through the shared implementation. Its exported
  signature is unchanged.
- `src/blr_csr_bayesc_types.h`: defines borrowed data, priors, initial state,
  controls, execution input, chain result, and aggregate typed result.
- `src/blr_csr_bayesc_core.h`: contains the moved binding-neutral sampler and
  aggregation implementation. It is an implementation header so it inherits
  the package's single established Armadillo configuration; a separate native
  translation unit would instantiate an incompatible alternate-RNG ABI.
- `tests/testthat/fixtures/blr-phase2-reference.R`: defines the compact
  deterministic fixture, metadata, normalizer, and raw/formatted runners.
- `tests/testthat/fixtures/blr-phase2-reference-hashes.R`: stores compact MD5
  identities frozen before the sampler source was edited.
- `tests/testthat/test-blr-framework-phase2.R`: adds exact comparisons,
  reproducibility, schema, ownership, neutrality, and allocation checks.
- `tools/benchmarks/blr_phase2_csr_bayesc.R`: adds warm-up-excluded repeated
  timings and optional sampled peak RSS for the Phase 1 workloads.
- `docs/dev/blr_framework_phase2_report.md`: this report.

`R/RcppExports.R`, `src/RcppExports.cpp`, and `NAMESPACE` are unchanged after
regeneration. No protected backend source outside `src/st_cpg_omp_csr.cpp`
changed.

## 5. Typed execution input

`CsrBayesCExecutionInput` contains:

- the explicit Phase 1 `ResolvedSpec`;
- `CsrBayesCDataView`, which borrows immutable CSR row pointers, column
  indices, float values, marker diagonals, aligned scores, trait sums of
  squares, and sample sizes;
- `CsrBayesCPriors` and `CsrBayesCInitialState`;
- `CsrBayesCControls`, preserving iterations, burn-in, thinning, chains,
  cores, seed and explicit chain seeds, update flags, priors, residual rebuild,
  LD-swap, `selection_s`, and output controls;
- borrowed LD-friend and marker-order views.

Validation runs before OpenMP execution and rejects missing/mismatched views,
dimensions, MCMC controls, chain-seed length, unsupported specifications, and
invalid LD-swap or `selection_s` controls.

## 6. Shared implementation

The existing chain seed construction, static trait-chain task mapping, marker
traversal, BayesC state probability and effect draw, residual updates, variance
and pi updates, LD swap, `selection_s`, retained-sample accumulation,
diagnostics, and chain aggregation were moved operation-for-operation.

No arithmetic order, loop order, marker traversal order, RNG invocation order,
distribution ownership, residual rebuild timing, OpenMP schedule, task mapping,
or aggregation order changed. The helper functions retain their original
inline declarations. Marker-loop source guards confirm that no explicit heap
allocation, resize, vector construction, string work, logging, virtual call,
or runtime dispatch was introduced inside either marker traversal branch.

## 7. Chain and RNG ownership

Each trait-chain task owns effects, residual, inclusion state, variance and pi
state, its `std::mt19937`, stateful distribution use, posterior accumulators,
diagnostics, and scratch storage. Seeds are still derived by the existing
`stblr_trait_seed`, `stblr_chain_seed`, and
`stblr_seed_with_chain_base` helpers.

CSR values, row pointers, column indices, marker diagonals, fixed statistics,
sample sizes, priors, LD friends, and marker order are shared immutable views.
The `CsrOperator` owner in the binding outlives all OpenMP work. No chain copies
CSR storage.

## 8. Typed result integration

`CsrBayesCResult` extends the Phase 1 `BlrResult` vocabulary and contains the
current BayesC marker summaries and final state, variance/pi/selection traces,
final variances and pi, retained-sample and diagnostic vectors, and optional
per-chain results. Native storage remains trait-major for the existing hot
path; only the thin adapter converts to canonical markers-by-traits and
samples-by-traits R matrices.

The adapter populates every current `stblr_raw_v1` BayesC field: marker
`bm/dm/wy/r/b/state`, traces, covariance/final variance fields, final/mean pi,
metadata, inputs, diagnostics, optional chains, `selection_s`, and LD-swap.
No later-model payload was added.

## 9. Public API and schema

- Public `stblr_csr()` arguments, defaults, routing, and formatted fields are
  unchanged.
- The exported native `stblr_cpg_omp_csr()` signature is unchanged.
- `stblr_raw_v1`, including actual `NULL` and empty-list behavior, is unchanged.
- `NAMESPACE` and generated R/native wrappers are unchanged.
- No resolved specification is attached to a public raw or fit object.
- Scheduled CSR, BayesR, BayesRC/SBayesRC, annotation/group/prior, BED,
  multivariate, and block-eigen sources were not modified.

## 10. Reference comparison

Seven compact configurations were frozen before editing the sampler source:
one trait/one chain/one core; one trait/two chains with one and two cores;
multiple independent traits; explicit seeds 401 and 402; `keep_chains = TRUE`;
and fixed `selection_s = -0.5`. The common seed was 31, with 8 iterations, 2
burn-in iterations, thinning 1, and LD swap disabled.

All seven complete normalized raw hashes and all seven complete normalized fit
hashes match their pre-refactor values exactly. Normalization changes only
temporary CSR paths and measured runtime seconds. Values, types, dimensions,
names, class, schema/version, actual-`NULL` behavior, chain ordering, and all
other stable fields are included in the identities.

## 11. Reproducibility

The Phase 2 tests passed exact repeated-seed comparisons, one-core/two-core
comparisons, the sequence one/two/two/one cores, an intervening multiple-trait
fit, explicit chain-seed mapping, multiple-trait column order, chain retention,
fixed and disabled `selection_s`, and disabled LD-swap behavior. Results depend
only on inputs and seeds in all protected configurations.

## 12. Performance and memory

Post-migration commands were:

```text
Rscript tools/benchmarks/blr_phase2_csr_bayesc.R --fixture=tiny --nrep=3
Rscript tools/benchmarks/blr_phase2_csr_bayesc.R --fixture=moderate --nrep=5 --peak=FALSE
Rscript tools/benchmarks/blr_phase2_csr_bayesc.R --fixture=moderate --nrep=1
```

Each configuration has one excluded warm-up. The timing-only five-run moderate
medians were 0.15, 0.24, and 0.14 seconds for minimal output (1/1, 2/1, and 2/2
chains/cores), and 0.14, 0.25, and 0.14 seconds for ordinary output. Means were
0.248, 0.242, 0.140, 0.138, 0.242, and 0.136 seconds; the 1/1 minimal mean
contains a 0.67-second outlier. These results show no material regression
against the Phase 1 medians, but are not evidence of a performance improvement.

The current moderate peak-RSS run measured 149.4, 149.8, and 150.2 MB for
minimal output and 150.0, 149.5, and 150.1 MB for ordinary output. Phase 1 used
the same whole-child-process method and measured 144.1--149.5 MB. Package load
dominates this measure; the observed maximum increase is below one percent of
the approximately 150 MB process baseline and is not a meaningful memory
regression. Tiny peak RSS was 144.3--147.4 MB and sub-centisecond timings were
below the timer's useful resolution after warm-up.

Peak RSS is the maximum child RSS sampled every 10 ms using optional installed
`processx`/`ps`; no dependency was added. Timed public calls include validation,
CSR loading, typed conversion, execution, and result conversion. The conversion
boundary is not separately instrumented, so initialization/conversion time is
reported only as part of total elapsed time. Short Windows timings and RSS are
variable and are not used to claim a speedup.

## 13. Test results

- Full baseline: 3,283 passed; 0 failed, 0 warned, 0 skipped.
- Baseline focused set: 1,222 passed; 0 failed, 0 warned, 0 skipped.
- Phase 1 post-change: 102 passed; 0 failed, 0 warned, 0 skipped.
- New Phase 2 file: 63 passed; 0 failed, 0 warned, 0 skipped.
- Combined focused regression set: 1,285 passed; 0 failed, 0 warned, 0 skipped.
- Full post-change suite: 3,346 passed; 0 failed, 0 warned, 0 skipped.

The final clean native load completed successfully in 149.9 seconds. Compiler
output contained only pre-existing unused-function warnings; testthat emitted
its package-built-under-R-4.4.3 startup warning, while the test suite itself
reported `WARN 0`.

## 14. Deviations and blockers

The shared implementation is a binding-neutral implementation header rather
than a separate `.cpp`. A standalone Armadillo translation unit instantiated a
different alternate-RNG/configuration ABI from the existing RcppArmadillo
translation units and caused an access violation in an untouched scheduled
backend. Including the neutral header after the established package Armadillo
configuration removes the duplicate ABI, keeps the shared source free of R or
Python types, and restores the scheduled and full suites. This is a resolved
implementation deviation, not a remaining blocker.

Peak RSS initially required an unsandboxed child-process permission and the
benchmark now drains child output while sampling. Conversion time is not
separately measurable without adding intrusive instrumentation. No acceptance
blocker remains.

## 15. Temporary comparison code

The compact reference definitions and pre-refactor hashes are test-only
comparison scaffolding. No second CSR sampler or legacy public/native entry was
retained. The existing generic loop remains only because block-eigen still uses
that operational path. Phase 3 may remove the frozen comparison scaffolding
after stabilization; it must not remove the block-eigen path as unrelated
cleanup.

## 16. Recommended Phase 3 boundary

Switch all unscheduled CSR BayesC internal execution to the migrated canonical
path, remove any temporary comparison mechanism, stabilize the public route,
and then remove superseded duplicate implementation code. Because the native
entry already uses the canonical path, Phase 3 should focus on stabilization
and removal of temporary reference scaffolding.

## 17. Readiness marker

PHASE 2 COMPLETE — CSR BAYESC MIGRATED WITH BEHAVIOR PRESERVED
