# BLR Framework Phase 10C1 Report

## 1. Executive summary

The corrected scheduled ordinary-CSR BayesC execution block was mechanically extracted into one guarded implementation header without changing Phase 10B behavior. The production route still performs the same validation and preparation, executes the same scheduler and sampler statements through a lexical include, and constructs the same R result below the seam.

## 2. Repository baseline

- Branch: `master`
- Starting commit and Phase 10B commit: `e8fa6e4 Fix scheduled CSR chain RNG ownership`
- Initial status: clean
- Initial `git diff --check`: clean
- Toolchain: R 4.4.1 on Windows, Rtools44 GNU C++17/OpenMP toolchain
- Baseline native compile/load: passed
- Baseline focused Phase 10A and Phase 10B tests: 79 expectations passed, 0 failures, 0 warnings, 0 skips
- Phase 10B checkpoint full suite: 4,441 expectations passed

## 3. Extraction seam

In the Phase 10B source, the moved block began at line 534 with:

```cpp
const int ntasks = stblr_num_chain_tasks(nt, nchains);
```

This follows R decoding and validation, Armadillo input conversion, CSR loading and pre-scaling, friend-list preparation, marker ordering, scheduler-control preparation, and trait-chain seed preparation.

The moved block ended at line 1125 after cross-chain aggregation and immediately before line 1126:

```cpp
// Build named raw schema v1
```

The replacement include occupies the same lexical seam beneath the existing “Output storage” heading.

## 4. Files changed

- `src/st_cpg_omp_csr_scheduled.cpp`: replaced the 592-line execution block with one implementation-header include; retained helpers, validation, preparation, and R conversion.
- `src/blr_csr_scheduled_bayesc_core_impl.h`: added guarded implementation detail containing the mechanically moved execution body.
- `tests/testthat/test-blr-framework-phase10c1.R`: added permanent architecture, corrected-reference, reproducibility, dense-reduction, and protected-source tests.
- Phase 10A and 10B tests: replaced assumptions that execution statements remained directly in the production source with source-plus-implementation-header structural checks.
- Historical Phase 9 structural tests: removed the intentionally obsolete scheduled-production source MD5 entry; Phase 10C1 now supplies permanent scheduled structural protection.
- `docs/dev/blr_framework_phase10c1_report.md`: added this report.

## 5. Number of lines moved

Exactly 592 source lines were moved from the Phase 10B production source. Lines 7–598 of the new 600-line header are the moved body; the remaining lines are only the include guard, implementation-detail comment, spacing, and closing guard. The numerical and scheduling statements were not rewritten.

## 6. Extracted execution content

The header now contains:

- task-indexed output and failure allocation;
- OpenMP static trait-chain execution;
- task seed resolution and one `ScheduledChainRng` construction per logical trait-chain;
- chain-local effects, states, residuals, variances, scheduler state, buckets, candidates, and active-marker lists;
- the MCMC iteration loop and scheduled marker traversal;
- full sweeps, active markers, candidates, due buckets, skip growth/reset, and neighbor wake-up;
- marker inclusion/effect draws, residual changes, variance/global updates, and posterior accumulation;
- scheduler-adjacent timing, failure propagation, LD-swap values, chain summaries, and cross-chain aggregation.

## 7. RNG ownership preservation

The Phase 10B ownership model is unchanged: each logical trait-chain constructs exactly one `ScheduledChainRng` after final task-seed resolution. Its `std::mt19937`, normal distribution, and uniform distribution live for one chain execution, belong to no worker thread, are shared with no other chain, and do not survive a fit. Variable-parameter chi-square and gamma objects remain locally constructed at their existing draw sites. No `static` or `thread_local` stateful distribution was introduced.

## 8. Scheduler preservation

The extraction preserves the iteration-zero full sweep, `full_sweep_every`, adaptive null skipping, skip reset/growth, burn-in-only skipping, candidate entry/lifetime/expiry, friend ordering, neighbor wake-up threshold/limit, and the full-sweep versus active/candidate/due traversal order. Skipped markers still bypass the marker-update call and therefore consume no marker-update RNG.

## 9. Existing R conversion

The named `stblr_raw_v1` construction and its Rcpp helper lambdas remain in `src/st_cpg_omp_csr_scheduled.cpp` below the include. No R result construction was moved into the implementation header, and public names, types, dimensions, `NULL` values, schema, and routing remain unchanged.

## 10. Corrected frozen references

All three Phase 10B corrected configurations passed exactly:

- raw references: 3/3 exact;
- formatted references: 3/3 exact.

The Phase 10A defective fresh-process references remain untouched historical audit artifacts and are not used as Phase 10C1 expected production outputs.

## 11. Reproducibility

The corrected deterministic behavior remains exact for repeated `A; A`, intervening scheduled `A; B; A`, normalized `1,2,2,1` core ordering, an intervening unscheduled fit, a different-chain-count intervening fit, and explicit trait-chain seeds. The opt-in isolated-process Phase 10B check was also run after extraction: all three fresh-process configurations matched reused-process results exactly (39 Phase 10B expectations passed).

## 12. Dense reduction

Dense scheduled controls still do not reduce byte-for-byte to canonical unscheduled BayesC. This is the documented pre-existing scheduled-versus-unscheduled implementation difference and is unrelated to the extraction; no formula or scheduler behavior was changed to force equality.

## 13. Protected backends

MD5 audits confirm no changes to canonical CSR BayesC, BayesR, SBayesRC, fixed-prior BayesC, group BayesC, learned-annotation BayesC, block-eigen, individual scheduled backends, or multivariate CSR. BED scheduled sources, generated wrappers, `NAMESPACE`, public native signatures, and public schemas are unchanged.

## 14. Tests

- Baseline Phase 10A/10B focused suite: 79 passed.
- Phase 10C1 structural/reference/reproducibility suite: 44 passed.
- Final Phase 10A/10B/10C1 focused suite: 124 passed.
- Full `devtools::test('.')`: 4,486 passed.
- Failures: 0; warnings: 0; skips: 0.

## 15. Deviations and blockers

No numerical, scheduler, RNG-ownership, API, or schema deviation was introduced. The only additional test maintenance beyond Phase 10A/10B was removal of obsolete scheduled-source byte hashes in older protection tests; the new permanent structural tests cover the extracted source/header architecture. No blockers remain.

## 16. Recommended Phase 10C2

> replace the lexically dependent corrected scheduled execution include with an explicit typed scheduled execution context and callable core using the Phase 10A contracts and Phase 10B chain-owned RNG state, while retaining current R conversion until corrected references pass again.

## 17. Readiness marker

PHASE 10C1 COMPLETE — CORRECTED SCHEDULED CSR EXECUTION BLOCK EXTRACTED
