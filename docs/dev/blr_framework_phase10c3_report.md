# Unified BLR Framework Phase 10C3 report

## 1. Executive summary

Scheduled ordinary-CSR BayesC migration is complete with one typed numerical
core and one named binding-layer result converter. Phase 10B deterministic RNG
ownership and scheduler behavior are preserved.

## 2. Repository baseline

The clean baseline was branch `master` at `d9e96dc` (`Activate typed scheduled
CSR BayesC execution boundary`), also the Phase 10C2 commit. Initial status and
`git diff --check` were clean. R 4.4.x with Rtools44/GCC and OpenMP compiled the
baseline. Phase 10A--10C2 focused tests passed 186 expectations. The Phase 10C2
benchmark baseline used the Phase 10B workloads; medians were 0.01 s (tiny),
0.03 s (dense 1 chain/1 core), 0.19 s (dense 2/1), 0.10 s (dense 2/2), 0.03 s
(aggressive skip), 0.09 s (conservative skip), and 0.03 s (neighbor wake-up).
Completed-fit whole-process RSS was approximately 128.0--137.6 MiB.

## 3. Phase 10C2 conversion inventory

| Item | Phase 10C2 state | Classification |
|---|---|---|
| Typed result fields | Complete numerical output | retain permanently |
| 23 result aliases | Former lexical names for marker, trace, variance, state, timing, and failure data | remove now |
| Marker/trace/diagonal lambdas | Inline binding conversion | centralize now |
| Raw list assembly | Inline after core call | centralize now |
| `wy`, sample sizes, iteration controls | Binding-only schema/input data | binding-only metadata |
| Native task/chain aggregation | Inside callable core | retain permanently |
| Inline-converter comments/assertions | Migration-era boundary | rename for stable use |

No unused typed-result field or duplicate numerical aggregation path was found.

## 4. Files changed

`src/st_cpg_omp_csr_scheduled.cpp` gains the named converter and removes the
compatibility aliases and inline duplicate. Phase 10C2 and new Phase 10C3 tests
protect the permanent boundary and corrected references. The Phase 10C3
benchmark preserves comparable workloads. The implementation plan and
capability matrix record migration closure. This report records evidence and
decisions.

## 5. Final execution path

Public scheduled `stblr_csr()` validation and alignment prepare CSR, friends,
seeds, and scheduler controls; the adapter constructs
`CsrScheduledBayesCExecutionContext`, calls `run_csr_scheduled_bayesc()`, passes
its typed result to the named converter, validates `stblr_raw_v1`, and uses the
unchanged canonical formatter.

## 6. Named result converter

`stblr_csr_scheduled_bayesc_result_to_raw()` consumes the typed result and a
small binding metadata view. It preserves list order, storage types,
dimensions, classes, actual `NULL`, traces, variance, pi, diagnostics,
selection, and schema version. It performs no sampling or scheduler work.

## 7. Temporary aliases removed

Aliases for `bm`, `dm`, their chain summaries, `b`, `r`, state, all six traces,
final variances and pi, retained counts, and task seconds were removed. The
binding metadata view is retained because `wy`, sample sizes, and iteration
controls are intentionally outside the numerical result.

## 8. Native adapter boundary

The adapter decodes and validates R inputs, aligns data, prepares CSR/friends,
seeds and controls, constructs the typed context, calls the core once, and
calls the converter once. It contains no MCMC, marker update, scheduler
transition, RNG draw, residual/variance update, or posterior formula.

## 9. Numerical aggregation

Task and chain aggregation occur once inside the typed core. The converter only
transposes native outputs and forms binding-level list objects and timing/pi
summaries exactly as before; it does not reaggregate chain effects.

## 10. Numerical core

One callable core, one MCMC loop, and one scheduler implementation remain in
the guarded implementation header. The core and typed headers remain free of
Rcpp and Python binding types.

## 11. RNG ownership

Each logical trait-chain constructs one `ScheduledChainRng` containing its own
`std::mt19937`, normal distribution, and uniform distribution. Lifetime is one
chain execution. No worker-owned, static, thread-local, or fit-persistent
distribution state exists.

## 12. Scheduler preservation

Iteration-zero full sweeps, active/candidate/due traversal, skip growth/reset,
candidate lifetime/expiry, ordered neighbor wake-up, and skipped-marker no-RNG
behavior remain in the numerical core with OpenMP static task scheduling.

## 13. Corrected frozen references

All three Phase 10B corrected raw references and all three formatted references
match exactly. Phase 10A defective references remain untouched historical
artifacts.

## 14. Reproducibility

`A; A`, `A; B; A`, normalized `1,2,2,1`, intervening unscheduled and scheduled
fits, different chain counts, and explicit seeds remain exact. Fresh-process
artifacts retain the declared fresh-process provenance and match reused-process
sampled output.

## 15. Dense reduction

Dense scheduled execution remains non-identical to canonical unscheduled
BayesC, as before Phase 10C3. Converter extraction cannot affect this known
scheduler/implementation distinction.

## 16. Public API and schema

Arguments, defaults, native signature, routing, generated wrappers,
`NAMESPACE`, `stblr_raw_v1`, formatted fields, and actual-`NULL` behavior are
unchanged.

## 17. Unsupported behavior

Scheduled ordinary-CSR BayesR, BayesRC, and SBayesRC remain unsupported, as do
the existing unsupported trait/shared-`ww`, chain-retention, and scheduler
counter cases.

## 18. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior, group, and learned-annotation
sources are hash-protected. Block-eigen, BED/individual scheduled, and
multivariate sources remain unchanged.

## 19. Performance and memory

`Rscript tools/benchmarks/blr_phase10c3_scheduled_csr.R` retains the
Phase 10B--10C2 2,000-marker workloads, warm-up, five repetitions, and all
scheduler modes. Phase 10C3 medians were 0.02 s (tiny), 0.03 s (dense 1/1),
0.07 s (dense 2/1), 0.03 s (dense 2/2), 0.01 s (aggressive), 0.03 s
(conservative), and 0.01 s (wake-up). Completed-fit RSS after representative
fits was 127.9--129.3 MiB (tiny 137.5 MiB). This is whole-process post-fit RSS,
not peak memory. Against the Phase 10C2 medians and 128.0--137.6 MiB RSS, no
material regression is visible; noisy short Windows timings are not treated as
speed improvements.

## 20. Tests

Phase 10A, 10B, 10C1, 10C2, and 10C3 passed 41, 39, 44, 62, and 31
expectations respectively. The opt-in fresh-process Phase 10B matrix passed 39.
The full suite passed 4,579 expectations with zero failures, warnings, or
skips. Native compilation and package loading succeeded.

## 21. Deviations and blockers

No model, scheduler, RNG, API, or schema deviation occurred. An initial
sandboxed fresh-process check could not launch a child process; the required
approved unsandboxed rerun passed all 39 expectations. There are no blockers.

## 22. Recommended Phase 10D

> canonicalize and stabilize scheduled ordinary-CSR BayesC, remove remaining migration-era wording, retain the Phase 10B corrected fixtures permanently, establish the Phase 10C3 benchmark as the canonical baseline, and then begin the next noncanonical backend audit.

## 23. Readiness marker

PHASE 10C3 COMPLETE — SCHEDULED CSR BAYESC MIGRATED WITH DETERMINISTIC BEHAVIOR PRESERVED
