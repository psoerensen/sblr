# Unified BLR Framework Phase 10D report

## 1. Executive summary

Scheduled ordinary-CSR BayesC is canonicalized and stabilized. Its corrected
deterministic typed architecture is the sole production implementation.

## 2. Repository baseline

The repository began clean on branch `master` at Phase 10C3 commit `dfd4c9b`
(`Complete scheduled CSR BayesC migration`). `git diff --check` was clean.
R 4.4.1, Rtools44/GCC, Armadillo, and OpenMP compiled successfully. The baseline
full suite passed 4,579 expectations with no failures, warnings, or skips. The
Phase 10C3 benchmark baseline used seven tiny/2,000-marker configurations; its
representative medians ranged from 0.03 to 0.20 seconds and completed-fit RSS
was approximately 128--129 MiB.

## 3. Phase 10C3 structure inventory

| Item | Classification | Decision |
|---|---|---|
| Typed context/result, core, converter, native aggregation | stable architecture | retain permanently |
| `ScheduledChainRng` and scheduler contracts | corrected production contract | retain permanently |
| Phase 10B corrected fixtures | deterministic post-correction evidence | rename conceptually as canonical references |
| Phase 10A defective fixtures/diagnostic | original-defect evidence | historical audit artifact |
| Phase-labelled tests and reports | provenance and regression history | retain permanently |
| Phase 10C3 benchmark comparison wording | migration checkpoint language | replace with canonical baseline wording |
| “ready for canonicalization” documentation | obsolete status | remove now |
| Commented pre-framework source below the active entry | noncompiled historical source, not a fallback | defer; deletion is outside bounded canonicalization |

No live migration alias, unused typed field, duplicate converter, duplicate
aggregation, selector, or alternate RNG path was found. Therefore native files
were intentionally left byte-identical.

## 4. Files changed

The implementation plan and capability matrix mark the backend canonical. A
Phase 10D permanent architecture/reference test, canonical benchmark wrapper,
and this report were added. No native, public R, wrapper, namespace, fixture, or
protected-backend file changed.

## 5. Canonical execution path

Public scheduled `stblr_csr()` validation and alignment prepare CSR, friends,
seeds and scheduler controls; construct `CsrScheduledBayesCExecutionContext`;
call canonical `run_csr_scheduled_bayesc()`; return a typed result; call
`stblr_csr_scheduled_bayesc_result_to_raw()`; validate unchanged
`stblr_raw_v1`; and format through the canonical R formatter.

## 6. Cleanup and naming

Stable architecture names were retained. Canonical documentation replaces
“migration complete/ready” wording. Phase-labelled reports/tests and Phase 10A
defect evidence remain historical provenance. Reproducibility checks are
consolidated in the permanent Phase 10D suite without deleting earlier gates.

## 7. Typed execution context

The templated context borrows immutable operator/CSR, statistics, sample-size,
marker-order, friend, initial-value, prior and explicit-seed storage. It copies
dimensions, iteration/update/seed/output controls via typed scheduler controls.
Borrowed storage outlives the call. It contains no R objects, schema metadata,
Python objects, mutable chain state, or worker-owned RNG.

## 8. Numerical core and implementation header

The guarded implementation-detail header is included only by the scheduled
binding source. It holds one callable core, one active MCMC loop, one scheduler
transition implementation and one native aggregation path. It is binding
neutral and does not redefine compiler or Armadillo configuration.

## 9. Typed execution result

The result owns marker posterior and multichain summaries, final effects and
residuals, state, variance/pi/VLE/VLD traces and final values, retained counts,
task timings, and conversion dimensions. It exposes no RNG/distribution state,
mutable scheduler internals, public scheduler counters, or R metadata.

## 10. Native adapter

The adapter decodes and validates R input, aligns markers/traits, prepares CSR,
friends, seeds and controls, constructs the context, calls the core once, calls
the converter once, and translates exceptions. It contains no active sampler,
scheduler transition, RNG draw, numerical update, or aggregation formula.

## 11. Result converter

The sole named converter preserves field order, R types, dimensions, classes,
actual `NULL`, marker/trait/chain order, schema version, posterior/variance/pi,
diagnostic, LD-swap, selection, timing, failure and input fields. It performs no
sampling, scheduling, or numerical reaggregation.

## 12. RNG ownership

Each logical trait-chain constructs one fit-bounded `ScheduledChainRng` with
one `std::mt19937`, normal distribution and uniform distribution. No worker,
static, thread-local, shared, or cross-fit state exists. Variable-parameter
gamma/chi-square distributions remain at unchanged draw sites. Phase 10A's
cache diagnostic remains historical proof of the corrected defect.

## 13. Scheduler semantics

Iteration zero is a full sweep in canonical marker order. Non-full traversal is
active markers, candidates, then due buckets. Skip base/max/growth/reset and
burn-in-only behavior, candidate entry/lifetime/decrement/expiry, ordered and
bounded neighbor wake-up, and skipped-marker no-RNG behavior are unchanged and
chain-local.

## 14. Permanent regression fixtures

The three Phase 10B configurations permanently cover dense one-chain/one-core,
nontrivial two-chain/one-core skipping/candidates/wake-up, and two-chain/two-
core burn-in-only scheduling with explicit chain seeds. Raw and formatted
objects, LD-swap/output controls, types, dimensions and actual `NULL` are
protected. Fixture generation remains manual maintenance tooling. Phase 10A
pre-correction artifacts remain historical and are never canonical expected
trajectories.

## 15. Exact reference results

Corrected raw references: 3/3 exact. Corrected formatted references: 3/3 exact.

## 16. Reproducibility

Repeated `A; A`, intervening `A; B; A`, scheduled/unscheduled intervention,
different chain counts, normalized `1,2,2,1`, explicit seeds, and fresh versus
reused process comparisons are exact. Explicit trait-chain seeds are independent
of worker assignment.

## 17. Dense reduction

Dense scheduled execution remains nonidentical to canonical unscheduled BayesC
for the previously documented implementation-specific reason. Phase 10D did
not alter formulas or force a reduction.

## 18. Public API and schema

Arguments, defaults, native signature, routing, generated wrappers,
`NAMESPACE`, `stblr_raw_v1`, actual-`NULL` behavior, and formatted fit are
unchanged.

## 19. Unsupported behavior

Scheduled ordinary-CSR BayesR, BayesRC and SBayesRC remain unavailable.
Existing trait/shared-`ww`, `keep_chains`, and public scheduler-counter
limitations remain unchanged.

## 20. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior, group, and learned-annotation
sources remain hash-identical. Block-eigen, packed-BED scheduled,
individual-level scheduled, multivariate, generated-wrapper and namespace
protections pass.

## 21. Performance and memory baseline

Run `Rscript tools/benchmarks/blr_phase10d_scheduled_csr.R`. The canonical
baseline reuses the warm-up and five-repetition tiny and 2,000-marker Phase
10C3 workloads for exact comparability. It records all controls, individual
times, mean/median/min/max/range, R/toolchain information, and whole-process RSS
after completed fits. RSS is not peak memory; interval sampling is unavailable.
Phase 10D medians were 0.01 s (tiny), 0.03 s (dense 1/1), 0.04 s (dense
2/1), 0.03 s (dense 2/2), 0.02 s (aggressive skip), 0.03 s (conservative
skip), and below the 0.01 s timer resolution (wake-up). Completed-fit RSS was
127.6--129.1 MiB for representative workloads and 137.9 MiB for tiny. Against
the Phase 10C3 medians (0.03--0.20 s representative) and approximately
128--129 MiB RSS, there is no material regression. Short Windows timings are
not interpreted as a speed improvement.

## 22. Tests

Phase 10A, 10B, 10C1, 10C2, 10C3 and new 10D passed 41, 39, 44, 62, 31
and 44 expectations respectively, with no failures, warnings or skips. The
opt-in fresh-process matrix passed 39 expectations. Native compilation and
package loading succeeded; the final full suite passed 4,623 expectations with
zero failures, warnings or skips.

## 23. Deviations and blockers

No numerical, scheduler, RNG, API or schema deviation is intended. No blocker
was identified during inventory.

## 24. Recommended next phase

> begin a contract, reference, scheduler, and RNG-ownership audit of the individual-level and packed-BED scheduled backends, comparing their execution models and identifying whether shared typed scheduled infrastructure is appropriate before migration.

## 25. Readiness marker

PHASE 10D COMPLETE — SCHEDULED CSR BAYESC CANONICALIZED AND STABILIZED
