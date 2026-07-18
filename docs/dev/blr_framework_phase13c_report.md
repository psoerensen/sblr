# Unified BLR Framework Phase 13C Report

## 1. Executive summary

Packed-BED BayesR per-chain execution now runs behind a typed binding-neutral
context and returns a typed per-chain result. Phase 13A numerical behavior is
preserved.

## 2. Repository baseline

Branch `master`; starting/Phase 13B commit `aa1e24d`; initial tree and
`git diff --check` clean. Toolchain: R 4.4.1, Rtools44 GNU C++17/OpenMP.
The Phase 13B baseline was 5,041 passes, no failures or warnings, and four
opt-in skips.

## 3. Lexical dependency inventory

| Symbol/category | Type/dimensions | Ownership/lifetime | Typed destination |
|---|---|---|---|
| `G`, packed bytes/stride | `FastPackedBedMatrixBR`, m x packed-n | borrowed immutable, fit | `genotype` view |
| marker maps/order | vectors, m | borrowed immutable, fit | `marker_maps`, `marker_order` |
| `y_mat` | n x traits | borrowed immutable, fit | `phenotype` |
| `b_init` | traits x m | borrowed immutable, fit | `initial_effects` |
| B/E and priors | trait x trait matrices | borrowed immutable, fit | typed matrix references |
| `c`, `pi`, `alpha` | K vectors | borrowed immutable, fit | `components` |
| nub/nue/adjE/update flags | scalars | copied, chain call | context controls |
| nit/nburn/nthin/rebuild | scalars | copied, chain call | context controls |
| sweep/skip/candidate inputs | scalars | copied, chain call | `scheduler` |
| resolved seed, trait, chain | scalars | copied, logical chain | context metadata |
| effects/residual/state/variances/pi | native vectors/scalars | mutable, chain-owned | core locals |
| active/candidate/due state | native vectors | mutable, chain-owned | core locals |
| engine/uniform/normal/jitter/gamma | STL RNG objects | mutable, chain-owned | core locals |
| posterior/final/trace/CPO/timing | Armadillo/STL values | result-owned | typed result |
| progress message values | compact events | result-owned then adapter-read | `progress_events` |
| job index, task vector, failures | adapter containers | task/adapter-owned | excluded |
| R names/schema/conversion | Rcpp metadata | adapter-owned | excluded |

## 4. Files changed

The BayesR adapter constructs contexts and renders progress events; the core
header implements the typed callable; the new types header defines and validates
native contracts. Phase 13A/13B structural expectations were updated, Phase 13C
tests and benchmark were added, and the plan/matrix status was synchronized.

## 5. Component specification

`BedBayesRComponentSpec` borrows ordered scale, initial-probability and
Dirichlet vectors. Component zero is validated as null with scale zero; active
scales and priors remain positive. No reordering or execution-time normalization
was introduced.

## 6. Packed-genotype view

`BedBayesRPackedGenotypeView` borrows `FastPackedBedMatrixBR` plus immutable byte
pointer, packed size, marker/sample counts, bytes per marker and stride. The
fit owns storage for all calls; no chain copies bytes.

## 7. Typed per-chain context

`BedBayesRChainExecutionContext` contains borrowed genotype/statistical/component
inputs and copied model, iteration, scheduler, logical-seed and trait/chain
controls. It contains no R object, schema metadata, path, handle or worker ID.

## 8. Scheduler-contract reuse

`ScheduledSweepControl`, `NullSkipControl`, and `CandidateControl` are reused
inside `BedBayesRSchedulerControl`. BayesR retains probability-adaptive growth,
initial jitter, zero-as-every-iteration sweep semantics, burn-in control and no
neighbor wake-up.

## 9. Callable numerical core

`run_bed_bayesr_chain(context)` validates once, binds hot members locally,
constructs chain state and RNG, runs the sole unchanged MCMC/scheduler loop, and
returns the typed result. It performs no R decoding, dispatch, aggregation,
conversion or disk access.

## 10. Typed per-chain result

`BedBayesRChainExecutionResult` owns marker means/PIP/component summaries, final
effects/states, variance traces/finals, probability summaries, CPO, retained
count, timing/failure data and compact progress events. It owns no genotype,
RNG, scheduler internals or binding metadata.

## 11. RNG ownership

Each logical chain constructs one `std::mt19937`, uniform, normal and scheduler
jitter distribution after resolved seed mapping. Variable-shape gamma objects
remain local to Dirichlet updates. No worker or fit-persistent stochastic state
exists.

## 12. Scheduler preservation

Initialization jitter, due buckets, reservations/compaction, iteration-zero and
periodic full sweeps, and adaptive active -> candidates -> due traversal are
unchanged. Duplicate prevention, candidate expiry, null skipping, and skipped
marker no-RNG/no-genotype-work behavior are unchanged.

## 13. Component sampling

Log weights, stabilization, exponentiation, normalization, uniform draw,
component-order cumulative inverse-CDF, conditional non-null normal draw and
residual update retain their original order. No `std::discrete_distribution`
is used.

## 14. Genotype and I/O

BED/BIM/FAM validation and blocked decoding remain adapter-side and once per
fit. The core borrows packed marker-major storage and has no file/path API.

## 15. Progress boundary

The core no longer calls `Rcpp::Rcout`. At the same iteration conditions it
captures the same message inputs as binding-neutral events. The adapter emits
them after the OpenMP region in task order. Message contents and numerical
behavior are preserved; output is intentionally delayed until task completion.

## 16. Native adapter

The adapter retains public validation, BED decoding, task enumeration/static
OpenMP dispatch, seed resolution, context construction, progress rendering,
inline aggregation and inline R conversion.

## 17. Temporary aliases

No legacy `ChainResultBayesR` alias remains. The task-result vector directly
uses the typed result; inline aggregation matrices and conversion aliases remain
intentionally for Phase 13D.

## 18. Exact references

Phase 13A corrected raw references: 3/3 exact. Formatted references: 3/3 exact.

## 19. Reproducibility

`A; A`, `A; B; A`, normalized `1,2,2,1`, intervening BayesC/BayesRC and
different-chain-count sequences remain exact. The enabled fresh-process matrix
remains the authoritative fresh/reused check.

## 20. Component identities

Final and posterior-mean probabilities sum to one, marker component rows sum to
one, `dm = 1 - P(null)`, assignments remain in `[0,K-1]`, active variance is
`vb*c[k]`, and fixed probabilities retain existing normalization behavior.

## 21. Reductions and nonreductions

The documented BayesR-vs-BayesC nonreduction and dense/full-sweep skip-base
nonreduction remain unchanged; scheduler initialization jitter remains the first
relevant mechanism.

## 22. Performance, memory, and I/O

`tools/benchmarks/blr_phase13c_bed_bayesr.R` repeats 1x1, 2x1, 2x2, aggressive,
conservative and progress-enabled workloads with five repetitions. Timings are
too short for speed claims: observed medians were 0.00--0.02 seconds and
completed-fit RSS was 121.0--129.0 MiB. Completed-fit RSS is explicitly not peak
RSS; the tiny BED was 7 bytes and page-cache caveats are recorded. No genotype
or I/O path changed.

## 23. Public API and schema

Public arguments, defaults, routing, native export, generated wrappers,
`NAMESPACE`, raw schema, formatted schema, component ordering and actual `NULL`
fields are unchanged.

## 24. Protected backends

BayesRC, public/experimental/sparse BED BayesC, CSR backends, block-eigen and
multivariate production sources were not edited.

## 25. Tests

Phase 13A/13B structural tests were updated only for the typed symbol names.
Phase 13C adds architecture, exact-reference and progress tests. Final focused,
fresh-process and full-suite results were: Phase 13A--C focused 136 passes with
one opt-in skip; enabled fresh-process 61/61; full suite 5,077 passes with four
opt-in skips, no failures or warnings.

## 26. Deviations and blockers

Progress text is emitted after each task completes rather than from within the
worker while sampling. This removes unsafe binding access from OpenMP execution
without changing the message condition, fields, sampler control flow or output.
No numerical blocker was found.

## 27. Recommended next phase

Introduce one typed aggregate BayesR result, centralize cross-chain aggregation
and R conversion around typed chain results, remove temporary adapter aliases,
preserve all Phase 13A references, and complete packed-BED BayesR migration.

## 28. Readiness marker

PHASE 13C COMPLETE — PACKED-BED BAYESR TYPED CHAIN BOUNDARY ACTIVE
