# Unified BLR Framework Phase 13D Report

## 1. Executive summary

Packed-BED BayesR migration is complete with one typed per-chain core, one
binding-neutral aggregation path, and one named R converter. Deterministic Phase
13A behavior is unchanged.

## 2. Repository baseline

Branch `master`; starting/Phase 13C commit `51c1652`; initial tree and
`git diff --check` clean. R 4.4.1, Rtools44 GNU C++17/OpenMP. Phase 13C compiled
and passed 5,077 tests with four opt-in skips. Its tiny benchmark medians were
0.00--0.02 seconds with completed-fit RSS of about 121--129 MiB.

## 3. Aggregation/conversion inventory

| Item | Phase 13C location | Phase 13D classification/location |
|---|---|---|
| task vector and trait/chain mapping | adapter | retained in dispatch; mapping consumed by aggregation |
| marker means/PIP/component probability | inline | native aggregation |
| final effects/states | inline | native aggregation, arithmetic mean |
| trace/final variance and pi | inline | native aggregation, arithmetic mean |
| SD/min/max | inline | native aggregation; sample SD (`n-1`) |
| CPO/retained samples | inline | native aggregation, arithmetic mean |
| timing | inline | native mean and maximum |
| failures | adapter | adapter exception handling plus aggregate count |
| optional `wy/r` scores | inline | native post-aggregation diagnostic fill; genotype-dependent, not chain aggregation |
| orientations/component labels | inline | named converter |
| R matrices/lists/classes/NULL/schema | inline | named converter |

## 4. Files changed

The types header gains aggregation/result contracts; a guarded aggregation
header owns formulas; the adapter calls aggregation and conversion once. Phase
13B/C structural assertions were updated, and Phase 13D tests, benchmark, report,
plan and matrix status were added.

## 5. Final execution path

`stblr_bed(method="bayesr")` -> R validation -> native BED validation/decoding
-> static OpenMP task dispatch -> typed contexts -> `run_bed_bayesr_chain()` ->
typed chain vector -> `aggregate_bed_bayesr_results()` -> optional native
`wy/r` fill -> `stblr_bed_bayesr_result_to_raw()` -> unchanged raw/fit.

## 6. Typed aggregate result

`BedBayesRExecutionResult` owns dimensions; marker means/PIP/component means;
final effects/states; SD/min/max; optional scores; all variance traces/finals;
per-trait component probabilities; final/mean pi; CPO; retained counts; timing;
and failure count. It owns no RNG, scheduler, genotype, path or R metadata.

## 7. Native aggregation

The singular callable preserves `job=chain*traits+trait`, arithmetic means,
sample SD with `nchains-1`, elementwise minima/maxima, component order, trace and
pi averaging, CPO/retained/timing means, and timing maxima. It consumes no RNG,
genotype, scheduler, file or R object.

## 8. Task dispatch

Trait-chain enumeration, `schedule(static)`, seed resolution, context
construction, exception collection and one result per task remain adapter-owned.

## 9. Progress handling

The core captures binding-neutral events at unchanged progress points. The
adapter renders them after task completion in task order. Aggregation and the
converter neither consume nor emit progress.

## 10. Named result converter

`stblr_bed_bayesr_result_to_raw()` consumes the aggregate result and binding
metadata. It performs orientation, naming, schema/list/class construction and
preserves all actual `NULL` fields. It performs no numerical aggregation.

## 11. Temporary aliases removed

Inline aggregate matrices, duplicate chain-index loops and inline summary
formulas were removed. The typed task vector remains. Binding metadata and small
orientation lambdas are retained solely inside the converter. Optional score
calculation remains a documented native post-aggregation diagnostic seam.

## 12. Final native adapter

The adapter is bounded to validation/decoding, preparation, dispatch, progress,
one aggregation call, optional requested genotype scores, one converter call and
exception translation. It contains no MCMC, scheduler, RNG, posterior or CPO
aggregation formula.

## 13. Per-chain core

There remains one binding-neutral core, one MCMC loop, one scheduler and one
logical-chain RNG construction path.

## 14. Component semantics

Component zero remains null/scale zero; active variance remains `vb*c[k]`.
Scale, pi, prior, state and output ordering and custom inverse-CDF sampling are
unchanged.

## 15. RNG ownership

One `std::mt19937`, uniform, normal and jitter distribution belongs to each
logical chain. Variable-shape Gamma distributions remain draw-local. No worker
or fit-persistent stochastic state exists.

## 16. Scheduler preservation

Initialization jitter, full sweeps, active -> candidates -> due traversal,
probability-adaptive skipping, candidate lifetime, compaction, and skipped-marker
no-RNG/no-genotype-work behavior are unchanged.

## 17. Genotype and I/O

Fit-owned packed data remain borrowed immutable. Decoding stays adapter-side;
neither core nor aggregation accesses files. No per-chain copy or extra decode
was introduced.

## 18. Exact references

Raw references: 3/3 exact. Formatted references: 3/3 exact.

## 19. Reproducibility

Repeated, intervening-fit, normalized core-order, fresh/reused, different-chain
count and worker-assignment checks remain exact.

## 20. Component identities

Final/mean pi and marker probability sums, `dm=1-P(null)`, state range,
`vb*c[k]`, and fixed-probability behavior remain valid.

## 21. Nonreductions

BayesR-vs-BayesC and full-sweep skip-base nonreductions retain their Phase 13A
first differences; no equality was forced.

## 22. Public API and schema

Arguments, defaults, routing, native export, wrappers, `NAMESPACE`, raw field
order/types/dimensions/classes/NULL and formatted fit are unchanged.

## 23. Protected backends

BayesRC, every BED BayesC route, CSR backends, block-eigen and multivariate
production sources were not modified.

## 24. Performance, memory, and I/O

The Phase 13D script repeats six four-component workloads five times and records
completed-fit RSS, BED size and cache caveats. Tiny timings do not support speed
claims: medians were 0.00--0.02 seconds and completed-fit RSS was 121--134 MiB,
consistent with Phase 13C noise. Peak RSS and moderate/larger workloads remain
opt-in resource runs.

## 25. Tests

Phase 13D adds permanent architecture, aggregation, references, identities and
reproducibility protections. Phase 13A--D focused tests passed 183 checks with
one opt-in skip; enabled fresh-process tests passed 61/61; the full suite passed
5,124 checks with four opt-in skips and no failures or warnings.

## 26. Deviations and blockers

Optional `wy/r` diagnostics require prepared genotype access and therefore are
filled after aggregation but before conversion; they are not cross-chain summary
formulas. No numerical blocker was found.

## 27. Recommended Phase 13E

Canonicalize and stabilize public packed-BED BayesR, retain the Phase 13A
fixtures permanently, remove remaining migration-era wording, establish the
Phase 13D benchmark as the canonical baseline, and leave BayesRC and all BayesC
routes unchanged.

## 28. Readiness marker

PHASE 13D COMPLETE — PACKED-BED BAYESR MIGRATED WITH DETERMINISTIC BEHAVIOR PRESERVED
