# Unified BLR Framework Phase 15B Report

## 1. Executive summary

Only proven nonnumerical packed-BED common infrastructure was consolidated and
canonicalized: task indexing, logical-chain seed resolution, and the compatible
BayesC/R packed-genotype view. Numerical implementations remain model-specific.

## 2. Repository baseline

The clean `master` baseline was `3b5fbe9` (Phase 15A). Its canonical references
were 9/9 raw and 9/9 formatted exact, the enabled fresh-process selection passed
157 assertions, and the full suite passed 5,497. R 4.4.1/Rtools44, sblr 0.1.0,
and the Phase 11D/13E/14E model-specific benchmarks form the baseline.

## 3. Phase 15A decisions

Task/seed logic was exact-shareable; BayesC/R views were representation-equivalent;
BayesRC was ownership-equivalent only. Complete dispatch, timing, progress,
diagnostics, contexts/results, numerical cores, aggregators, converters,
schedulers, and probability policies were retained model-specific.

## 4. Files changed

`blr_bed_family_types.h` adds the narrow production vocabulary. BayesC/R type
headers adopt the common view; the BayesRC type header imports only task/seed
helpers. Three canonical adapters adopt task/seed helpers. Shared test and
benchmark helpers, Phase 15B tests/tooling, plan/matrix, and this report are added
or updated. Numerical core, aggregate implementation, converter body, wrapper,
and schema files are not changed.

## 5. Shared task-index helper

`BedFamilyTaskIndex { int job, trait, chain; }` and
`make_bed_family_task_index(int job,int trait_count)` implement exactly
`trait=job%trait_count`, `chain=job/trait_count`, and preserve `job` as the result
slot. All three adapters use it. Representative 1x1, 2x1, 1x3, and 3x4 mappings
are permanently tested.

## 6. Shared seed resolver

`resolve_bed_family_logical_chain_seed(int seed,int trait,int chain)` preserves
signed expression evaluation followed by `static_cast<unsigned int>` and return
as `uint64_t`: `seed+1000003*(trait+1)+9176*(chain+1)`. All adapters resolve it at
context construction; worker identity is absent. Representative zero, ordinary,
large supported seed, trait, and chain inputs match the old formula.

## 7. Shared BayesC/R genotype vocabulary

`BedPackedGenotypeView<PackedGenotype>` owns no bytes and contains the immutable
storage reference, byte pointer, packed size, marker/sample counts, bytes per
marker, and stride in the original aggregate-initialization order. BayesC uses it
directly; BayesR retains the readable `BedBayesRPackedGenotypeView` alias. Storage
remains fit-owned and no chain copy is introduced.

## 8. BayesRC genotype decision

`BedBayesRCPackedGenotypeView` remains separate with storage, marker/sample counts,
and bytes per marker only. Unused pointer/size/stride fields were not added. It
shares the immutable fit-owned/no-copy/no-reread ownership contract.

## 9. Failure-envelope decision

Final decision: **model-specific production failure payloads retained**. Although
all chain results have `failed/error`, exception timing, rendering, aggregate
storage, and public presentation differ. Shared structural test vocabulary covers
the common contract without introducing ambiguous partial production sharing.

## 10. OpenMP dispatch

Three separate static-schedule loops remain. Context construction, typed result
assignment, exceptions, and progress rendering remain in their adapters. Only
task/seed helpers are shared; no executor, callback, polymorphism, or loop movement
was introduced. Task enumeration and result order are unchanged.

## 11. Timing and progress

BayesC timing/progress, BayesR per-chain timing and typed progress rendering, and
BayesRC behavior remain unchanged and model-specific. No new progress event or
common timing production type was introduced.

## 12. Optional diagnostics

BayesC retains its existing diagnostic boundary; BayesR/BayesRC retain their
post-aggregation fills. Genotype ownership is not moved into aggregate types and
actual R `NULL` behavior is unchanged.

## 13. Converter and orientation helpers

No converter or binding helper was shared: the small apparent duplication is
order/schema sensitive. Common two-dimensional conventions remain documented;
binary, component, stick, and annotation orientations remain separate.

## 14. Test helper consolidation

The Phase 12A helper now provides zero-safe occurrence counting, forbidden-token
checks, and protected-source hash checks. Model-specific semantic assertions are
retained in their permanent tests.

## 15. Reference helper consolidation

The shared test helper adds a fresh-process launcher and recursive first-difference
reporting for type, dimensions, names/dimnames, class, path, index, expected, and
observed values. Capture and normalization remain model-specific; no normalization
was broadened and fixtures were not regenerated.

## 16. Benchmark convention

The common reporting helper standardizes versions, dimensions, MCMC controls,
individual/summary timings, completed-fit RSS, BED size, and notices. Scheduler,
mixture, and annotation controls remain model-specific. It states no cross-model
ranking, model-specific workloads, completed RSS is not peak, page-cache effects,
and tiny timings are regression signals.

## 17. Model-specific architecture protection

BayesC retains binary inclusion and adaptive null/candidate scheduling. BayesR
retains ordered global mixtures, inverse-CDF sampling, and adaptive scheduling.
BayesRC retains marker annotation probit sticks, latent/sequential-alpha updates,
and full sweeps. Contexts, chain/aggregate results, cores, aggregators, converters,
and schemas remain distinct.

## 18. Exact references

BayesC, BayesR, and BayesRC each retain 3/3 raw and 3/3 formatted exact fixtures:
9/9 raw and 9/9 formatted in total.

## 19. Reproducibility

Same-process, fresh-process, core-order, chain-count, worker-assignment, and
intervening-fit checks remain exact for all three families.

## 20. Fixed-alpha reduction

The intercept-only fixed-alpha BayesRC to matched fixed-pi BayesR reduction remains
exact; general model nonreductions remain protected.

## 21. Public API and schema

Arguments, defaults, routes, native signatures, raw/formatted field order and
types, actual `NULL`, generated wrappers, and `NAMESPACE` are unchanged.

## 22. Protected backends

Experimental/sparse BayesC, CSR, block-eigen, multivariate, wrappers, and namespace
remain unchanged from the Phase 15B baseline.

## 23. Runtime, memory, and I/O

Phase 11D/13E/14E are the within-model baselines. Phase 15B reruns use the common
reporting convention and compare only within each model. Completed-fit RSS is not
peak, page cache affects reads, and tiny timings are not performance claims. The
common layer adds no decoding, genotype copy, MCMC I/O, or hot-loop abstraction.

The final Windows/Rtools44 rerun found no unexplained material regression. BayesC
tiny medians were 0--0.01 s (one first-repetition outlier reached 0.54 s); the
2,000-marker/200-sample dense chain took 0.045--0.049 s per retained run, and
completed-fit RSS was approximately 124.5--135.8 MiB. BayesR tiny medians were
0--0.01 s, maxima excluding the first-load 0.29 s observation were at most 0.03 s,
and completed-fit RSS was 121.2--134.2 MiB. BayesRC tiny medians were 0--0.02 s,
maxima excluding the first-load 0.42 s observation were at most 0.02 s, and
completed-fit RSS was 122.7--134.0 MiB. Moderate/larger BayesR and BayesRC runs and
sampled peak RSS remain opt-in in their canonical scripts. BED files are decoded
once per fit and never reread during MCMC; repeated reads remain subject to the OS
page cache. Completed-fit RSS is not peak RSS, and these short timings are
regression signals rather than performance claims.

## 24. Tests

Phase 15A passed 64 focused expectations and Phase 15B passed 127. The enabled
canonical fresh-process selection passed 157 expectations. Canonical fixture
checks passed 3/3 raw and 3/3 formatted for each family (9/9 and 9/9 total), and
same-process, core-order, chain-count, worker-assignment, intervening-fit, and the
fixed-alpha BayesRC-to-BayesR reduction checks were exact. `compileAttributes()`
and a fresh native `load_all(..., compile=TRUE)` succeeded. The final full suite
reported 5,621 passes, zero failures, zero warnings, and eight intentional opt-in
skips.

## 25. Deviations and blockers

The proposed failure envelope was deliberately not added because production
exception/presentation semantics are not exact. This is a finalized conservative
decision, not a blocker. `git diff --check` passes. Normal Git-for-Windows
LF-to-CRLF checkout notices remain, but `git diff --stat` and
`git diff --ignore-space-at-eol --stat` are identical; no line-ending-only churn
or generated-wrapper change is present.

## 26. Canonicalization statement

The narrow common infrastructure is canonical. No separate Phase 15C is required.

## 27. Recommended next phase

Audit and dispose of the remaining experimental and sparse packed-BED BayesC
routes, combining the support decision and implementation when the evidence is
unambiguous, while leaving all canonical packed-BED family implementations unchanged.

## 28. Readiness marker

PHASE 15B COMPLETE — PACKED-BED COMMON INFRASTRUCTURE CONSOLIDATED AND CANONICALIZED
