# Unified BLR Framework Phase 14D Report

## 1. Executive summary

Packed-BED BayesRC migration is complete with one typed chain core, one typed
aggregate result, one binding-neutral aggregation path, one final marker-prior
recomputation path, and one named R converter. Phase 14A behavior is unchanged.

## 2. Repository baseline

The baseline was clean `master` at `68d085a` (`Activate typed packed-BED BayesRC
chain boundary`), which is the Phase 14C commit. R 4.4.1 and Rtools44 g++17 with
OpenMP were used. Phase 14C ended with 5,328 passing assertions, zero failures or
warnings, six opt-in skips, and exact 3/3 raw plus 3/3 formatted references. Its
tiny benchmark medians were 0.00--0.02 seconds and completed-fit RSS was
approximately 122.0--123.4 MiB.

## 3. Aggregation/conversion inventory

| Item | Phase 14D classification |
|---|---|
| typed task-result vector and job-to-trait/chain mapping | retain in task dispatch |
| marker means, PIP, component probabilities and final state | move to binding-neutral aggregation |
| trace and final variance summaries | move to binding-neutral aggregation |
| alpha means/finals and step-variance means/finals | move to binding-neutral aggregation |
| final marker priors and component counts | move to final-prior/aggregation path |
| CPO, retained counts and failure collection | move to binding-neutral aggregation |
| optional `wy` and residual marker score | retain as post-aggregation genotype diagnostic |
| component labels, R orientations, lists, classes and actual `NULL` | move to converter |
| model/backend/schema and MCMC metadata | binding-only metadata |
| inline aggregation and raw-list fragments | remove now |

The production route has no timing/progress fields beyond its existing diagnostic
vocabulary, so no new public fields were introduced.

## 4. Files changed

`blr_bed_bayesrc_types.h` adds aggregation and aggregate-result contracts;
`blr_bed_bayesrc_aggregate_impl.h` adds the singular native aggregation;
the production adapter now dispatches, aggregates once, fills optional genotype
diagnostics, and converts once. Phase 14 tests, benchmark, plan, matrix, and this
report were updated. Earlier tests changed only where inline-boundary assumptions
or the intentional adapter hash became obsolete.

## 5. Final execution path

`stblr_bed(method="bayesrc")` performs R annotation preprocessing, calls the
native BED adapter, decodes fit-local genotype storage, dispatches typed chain
contexts to `run_bed_bayesrc_chain()`, passes typed chain results to
`aggregate_bed_bayesrc_results()`, optionally fills genotype-derived diagnostics,
and calls `stblr_bed_bayesrc_result_to_raw()` before the unchanged formatter.

## 6. Typed aggregate result

`BedBayesRCExecutionResult` owns marker means/PIP/component probabilities/final
states, traces and final variances, alpha and step-variance summaries, final and
posterior marker priors, component counts, CPO/retained counts, failures,
dimensions, optional diagnostics, and retained chain results when requested. It
owns no RNG, scheduler, genotype, annotation, file, R, or schema object.

## 7. Native aggregation

`aggregate_bed_bayesrc_results()` consumes trait-major slots indexed by
`chain * ntraits + trait`. It preserves the original trait-then-chain addition
order, averages every former averaged field by `1 / nchains`, derives component
counts from averaged component probabilities, and collects failures. No R API,
RNG, decoding, latent update, scheduler, progress output, or I/O is present.

## 8. Final marker-prior recomputation

For each chain final alpha (`annotations x (K-1)`), aggregation invokes the sole
`st_bayesrc_compute_snp_pi()` helper on the prepared marker-by-annotation matrix.
That helper uses `StandardNormalProbability::cdf()`, unchanged annotation/stick
order, caller-side clipping, flooring, residual final stick and row
renormalization. The resulting marker-by-component matrix and its column means
are aggregated in the original chain order.

## 9. Task dispatch

OpenMP `schedule(static)`, task count, `job % ntraits`, `job / ntraits`, logical
seed resolution, typed context construction, exception boundary, and one core
call per task remain adapter-owned and unchanged.

## 10. Progress handling

This production route has no chain progress-event stream. Any console/error
translation remains outside numerical aggregation; aggregation emits nothing and
has no numerical side effects.

## 11. Optional genotype-derived diagnostics

`wy` and the residual marker score remain after typed aggregation because they
require packed-genotype access and binding flags rather than cross-chain
statistics. They fill optional fields on the aggregate result before conversion;
the aggregate type does not own genotype storage.

## 12. Named result converter

`stblr_bed_bayesrc_result_to_raw()` is the sole converter. It consumes the typed
aggregate plus binding metadata and preserves field order, storage types,
orientations, labels, classes, placeholders, actual R `NULL`, and
`stblr_raw_v1` compatibility. It performs no aggregation or prior recomputation.

## 13. Temporary aliases and duplicate code

Inline marker/trace/alpha aggregation, duplicate final-prior evaluation,
component-count calculation, diagnostic scaling, retained-chain conversion, and
raw-list construction were removed from the exported function. The `jobs` name
is retained as adapter-owned typed task storage; `BedBayesRCBindingMetadata` is a
stable binding-only vocabulary, not an execution alias.

## 14. Final native adapter

The adapter is limited to binding validation, annotation/BED preparation,
component/prior/task construction, dispatch, seed resolution, one aggregate
call, optional genotype diagnostics, one converter call, and exception
translation. It contains no MCMC, latent/alpha update, RNG draw, cross-chain
formula, final-prior formula, or raw-list construction outside the converter.

## 15. Per-chain numerical core

One `run_bed_bayesrc_chain()` implementation, one MCMC loop, one exact full
sweep, one component sampler, one latent truncated-normal path and one sequential
alpha update remain. The numerical headers remain binding neutral.

## 16. Component and annotation semantics

Component zero is the point-mass null with `gamma[0]=0`; active variance is
`vb*gamma[k]`; `K-1` ordered sticks leave component `K-1` as the residual stick.
Assignments remain zero-based, alpha remains `P x (K-1)`, intercept index remains
zero, marker priors remain marker-by-component, and `dm=1-P(null)`.

## 17. Normal-probability boundary

The stateless C++ `StandardNormalProbability` interface remains backed by the
same `R::pnorm`/`R::qnorm` conventions. Core and aggregation code call no R
namespace probability primitive directly.

## 18. RNG ownership

Each logical chain owns one `std::mt19937`, uniform and standard-normal
distributions, plus parameterized local update distributions. Seed mapping and
draw sites are unchanged. Workers, fits, aggregation and conversion own no RNG.

## 19. Genotype, annotation, and I/O ownership

Genotype and prepared annotation matrices are decoded/prepared once per fit and
borrowed immutably by chain contexts. No full per-chain copy or numerical disk
access exists. Alignment and factor/intercept handling remain adapter-side.

## 20. Exact references

Phase 14A raw references: 3/3 exact. Phase 14A formatted references: 3/3 exact.
No fixture was regenerated.

## 21. Reproducibility

Repeated A, A/B/A, normalized 1/2/2/1 core order, one/multiple/one chain order,
fresh versus reused process, intervening BayesR/BayesC, and worker assignment are
exact after only the declared execution-metadata normalization.

## 22. Annotation and probability identities

Stick and component probabilities are finite/bounded/nonnegative; rows sum to
one; the residual-stick and `dm` identities hold; assignments are valid; alpha
dimensions/intercept order, factor expansion, fixed-alpha/fixed-prior and
zero-alpha ordered-stick policies remain exact. Constant non-intercept and
missing annotations remain rejected.

## 23. Reduction and nonreductions

The intercept-only fixed-alpha configuration remains exactly equal to matched
fixed-pi packed-BED BayesR. General BayesR and CSR SBayesRC remain documented
policy/execution nonreductions; their first differences are unchanged.

## 24. Public API and schema

Public arguments, defaults, dispatch, native export signature, raw field order,
types, dimensions, actual `NULL`, formatted fields, wrappers, and `NAMESPACE`
are unchanged.

## 25. Protected backends

Canonical packed-BED BayesR, all packed-BED BayesC routes, CSR SBayesRC/BayesR
and other CSR backends, block-eigen, multivariate sources, generated wrappers,
and `NAMESPACE` remain byte-identical to the starting commit.

## 26. Performance, memory, and I/O

The Phase 14D benchmark repeats the Phase 14C tiny workloads after warm-up with
five repetitions. Medians were 0.00, 0.02, 0.02 and 0.03 seconds; completed-fit
RSS ranged from 121.66 to 122.42 MiB. One 0.75-second Windows timer outlier made
the intercept-only mean 0.156 seconds; the other means were 0.020--0.024 seconds.
Completed-fit RSS is not sampled peak RSS; Windows timer resolution and page
cache dominate tiny runs. No extra decoding or numerical I/O was introduced and
no speedup is claimed.

## 27. Tests

Phase 14A--C references/architecture tests, new Phase 14D architecture and
aggregation tests, fresh-process matrix, identities, reduction, and protected
backends pass. The enabled Phase 14A--D fresh-process run has no skip or failure.
The full suite reports 5,371 passes, zero failures, zero warnings and seven
declared opt-in skips.

## 28. Deviations and blockers

No numerical or schema deviation is accepted. The aggregate result retains chain
results only when `keep_chains=TRUE`; this preserves the public chain schema while
avoiding the copy in the default route. No blocker remains.

## 29. Recommended Phase 14E

Canonicalize and stabilize public packed-BED BayesRC, retain the Phase 14A
fixtures permanently, remove remaining migration-era wording and aliases,
establish the Phase 14D benchmark as the canonical runtime/completed-fit-RSS/I/O
baseline, and leave all canonical BayesC and BayesR implementations unchanged.

## 30. Readiness marker

PHASE 14D COMPLETE — PACKED-BED BAYESRC MIGRATED WITH DETERMINISTIC BEHAVIOR PRESERVED
