# Unified BLR Framework Phase 14E Report

## 1. Executive summary

Public packed-BED BayesRC is canonicalized and stabilized. Its typed chain core,
aggregation, final-prior calculation, converter, fixtures, and resource baseline
are permanent without changing numerical behavior or public schemas.

## 2. Repository baseline

The clean baseline was `master` at Phase 14D commit `4083763` (`Complete
packed-BED BayesRC migration`). R 4.4.1, Rtools44 g++17, and OpenMP were used.
The baseline full suite had 5,371 passes, zero failures/warnings and seven opt-in
skips. Phase 14D medians were 0.00--0.03 seconds and completed-fit RSS was
121.66--122.42 MiB; one Windows timing outlier was 0.75 seconds.

## 3. Phase 14D structure inventory

The baseline contained one component, annotation, coefficient-prior, packed
genotype, chain-context, chain-result and aggregate-result vocabulary; one chain
core; adapter-owned static task dispatch; one Rmath-backed probability boundary;
one aggregate; one final-prior helper call site; optional post-aggregate genotype
diagnostics; one converter; and three raw/formatted fixture pairs.

## 4. Migration-artifact classification

| Classification | Items |
|---|---|
| retain permanently | all stable production type names, core, aggregation, probability helper, converter, task mapping, fixtures |
| remove now | unused aggregate `has_wy` and `has_residual_marker_score` flags and assignments |
| rename for canonical use | benchmark/report status wording only |
| consolidate | Phase 14E permanent structural/reference coverage; no numerical helper consolidation needed |
| historical audit artifact | Phase 14A audit types/reports, Phase 14B extraction report, Phase 14C/D migration reports |
| defer with justification | broader packed-BED family consolidation; outside this phase and risks canonical cores |

No fallback, selector, duplicate aggregation, duplicate converter, or competing
final-prior implementation was found.

## 5. Files changed

Native edits remove only two proven-unused boolean fields and assignments.
Phase 14E adds permanent tests, a canonical benchmark and this report; the plan,
matrix and intentionally obsolete structural hashes/status assertions are updated.

## 6. Diff and line-ending hygiene

The starting tree had no diff. `git diff --check`, `--numstat`, and
`--ignore-space-at-eol --stat` were inspected. Intentional content changes are
limited to the files reported here; no line-ending-only file is retained.
Windows reports prospective LF-to-CRLF checkout notices for intentionally edited
text files, but no whitespace error or whole-file EOL churn remains.

## 7. Canonical execution path

`stblr_bed(method="bayesrc")` validates and prepares annotations, validates and
decodes BED storage, statically dispatches typed chain contexts, runs
`run_bed_bayesrc_chain()`, aggregates typed results, fills optional genotype
diagnostics, converts through `stblr_bed_bayesrc_result_to_raw()`, and invokes the
unchanged formatter.

## 8. Component specification

Component zero is the point-mass null with exact zero scale; active variance is
`vb*gamma[k]`; `K-1` ordered sticks leave component `K-1` as the residual stick.
Ordering and finite-positive active-scale validation remain exact.

## 9. Annotation specification

The numerical view is immutable marker-by-annotation storage with explicit
marker/annotation counts and intercept index zero. Alignment, factor expansion,
names, missingness checks and reordering are excluded and remain adapter-side.

## 10. Coefficient-prior specification

Alpha remains `P x (K-1)`. Initial step variances, flat-intercept policy,
inverse-chi-square controls, update flag/frequency, ordering, positivity and
finite-value validation remain unchanged.

## 11. Immutable views and ownership

Genotype and prepared annotation storage are fit-owned and borrowed immutably.
They outlive parallel chain calls and are never copied fully per chain. Chain
effects, residuals, assignments, probabilities, alpha, latent/workspace state,
variances and stochastic objects are private.

## 12. Per-chain context

The context borrows genotype, annotation, phenotype, marker maps/order, component
and coefficient-prior inputs and initial state; it copies scalar MCMC/update
controls, the resolved seed, and logical trait/chain indices. It excludes R,
schema, paths, handles, workers and cross-chain state.

## 13. Per-chain numerical core

One core contains one MCMC loop, exact full sweeps, component/probit/latent/alpha
and variance updates, posterior/CPO accumulation and finalization. It contains no
binding, alignment, decoding, task, aggregation, schema or I/O work.

## 14. Per-chain result

The result owns marker summaries, final effect/state/residual, component
probabilities, traces, alpha and step-variance means/finals, prior means, final
variances, CPO, retained count and failure text. It excludes RNG, mutable latent
state, genotype/annotation ownership and binding metadata.

## 15. Normal-probability boundary

`StandardNormalProbability::cdf()` and `quantile()` are stateless C++ signatures
backed by the original lower-tail, non-log `R::pnorm`/`R::qnorm` calls. Clipping
remains at established callers. No approximation, lookup table or alternate path
exists; core and aggregation contain no direct R probability call.

## 16. Task dispatch

Tasks remain `chain*ntraits+trait`, OpenMP scheduling remains static, and seeds
remain `seed+1000003*(trait+1)+9176*(chain+1)`. Worker assignment does not own or
alter stochastic state.

## 17. Typed aggregate result

The aggregate owns all numerical marker, trace, variance, component, annotation,
prior, diagnostic, failure and requested retained-chain data required by the
converter. Optional genotype-derived matrices are filled after aggregation.

## 18. Native aggregation

One binding-neutral callable preserves trait-then-chain order, original averaging,
marker/component/alpha/trace/CPO summaries, component counts and failures. The
current raw route defines no additional sample-SD/min/max or timing fields, so no
new public summaries were invented.

## 19. Final marker-prior recomputation

One aggregation call site passes prepared annotations and each typed final alpha
to the established probit-stick helper. Annotation/stick/component order,
Rmath CDF, clipping, flooring, residual stick, normalization and marker-by-
component orientation remain exact. No adapter recomputation remains.

## 20. Optional genotype diagnostics

After aggregation the adapter optionally computes `wy` and residual marker scores
using prepared genotype storage. These are not cross-chain statistics, consume no
RNG, do not mutate canonical summaries, and preserve actual-`NULL` behavior.

## 21. Result converter

The sole named converter consumes the aggregate and binding metadata and preserves
all raw fields, ordering, storage types, dimensions, labels, classes, nullable fields,
actual R `NULL`, schema/version and backend identifiers. It does no aggregation.

## 22. Final adapter

The adapter is bounded to R validation, annotation and BED preparation, task/prior
construction, static dispatch, seed/context creation, failure handling, one
aggregate call, optional diagnostics, one converter call and exception translation.

## 23. RNG ownership

Each logical chain owns one `std::mt19937`, one uniform and standard-normal
distribution, and local parameterized update distributions. No static,
thread-local, worker, fit-persistent, aggregate or converter RNG exists.

## 24. Full-sweep and probability semantics

Every iteration visits every marker once in production order. No scheduler, skip,
candidate, active, due or jitter state exists. Probit sticks preserve predictor,
CDF, complement/product, residual-stick, flooring and row-normalization order.

## 25. Latent and alpha updates

Indicator and truncated-normal branch order, CDF bounds, clipping, uniform and
inverse-CDF draw timing, marker/stick order, sequential alpha residualization,
flat intercept, conditional moments, draw sites and traces remain unchanged.

## 26. Genotype, annotation, and I/O

BED and annotation preparation occur once per fit. Immutable storage is shared;
there is no numerical reread, file/path handle, alignment pass, or full per-chain
copy in core or aggregation.

## 27. Permanent fixtures

Canonical fixtures cover one chain/one core fixed intercept alpha, two chains/one
core updated three-column annotations, and two chains/two cores with the same
updates. They cover four components, probabilities, alpha means/finals, CPO,
diagnostics, actual R `NULL`, and schema/version. Generation remains manual only.

## 28. Exact references

Raw references: 3/3 exact. Formatted references: 3/3 exact. No expected value was
regenerated or broadly normalized.

## 29. Reproducibility

A/A, A/B/A, normalized 1/2/2/1 cores, one/multiple/one chains, fresh/reused
processes, intervening BayesR/BayesC and worker assignment remain exact after only
declared execution-metadata normalization.

## 30. Identities

Stick bounds, component finiteness/nonnegativity, row sums, residual stick,
`dm=1-P(null)`, state range, alpha dimensions, intercept order, alignment,
factor expansion, fixed/zero-alpha policies and validation rejections pass.

## 31. Reduction and nonreductions

Intercept-only fixed-alpha BayesRC remains exactly reducible to matched fixed-pi
packed-BED BayesR. General global-probability BayesR and CSR SBayesRC remain
policy/execution nonreductions with unchanged first differences.

## 32. Public API and schema

Public arguments/defaults/routing, native signature, raw and formatted schemas,
actual R `NULL`, generated wrappers and `NAMESPACE` remain unchanged.

## 33. Unsupported behavior

Missing annotations, adaptive scheduling, public explicit chain seeds, public
scheduler counters, new annotation policies and alternate component order remain
unsupported exactly as before.

## 34. Protected backends

Canonical packed-BED BayesR/BayesC, experimental/sparse BayesC, CSR SBayesRC,
CSR BayesR and other CSR backends, block-eigen, multivariate sources, wrappers and
`NAMESPACE` are byte-identical to the Phase 14E starting commit.

## 35. Performance, memory, and I/O baseline

The Phase 14E benchmark repeats Phase 14D tiny fixed/updated-alpha 1x1, 2x1 and
2x2 workloads after warm-up for five repetitions, recording individual and
summary timings, completed-fit RSS, annotation size, versions and limitations.
Observed five-run timing summaries (seconds) were: intercept fixed 1x1
`0.98, 0.01, 0.03, 0.02, 0.02` (median 0.02; range 0.97), updated-alpha 1x1
`0.02, 0.01, 0.02, 0.03, 0.02` (median 0.02; range 0.02), updated-alpha 2x1
`0.02, 0.01, 0.02, 0.02, 0.03` (median 0.02; range 0.02), and updated-alpha
2x2 `0.02, 0.01, 0.03, 0.01, 0.03` (median 0.02; range 0.02). Completed-fit
RSS was 128.72, 121.79, 121.82, and 121.77 MB respectively. Annotation object
sizes were 232 bytes for intercept-only and 264 bytes for three annotations.
Moderate/larger workloads and sampled peak RSS remain explicit opt-ins.
Completed-fit RSS is not peak RSS; tiny timings are regression signals rather
than performance claims. Repeated BED reads are page-cache affected; the tiny
fixture does not expose a standalone BED byte count. No new decode pass or
MCMC-time I/O exists, and no unexplained material difference from Phase 14D was
observed.

## 36. Tests

The Phase 14E focused file passed 57 assertions with one opt-in fresh-process
skip. Phase 14A--E focused checks passed 259 assertions with four opt-in skips;
with fresh processes enabled they passed all 267 assertions with no skips.
The final full suite passed 5,428 assertions with zero failures, zero warnings,
and eight documented opt-in skips. Compilation and load succeeded. The three
raw and three formatted Phase 14A references were exact, as were the permanent
same-process, fresh-process, core-order, chain-count, worker-independence,
identity, and fixed-alpha reduction checks.

## 37. Deviations and blockers

The only native cleanup removes two unused availability flags; it changes no
operation, RNG draw, aggregation value, converter output or schema. No blocker
remains.

## 38. Recommended next phase

Begin a cross-backend packed-BED family consolidation audit covering shared
genotype views, task mapping, seeds, failure/timing vocabulary, converter
conventions, and permanent architecture contracts while leaving the canonical
BayesC, BayesR, and BayesRC numerical cores unchanged.

## 39. Readiness marker

PHASE 14E COMPLETE — PACKED-BED BAYESRC CANONICALIZED AND STABILIZED
