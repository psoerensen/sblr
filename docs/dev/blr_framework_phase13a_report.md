# Unified BLR Framework Phase 13A Report

## 1. Executive summary

Public packed-BED BayesR was audited without migrating or changing production
execution. Deterministic references, binding-neutral audit contracts, ownership
tests, reductions, benchmarks, and a future extraction seam are established.

## 2. Repository baseline

Branch `master`; starting/Phase 12A commit `8c41d18`; initial status and
`git diff --check` clean. R 4.4.1, Rtools44, GNU C++17/OpenMP. Baseline suite:
4,941 passes, no failures/warnings, three opt-in skips. Fast/extended CI exists.

## 3. Route and call graph

| classification | route |
|---|---|
| public production | `stblr_bed(method="bayesr")` |
| R helper | `.fit_stblr_bed_bayesr()` |
| native | `stblr_cpg_omp_bed_marker_scheduled_chains_bayesr()` |
| production source | `src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp` |
| per-chain helper | `run_one_bayesr_chain()` |
| shared helper | `st_bed_bayesr_common.h` decoding/math, also used by BayesRC |
| dispatch | adapter OpenMP `schedule(static)` over trait-chain jobs |
| aggregation/conversion | inline adapter blocks after dispatch |
| formatting | named `stblr_raw_v1` through canonical raw-to-fit formatter |

No second active BED BayesR sampler exists. Archived notes are historical; CSR
BayesR is separate. The future target is the public per-chain helper.

## 4. Statistical model

The mixture is a point-mass null plus K-1 zero-mean normals. Component k>0 has
variance `vb*c[k]`; `c[0]=0`. Stable log posterior weights combine chain-local
`pi[k]`, determinant and score terms. Inverse-CDF uniform sampling selects k;
a normal effect draw occurs only for k>0. Residuals update by effect difference.
`vb` uses `sum(b[j]^2/c[d[j]])`; `ve` uses residual SS. Scales are fixed.
When enabled, full sweeps update pi through Gamma(`alpha[k]+count[k]`,1) draws
normalized to a Dirichlet vector. Probabilities are chain/trait specific.

## 5. Component ordering and probability policy

| index/name | role | default scale | initial pi |
|---|---|---:|---:|
| 0/`component_0` | null | 0 | .95 |
| 1/`component_1` | non-null | .01 | .03 |
| 2/`component_2` | non-null | .1 | .015 |
| 3/`component_3` | non-null | 1 | .005 |

User order is preserved. `pi` is final, `pim` posterior mean, and
`comp_prob[[t]][j,k]` posterior assignment probability. `dm=1-P(k=0)`;
`component` is posterior mean index. Positive Dirichlet priors match K.

## 6. Genotype representation

The adapter validates magic/mode, SNP-major layout and BED/BIM/FAM ordering,
then blocked `fseek`/`fread` decoding creates fit-local packed marker-major
storage and marker maps. Missing values, centering and scaling remain in shared
lookup utilities. MCMC performs no file access or additional decoding.

## 7. Ownership and lifetime

| object | owner/lifetime | access |
|---|---|---|
| file handles | adapter decode | temporary, closed before MCMC |
| packed genotype/maps | fit | immutable shared, no chain copy |
| phenotype/prior/order | fit | borrowed immutable |
| task descriptors/results | adapter | one per trait-chain |
| effects/residuals/state | logical chain | mutable chain-local |
| scheduler state | logical chain | mutable chain-local |
| engine/distributions | logical chain | one chain execution |

## 8. Scheduler controls

`full_sweep_every>=0` (zero means every iteration), positive
`null_skip_base`, nonnegative `null_skip_max`, threshold in [0,1], nonnegative
candidate lifetime, and burn-in-only skipping flag. Skip length is a capped
piecewise function of non-null probability plus chain-RNG jitter. Sweep,
null-skip and candidate controls match canonical BED BayesC field semantics;
BayesR retains a wrapper because activity and initialization are component-based.

## 9. Mutable scheduler state

Chain-local due buckets, `scheduled_at`, `last_updated`, candidate/active flags
and lists, membership flags and `last_interesting` control traversal. Initial
nulls receive jittered due times; active assignments enter active/candidate
lists. Lists compact every 50 iterations. State is internal and not returned.

## 10. Marker traversal

Iteration zero and periodic full sweeps follow marker order. Adaptive traversal
is active, candidate, then due. `last_updated` prevents duplicates. Skipped
markers avoid genotype and categorical/effect RNG work. Visited markers draw one
uniform; normal draws are conditional on a non-null component.

## 11. RNG ownership

`std::mt19937`, uniform, cached normal and uniform-integer jitter distributions
are constructed per logical chain after its existing seed mapping. Variable-
shape Gamma distributions are local at pi updates; variance helpers use the
same chain engine. There is no static/thread-local, worker-owned, shared, or
cross-fit state and no `discrete_distribution` (selection uses inverse CDF).

## 12. Reproducibility

Exact after declared timing/core normalization: `A;A`, `A;B;A`, `1,2,2,1`,
one/multiple/one chains, intervening BayesC, and intervening BayesRC. Fresh
process 1x1 and 2x2 configurations also match exactly.

## 13. Frozen references

Three configurations: 1 chain/1 core/seed 71; 2 chains/1 core/seed 73; and 2
chains/2 cores/seed 73. Each stores raw, formatted fit, schema, components,
scheduler, toolchain and RNG-ownership metadata. Raw 3/3 and formatted 3/3 are
exact. Phase 11A references were not overwritten.

## 14. Component identities

Final/mean pi and marker component-probability rows sum to one; states lie in
`[0,K-1]`; `dm=1-P(null)`; effective non-null variance is `vb*c[k]`; fixed pi
remains supplied values within normalization rounding.

## 15. Reductions and nonreductions

Null-plus-one-non-null BayesR is not forced equal to BayesC: parameterization,
categorical/Bernoulli RNG order, probability policy and seeds differ. Two
full-sweep BayesR runs with skip bases 50 and 1 are nonidentical because initial
scheduler jitter consumes different RNG before traversal. Fixed-pi identity is
exact within floating normalization tolerance.

## 16. Aggregation

`ChainResultBayesR` owns marker summaries, component probabilities, effects/
state, traces/finals, CPO, retained count, timing/failure. Inline native code
maps `chain*traits+trait`, computes mean/sample-SD/min/max, traces, finals, CPO
and timings, then converts to R. Future aggregation should be one binding-neutral
callable consuming typed chain results; schema construction stays binding-side.

## 17. Raw and formatted schemas

Raw categories: marker bm/dm/wy/r/b/state and summaries; vb/vg/ve/vle/vld
traces/finals; final/mean pi; component probabilities/mean index/scales; CPO,
retained counts, timing/failure; stable placeholder categories and actual NULL.
Formatted output preserves bm/dm, component fields, pi aliases, traces,
diagnostics, marker/trait/component dimensions and stable fit fields.

## 18. Typed contracts

`blr_bed_bayesr_audit_types.h` defines component, scheduler, execution,
ownership, chain-result and aggregate-result vocabularies. It validates K/null/
scales, normalized pi, positive alpha, iteration/chain/core and scheduler
controls. It contains no Rcpp, SEXP, file handle, Python, sampler or RNG call.

## 19. Shared infrastructure decision

| BED BayesC infrastructure | decision |
|---|---|
| packed genotype view/ownership | reuse unchanged |
| sweep/null-skip/candidate controls | reuse unchanged in BayesR wrapper |
| task mapping/static scheduling | reuse conceptually, preserve behavior |
| seed/timing/failure vocabulary | reuse conceptually |
| BayesC chain/result types | do not reuse; semantic split |

## 20. Future extraction seam

Decoding ends after packed G, maps, marker order, y matrix and job allocation.
Execution begins at `run_one_bayesr_chain()` and ends at `ChainResultBayesR`.
Aggregation starts at `// Aggregate across chains`; conversion starts at
`// Build named raw schema v1`. Extract the per-chain helper first: storage is
fit-owned, mutable/RNG/scheduler state chain-local, and dispatch is separable.

## 21. Production behavior statement

BayesR production is byte-identical to `8c41d18`. Formulas, traversal, RNG,
seeds, decoding, OpenMP, aggregation, conversion, API and schemas are unchanged.

## 22. Performance, memory, and I/O baseline

The audit benchmark covers dense 1x1, adaptive 2x1/2x2, aggressive and
conservative skipping, warm-up and five repetitions. It reports distributions,
completed-fit RSS, BED size and page-cache caveats. Completed-fit RSS is not
peak; Phase 12A tooling provides sampled peak measurement. Tiny timings are not
speed claims. Dense 1x1 times were 0.28, 0, 0.01, 0.02 and 0.02 seconds
(median 0.02); adaptive 2x1 median was 0; the remaining tiny configurations had
medians 0--0.01 seconds. Completed-fit RSS ranged 121.75--131.42 MiB. A Windows
child-process peak attempt sampled 153,956,352 bytes but the child terminated
with status `-1073741819`; that sample is not promoted to a successful peak
baseline. The successful Phase 12A smoke validates the measurement route, while
Phase 13A records this platform limitation explicitly.

## 23. Protected backends

BayesR/BayesRC, canonical/experimental BED BayesC, canonical CSR, block-eigen,
multivariate, public routes, generated wrappers and `NAMESPACE` remain unchanged.

## 24. Tests

Contract/source, 3 raw/3 formatted references, reproducibility, component
identities, reductions/nonreductions, protected hashes and opt-in fresh-process
coverage are separate. Phase 13A focused: 57 passes and one opt-in skip; enabled
fresh-process: 61/61. Combined Phase 12A/13A focused: 105 passes and two opt-in
skips. Full suite: 4,998 passes, zero failures/warnings and four opt-in skips;
the Phase 13A fresh path was executed separately and passed.

## 25. Confirmed defects or risks

No correctness defect found. Migration risks: inline aggregation/conversion,
BayesR component state, scheduler-initialization jitter, and conditional
`Rcpp::Rcout` progress. Extraction must preserve them; console redesign is not
part of mechanical extraction.

## 26. Recommended Phase 13B

> mechanically extract the deterministic packed-BED BayesR per-chain numerical execution into one guarded implementation header while preserving Phase 13A references, component semantics, scheduler behavior, RNG ownership, genotype storage, aggregation, and public schemas.

## 27. Readiness marker

PHASE 13A COMPLETE — PACKED-BED BAYESR CONTRACT AND MIGRATION AUDIT READY
