# Unified BLR Framework Phase 11C3 report

## 1. Executive summary

The public scheduled multichain packed-BED BayesC migration is complete with
one typed per-chain numerical core, one binding-neutral native aggregation path,
one typed aggregate result, and one named R converter. Phase 11B deterministic
behavior and the public schema are unchanged.

## 2. Repository baseline

The clean baseline was branch `master` at Phase 11C2 commit `15965d7`. Initial
`git status --short` was empty and `git diff --check` passed. The Windows R 4.4
/ Rtools44 GNU C++17/OpenMP toolchain compiled the baseline. Phase 11C2 recorded
4,792 passing full-suite expectations, zero failures/warnings, and two opt-in
fresh-process skips. Its benchmark medians were 0.09 s for moderate dense 1x1
and 0.08 s for moderate aggressive 2x2, with completed-fit RSS about 122--139
MB; this RSS is not peak memory.

## 3. Phase 11C2 aggregation/conversion inventory

| Item | Phase 11C2 owner | Classification | Phase 11C3 owner |
|---|---|---|---|
| typed chain-result vector | adapter task loop | retain in task dispatch | adapter |
| trait/chain task mapping | adapter | retain in task dispatch | adapter |
| chain failure logging/exception translation | adapter | retain in task dispatch | adapter |
| marker means and final-state means | inline adapter | move to native aggregation | aggregator |
| marker SD/min/max | inline adapter | move to native aggregation | aggregator |
| variance, inclusion, VLE/VLD traces | inline adapter | move to native aggregation | aggregator |
| final variances and inclusion probability | inline adapter | move to native aggregation | aggregator |
| CPO, retained counts, timing means/maxima | inline adapter | move to native aggregation | aggregator |
| diagnostic `wy` and residual scores | inline adapter | move to native aggregation | aggregator |
| `ChainResultSTScheduled` alias | adapter shim | remove now | removed |
| marker-map alias | decoding/helper vocabulary | retain permanently | adapter helper |
| Rcpp matrix conversion lambdas | inline adapter | move to converter | converter |
| raw marker/trace/variance lists | inline adapter | move to converter | converter |
| schema, model/backend labels, MCMC metadata | inline adapter | binding-only metadata | converter metadata |
| actual R `NULL` placeholders and class | inline adapter | move to converter | converter |

No converter reached back into mutable execution state. The only numerical
data absent from the chain result were cross-chain summaries and optional
diagnostic genotype scores; both are now represented by the aggregate result.

## 4. Files changed

- `src/blr_bed_scheduled_bayesc_types.h`: adds typed aggregation context and
  aggregate execution result.
- `src/blr_bed_scheduled_bayesc_aggregate_impl.h`: adds the sole binding-neutral
  native aggregation callable.
- `src/st_cpg_omp_individual_scheduled_chains.cpp`: replaces inline aggregation
  and conversion with one aggregation call and one named converter.
- `tests/testthat/test-blr-framework-phase11c2.R`: replaces the obsolete
  inline-aggregation expectation with the permanent closed-boundary check.
- `tests/testthat/test-blr-framework-phase11c3.R`: adds permanent architecture,
  reference, reproducibility, schema, and protected-source checks.
- `tools/benchmarks/blr_phase11c3_bed_scheduled_bayesc.R`: establishes the final
  migrated benchmark workload.
- the implementation plan and capability matrix mark migration complete.
- this report records the closure evidence.

## 5. Final execution path

`stblr_bed(method="bayesc")` decodes and validates BED input, creates fit-owned
packed storage, enumerates trait-chain tasks, and dispatches them with unchanged
static OpenMP scheduling. Each task constructs a typed context and calls
`run_bed_scheduled_bayesc_chain()`. The typed chain results are aggregated once
by `aggregate_bed_scheduled_bayesc_results()`, converted once by
`stblr_bed_scheduled_bayesc_result_to_raw()`, and formatted through the existing
canonical raw-to-fit path.

## 6. Typed aggregate result

`BedScheduledBayesCExecutionResult` owns marker posterior means, PIP means,
SD/min/max summaries, averaged final effects/states, optional `wy` and residual
scores, marker/genetic/residual variance traces, inclusion traces, VLE/VLD,
final variance and inclusion values, CPO summaries, retained counts, timing
means/maxima, and marker/sample/trait/chain/task/trace dimensions. It contains
no R metadata, paths, genotype ownership, scheduler internals, or RNG state.

## 7. Native aggregation callable

`aggregate_bed_scheduled_bayesc_results()` validates dimensions, consumes the
chain-result vector in unchanged `chain * nt + trait` order, computes the same
arithmetic means, sample SD (`nchains - 1`), minima/maxima, traces, final-state
summaries, CPO, retained-count and timing summaries, then computes optional
diagnostic scores from immutable fit-owned genotype data. It runs once, does
not sample, and contains neither Rcpp nor RNG.

## 8. Task dispatch boundary

Trait-chain enumeration, logical seed resolution, context construction,
`schedule(static)` OpenMP dispatch, one chain-core call per task, result storage,
failure logging and exception translation remain in the adapter.

## 9. Named result converter

`stblr_bed_scheduled_bayesc_result_to_raw()` consumes only the typed aggregate
result and binding metadata (`nit`, `nburn`, `nthin`, `n_used`). It preserves
field order, types, dimensions, classes, actual R `NULL`, model/backend labels,
schema version, marker/trait/chain order, traces, CPO, diagnostics and selection
placeholders. It performs no aggregation, genotype access, I/O, or RNG work.

## 10. Temporary aliases and duplicate code removed

`ChainResultSTScheduled`, inline aggregate dimension aliases, duplicate
chain-index calculations, inline mean/SD/min/max and trace/CPO/timing formulas,
and inline raw-list construction were removed. `MarkerMapSTScheduledChains` is
retained as stable adapter/helper vocabulary because BED decoding and numerical
helper signatures still use it. No fallback or compatibility selector remains.

## 11. Final native adapter

The adapter is bounded to R decoding/public validation, BED/BIM/FAM validation,
fit-local genotype decoding, marker/trait preparation, task enumeration and
dispatch, logical seed mapping, typed context construction, one core call per
task, failure translation, one aggregation call and one converter call.

## 12. Per-chain numerical core

One `run_bed_scheduled_bayesc_chain()` implementation, one public-route MCMC
loop and one scheduler implementation remain in the guarded implementation
header. The core remains binding neutral and contains no disk access.

## 13. RNG ownership

Every logical trait-chain owns one `BedScheduledBayesCChainRng`, comprising one
`std::mt19937`, normal distribution and uniform distribution for exactly one
chain execution. There is no worker-owned, static, thread-local, shared or
fit-persistent distribution state. Aggregation/conversion consume no RNG and
seed formulas and draw sites are unchanged.

## 14. Genotype and I/O ownership

BED validation and blocked decoding remain in the adapter. Already-decoded
SNP-major packed storage is fit-owned and borrowed immutable by each chain and
by optional diagnostic aggregation. There is no per-chain genotype copy, core
disk access, extra decoding pass, or MCMC-time read.

## 15. Scheduler preservation

Iteration-zero/full-sweep timing, full-sweep frequency, adaptive null skipping,
candidate threshold/lifetime, due buckets, burn-in-only policy,
active/candidate/due order, marker order and skipped-marker no-RNG/genotype-work
behavior remain inside the unchanged per-chain core.

## 16. Corrected references

All three Phase 11B corrected raw references and all three formatted references
match exactly. Phase 11A BayesC remains historical and BayesR/BayesRC fixtures
remain unchanged.

## 17. Reproducibility

Exact `A; A`, `A; B; A`, normalized `1,2,2,1`, fresh/reused-process,
intervening BayesR, intervening BayesRC and different-chain-count behavior are
retained. The experimental single-chain and public multichain-one-chain routes
remain nonidentical because their established seed mappings differ.

## 18. Public API and schema

Public arguments/defaults/routing, native exported signature, generated
wrappers, `NAMESPACE`, `stblr_raw_v1`, formatted fields, types, dimensions,
classes, order and actual R `NULL` behavior are unchanged.

## 19. Unsupported behavior

There is still no public explicit-chain-seed argument, public scheduler-event
counter, scheduled packed-BED BayesR/BayesRC migration, or newly supported
configuration. The experimental route remains separate and nonidentical.

## 20. Protected backends

Packed-BED BayesR/BayesRC, experimental/sparse BayesC, canonical CSR,
block-eigen, multivariate sources, generated wrappers and `NAMESPACE` remain
unchanged from `15965d7`.

## 21. Performance, memory, and I/O

`Rscript tools/benchmarks/blr_phase11c3_bed_scheduled_bayesc.R` used the Phase
11C2 dense 1x1/2x1/2x2, aggressive, conservative and 2,000-marker/200-sample
workloads with warm-up and five timed repetitions. Moderate dense 1x1 had
times 0.05, 0.05, 0.04, 0.05 and 0.05 s (median 0.05, range 0.01); moderate
aggressive 2x2 had 0.04, 0.06, 0.07, 0.05 and 0.06 s (median 0.06, range 0.03).
These are consistent with Phase 11C2's 0.09/0.08-s medians and do not establish
an improvement. Completed-fit whole-process RSS was about 122--139 MB and is
not peak memory. The BED was 100,003 bytes, decoding/read strategy is unchanged,
and OS page-cache effects remain a limitation. No material runtime, memory or
I/O regression was observed.

## 22. Tests

The combined Phase 11A--11C3 focused set passes 215 expectations with two
expected opt-in skips; Phase 11C3 contributes 45. The enabled Phase 11B
fresh-process matrix passes 37/37 without skips. The final full suite passes
4,839 expectations with zero failures/warnings and two expected opt-in skips.
Native compilation and `Rcpp::compileAttributes()` succeeded with unchanged
generated wrappers and the existing non-fatal compiler warnings.

## 23. Deviations and blockers

No known numerical or architectural blocker. The marker-map alias is retained
because it is shared by adapter-side decoding helpers; removing it would be
unrelated cleanup. Seven Phase 9/10 historical hashes for the selected BED
source and the Phase 11A source hash were intentionally replaced/removed after
the migrated file changed; permanent typed-core/aggregation/converter checks
and hashes for genuinely protected sources replace that obsolete protection.

## 24. Recommended Phase 11D

> canonicalize and stabilize the public scheduled packed-BED BayesC route, retain the Phase 11B corrected fixtures permanently, remove remaining migration-era wording, establish the Phase 11C3 benchmark as the canonical baseline, and leave the experimental BayesC, BayesR, and BayesRC routes unchanged.

## 25. Readiness marker

PHASE 11C3 COMPLETE — PUBLIC SCHEDULED PACKED-BED BAYESC MIGRATED WITH DETERMINISTIC BEHAVIOR PRESERVED
