# Unified BLR Framework Phase 11D report

## 1. Executive summary

The public scheduled multichain packed-BED BayesC route is canonicalized and
stabilized. Its Phase 11C3 native architecture required no source changes:
canonicalization consists of permanent structural/reference protections,
canonical benchmark and documentation baselines, and explicit classification
of remaining historical material.

## 2. Repository baseline

The clean baseline was branch `master` at Phase 11C3 commit `2af74c7`. Initial
`git status --short` was empty and `git diff --check` passed. R 4.4.1 with the
Rtools44 GNU C++17/OpenMP toolchain compiled and loaded the package. The
baseline full suite passed 4,839 expectations with zero failures/warnings and
two opt-in fresh-process skips. Re-running the Phase 11C3 benchmark after the
debug rebuild produced moderate medians of 0.19 s for dense 1x1 and 0.19 s for
aggressive 2x2; completed-fit RSS was approximately 123--139 MB and the BED was
100,003 bytes. These order/build-sensitive timings are retained without a
speed claim.

## 3. Phase 11C3 structure inventory

| Item | State at baseline | Phase 11D classification |
|---|---|---|
| `BedPackedGenotypeView` | borrowed fit-owned storage | retain permanently |
| typed per-chain context/result | one each | retain permanently |
| callable chain core/MCMC/scheduler | one each | retain permanently |
| OpenMP trait-chain dispatch | adapter, `schedule(static)` | retain permanently |
| typed aggregate result/aggregator | one each | retain permanently |
| named R converter | one | retain permanently |
| `BedScheduledBayesCChainRng` | logical-chain owned | retain permanently |
| Phase 11B corrected fixtures | exact | canonical permanent fixtures |
| Phase 11A defective BayesC fixture | pre-correction | historical audit artifact |
| `MarkerMapSTScheduledChains` alias | decoding-helper vocabulary | retain; stable helper name |
| commented pre-extraction source snapshots | inactive historical text | defer to cross-cutting hardening |
| Phase 11C reports/tests | phase evidence | retain as historical phase gates |
| “migrated/ready” plan and matrix wording | stale status | rename for canonical use |
| Phase 11C3 benchmark wording | migration baseline | consolidate into canonical 11D baseline |

No unused context, chain-result or aggregate-result field, duplicate active
aggregation/conversion fragment, fallback, selector or temporary native alias
was found. Deleting the large inactive commented snapshots would create an
unrelated native-source diff and is explicitly deferred.

## 4. Files changed

- `tests/testthat/test-blr-framework-phase11d.R` adds permanent architecture,
  reference, reproducibility, header-safety and protected-source checks.
- `tools/benchmarks/blr_phase11d_bed_scheduled_bayesc.R` establishes the
  canonical runtime/completed-fit-RSS/I/O workload.
- the implementation plan and capability matrix mark only the public
  multichain route canonical.
- this report records the canonical ownership and validation baseline.

No native production source, public wrapper, generated wrapper or namespace
file was modified.

## 5. Canonical execution path

`stblr_bed(method="bayesc")` performs existing R validation and BED alignment,
decodes packed genotype once per fit, enumerates trait-chain tasks, dispatches
them with static OpenMP scheduling, constructs one typed context per logical
chain, calls `run_bed_scheduled_bayesc_chain()`, aggregates typed chain results
once, converts the typed aggregate once to `stblr_raw_v1`, and uses the existing
raw-to-fit formatter.

## 6. Cleanup and naming

Stable architecture names are retained unchanged. Canonical plan/matrix and
benchmark wording replaces migration-in-progress wording. No dead native alias
or duplicate active code existed to remove. Phase 11A worker-cache material is
historical; Phase 11C reports/tests remain migration evidence; inactive
commented source snapshots are deferred to the recommended hardening phase.

## 7. Packed genotype view

`BedPackedGenotypeView<PackedGenotype>` borrows the fit-owned storage object and
const packed-byte pointer plus packed size, marker/sample counts, bytes per
marker and stride. It owns no bytes, path or file handle and has no binding or
Python type. Its lifetime is bounded by the enclosing fit and outlives all task
calls.

## 8. Typed per-chain context

Borrowed immutable fields are packed genotype, marker maps/order, phenotype,
initial effects/state and prior matrices/vectors. Copied fields are scheduler,
MCMC/model controls, update flags, resolved chain seed, trait index and chain
index. It contains no binding metadata, file state, worker identity or mutable
cross-chain state.

## 9. Per-chain numerical core and implementation header

The guarded implementation-detail header is included only by the intended
binding translation unit. It contains one callable core, one active public MCMC
loop, one scheduler transition implementation and one chain RNG construction.
It contains no Rcpp/Python type, R conversion, file opening or disk access and
does not redefine compilation configuration.

## 10. Typed per-chain result

The result owns marker posterior/PIP means, final effects and indicators,
variance/inclusion/VLE/VLD traces and finals, CPO values, retained count,
timing and failure information used by aggregation. It exposes no genotype,
scheduler or RNG ownership and no binding metadata.

## 11. Task dispatch

The adapter retains `chain * nt + trait` task mapping, unchanged logical seed
formula, static OpenMP assignment, one context/core call per task, ordered
result storage, failure logging and exception translation.

## 12. Typed aggregate result

`BedScheduledBayesCExecutionResult` owns all converter-consumed numerical
marker means, SD/min/max, final states, diagnostic scores, traces, variance and
inclusion finals, VLE/VLD, CPO, retained counts, timing summaries and dimensions.

## 13. Native aggregation

Exactly one `aggregate_bed_scheduled_bayesc_results()` consumes typed chain
results in established task order. Arithmetic means, sample SD using
`nchains-1`, min/max, traces, final-state rules, diagnostic scores, CPO,
retained counts and timing formulas are unchanged. It contains no Rcpp, RNG,
scheduler transition, decoding or disk access.

## 14. Native adapter

The adapter is bounded to R/public validation, BED validation/decoding,
marker/trait/task preparation, dispatch/seed resolution, context/core calls,
failure translation, one aggregation call and one converter call. It contains
no active MCMC, marker update, RNG draw, posterior or aggregation formula.

## 15. Result converter

Exactly one `stblr_bed_scheduled_bayesc_result_to_raw()` remains binding-only.
It preserves field order/names, storage types, dimensions, classes, marker and
trait order, actual R `NULL`, schema/backend/model identifiers, posterior and
variance outputs, inclusion traces, VLE/VLD, CPO, timing, failures, chain
summaries, placeholders and input metadata.

## 16. RNG ownership

Each logical trait-chain owns one `BedScheduledBayesCChainRng`, one
`std::mt19937`, one normal distribution and one uniform distribution for one
chain execution. Worker ownership, sharing, fit persistence and active
static/thread-local distributions are absent. Local variable-parameter draws,
seed formula and logical draw sites are unchanged. Phase 11A documents the
historical worker-cache defect.

## 17. Scheduler semantics

The first/full-sweep timing and `full_sweep_every` are unchanged. Non-full
traversal remains active, candidate, then due markers. Null-skip base/max,
probability-adaptive growth/reset, burn-in-only option, candidate threshold,
lifetime/decrement/expiry and marker ordering are unchanged. Skipped markers
consume no marker-update RNG and retain established genotype-work avoidance.
Neighbor wake-up remains unsupported.

## 18. Genotype and I/O ownership

The adapter retains magic/mode and BED/BIM/FAM consistency checks, blocked
`fseek`/`fread`, SNP-major decoding, order, missing-value handling and scaling.
Packed storage is immutable and shared across chains. No per-chain full copy,
additional decoding pass, numerical-core file state or MCMC-time disk read is
present. Page-cache effects apply to benchmark interpretation.

## 19. Permanent regression fixtures

Permanent corrected configurations are single-chain 1x1, multichain 2x1 and
multichain 2x2 with deterministic mapping and adaptive controls, each carrying
raw/formatted schema, actual-NULL and CPO coverage. Fixture regeneration remains
manual maintenance only. Phase 11A BayesC stays historical; its trajectory is
not an ordinary expected value.

## 20. Exact reference results

Corrected raw references: 3/3 exact. Corrected formatted references: 3/3 exact.

## 21. Reproducibility

`A; A`, `A; B; A`, normalized `1,2,2,1`, fresh/reused process, intervening
BayesC/BayesR/BayesRC, different-chain-count and worker-assignment checks are
exact after normalizing only declared execution metadata. Experimental
single-chain versus public one-chain remains nonidentical because seed mappings
differ.

## 22. Public API and schema

Arguments, defaults, routing, native signatures, generated wrappers,
`NAMESPACE`, `stblr_raw_v1`, formatted fit structure and actual R `NULL`
behavior are unchanged.

## 23. Unsupported behavior

No public explicit chain seeds, scheduler counters or neighbor wake-up were
added. Experimental/sparse BayesC, BayesR and BayesRC remain separate,
noncanonical and otherwise unchanged. Unsupported models/configurations remain
unsupported.

## 24. Protected backends

Git/hash audits against `2af74c7` confirm no change to experimental or sparse
BayesC, packed-BED BayesR/BayesRC, canonical ordinary/fixed/group/annotation/
scheduled CSR, block-eigen, multivariate implementations, public routes,
generated wrappers or `NAMESPACE`.

## 25. Performance, memory, and I/O baseline

`Rscript tools/benchmarks/blr_phase11d_bed_scheduled_bayesc.R` uses warm-up and
five repetitions for dense 1x1/2x1/2x2, aggressive and conservative skipping,
plus 2,000-marker/200-sample moderate cases. Individual times and summary
statistics are printed with R/package/toolchain metadata. On R 4.4.1/sblr
0.1.0, the moderate dense 1x1 times were 0.17, 0.18, 0.17, 0.18 and 0.21 s
(mean 0.182, median 0.18, range 0.04 s); moderate aggressive 2x2 times were
0.19, 0.21, 0.21, 0.19 and 0.22 s (mean 0.204, median 0.21, range 0.03 s).
Tiny-workload results were dominated by warm-up and timer resolution. The
100,003-byte BED fixture used blocked reading. Whole-process completed-fit RSS
ranged from 122.6 to 136.0 MiB across reported configurations; this is not peak
RSS. OS page-cache and run-order effects remain explicit limitations. Compared
with the Phase 11C3 rerun on the same debug build (moderate medians 0.19 s), no
material runtime, completed-fit-RSS or I/O regression is evident, and no speed
improvement is claimed.

## 26. Tests

Phase 11A--11D focused tests passed 268 expectations with two opt-in
fresh-process skips; the separately enabled fresh-process Phase 11A/11B matrix
passed 111/111. Phase 11D itself passed 52/52. The final full suite passed 4,891
expectations with zero failures, errors or warnings and two opt-in skips (4,893
total expectations). Phase 11A, 11B, 11C1, 11C2 and 11C3 protections, BED
focused tests and protected-backend tests all passed within those runs.

## 27. Deviations and blockers

No canonicalization blocker. Peak/interval RSS sampling is unavailable, so the
baseline accurately reports completed-fit RSS. Historical commented source is
classified and deferred rather than changed in this bounded phase.

## 28. Recommended next phase

> perform a cross-cutting architecture and documentation hardening pass covering shared validation consistency, structural-test robustness, binding/core diagnostic separation, stale plan metadata, historical commented source, CI coverage, and peak-memory measurement before beginning the next BED-model migration.

## 29. Readiness marker

PHASE 11D COMPLETE — PUBLIC SCHEDULED PACKED-BED BAYESC CANONICALIZED AND STABILIZED
