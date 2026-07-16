# Unified BLR Framework Phase 9G Report

## 1. Executive summary

Learned-annotation CSR BayesC is canonicalized and stabilized with one typed numerical core, one marker loop, one binding converter, one wrapper aggregation path, permanent exact fixtures, and a canonical runtime/completed-fit-RSS baseline.

## 2. Repository baseline

- Branch/start/Phase 9F3 commit: `master`, `1b4a0e4` (`Migrate learned-annotation CSR BayesC to typed core`).
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1 UCRT, Rtools44 GCC 13.2.0.
- Baseline fresh-build suite: 4,327 passes, 0 failures, warnings, or skips.
- Phase 9F3 pre-cleanup benchmark, 2,000-marker both-learning means: 0.562 s (1 chain/1 core), 1.122 s (2 chains/1 core), 1.102 s (retained 2 chains/2 cores); completed-fit RSS about 128--129 MiB.

## 3. Phase 9F3 structure inventory

| Item | Classification | Phase 9G disposition |
|---|---|---|
| stable context/result/core/converter names | retain permanently | unchanged |
| guarded implementation header and TU-only inclusion macro | retain permanently | required build-safety mechanism |
| original-name aliases at top of numerical core | defer with explicit justification | they preserve the validated numerical body and expose every typed dependency clearly |
| context fields and result `raw` payload | retain permanently | every field is used by validation/execution; the result is used by aggregation/conversion |
| marker/trace/annotation conversion helpers | retain permanently | binding-only, nonduplicated shape conversion |
| wrapper aggregation aliases and accumulators | retain permanently | required for exact chain summaries, coefficient payloads, LD-swap counts, and `keep_chains` |
| duplicate converter/aggregation paths | remove now | none found |
| `cpg_annot_raw_v1()` compatibility wording/path | remove now | none remains |
| old/new selector, fallback, lexical-path comments | remove now | none remains |
| Phase-numbered historical tests/reports | retain permanently | historical audit trail and exact fixtures |
| post-migration benchmark wording | rename for canonical use | new Phase 9G canonical entry freezes the Phase 9F3 workload |
| stale documentation calling learned annotation unmigrated | rename for canonical use | replaced with canonical status |

No native cleanup was justified: removing aliases or fields would touch the validated numerical boundary without reducing an active duplicate or ambiguity.

## 4. Files changed

- `tests/testthat/test-blr-framework-phase9g.R`: permanent canonical architecture, policy, exact-reference, reproducibility, and protected-backend tests.
- `tools/benchmarks/blr_phase9g_csr_learned_annotation_bayesc.R`: canonical benchmark entry retaining the frozen Phase 9F3 workload.
- `docs/dev/blr_framework_implementation_plan.md`: canonical status and removal of stale unmigrated wording.
- `docs/dev/blr_model_capability_matrix.md`: canonical learned-annotation row.
- This report: inventory, validation, performance, and readiness record.
- Native source, typed headers, prior tests, public/generated files, and fixtures were not changed.

## 5. Canonical execution path

Public R validation/alignment prepares the native annotation policy and CSR state, constructs `CsrLearnedAnnotationBayesCExecutionContext`, calls canonical `run_csr_learned_annotation_bayesc()`, aggregates native chain payloads once, calls `stblr_csr_learned_annotation_bayesc_result_to_raw()`, validates unchanged `stblr_raw_v1`, and formats through the canonical R raw-to-fit path.

## 6. Cleanup and naming

Stable architecture names were retained. Stale documentation was renamed from migrated/ready to canonical. The canonical benchmark name replaces post-migration wording for future comparisons. Historical phase tests/reports and the Phase 9F3 benchmark remain as audit records. No dead native alias, unused field, duplicate converter, source-hash expectation, fixture-generation hook, or fallback was found.

## 7. Numerical core and implementation header

The native core was deliberately unchanged. Its conventional guard, translation-unit inclusion safeguard, inline callable definition, Armadillo/compiler configuration, binding neutrality, and sole active `isort` marker loop are permanently tested. No allocation was added inside the marker loop.

## 8. Native adapter

The adapter remains limited to R decoding, alignment, annotation/coefficient/prior/proposal/bound preparation, CSR/native-state preparation, typed-context construction, core invocation, wrapper aggregation, converter invocation, and exception translation. It contains no numerical sampling policy logic.

## 9. Result converter

`stblr_csr_learned_annotation_bayesc_result_to_raw()` remains the sole learned-annotation converter. It is binding-only and preserves field order, R types/dimensions/classes, actual `NULL`, marker/trait/annotation/coefficient/chain order, schema, annotation summaries, proposal/acceptance diagnostics, variance/global outputs, LD-swap, selection, timing, and optional chain payloads.

## 10. Wrapper-level multichain aggregation

Each core call returns one native payload. The wrapper aggregates marker/coefficient summaries and LD-swap counts once, in chain order, before conversion; effective probabilities/multipliers are not recomputed. Optional chain payloads preserve invocation order and are emitted only under existing `keep_chains` rules. This remains outside the core because it is schema- and chain-payload-specific.

## 11. Learned-annotation policy

Column-major marker-by-annotation layout, annotation order, explicit intercept handling, `K x nt` coefficient orientation/order/initialization, centered-logistic probabilities and centering, exponential multipliers, priors, random-walk proposals, flags/frequency, bounds/clipping order, proposal/acceptance diagnostics, global `pi`/`B`, and validation behavior are unchanged.

## 12. Ownership

CSR storage, marker/trait statistics, annotation matrix, marker order, priors, initial coefficients, and fixed proposal inputs are borrowed immutable resources that outlive the core call. Effects, inclusion states, residuals, variances, coefficient state, global parameters, RNG, accumulators, diagnostics, and workspaces remain chain-owned.

## 13. Logging

Retained `std::cout` messages are outside the marker loop and unsafe R APIs. They consume no RNG, do not alter control flow or returned state, and worker regions do not write per-marker diagnostics concurrently. Post-region trait reporting preserves existing diagnostic content.

## 14. Permanent regression fixtures

The three permanent configurations are `annot_fixed`, `annot_learned`, and `annot_explicit`. Together they protect fixed/learned coefficients, multiple annotations and explicit intercept behavior, both learning links, update frequency/bounds, diagnostics, core counts/order, chains/retention, disabled LD-swap, selection/schema fields, and validation/unsupported cases. Fixture generation remains manual and is not invoked by tests.

## 15. Exact reference results

Learned-annotation raw references: 3/3 exact. Learned-annotation formatted references: 3/3 exact.

## 16. Reproducibility

Exact tests pass for repeated calls, one/two cores, reversed `1,2,2,1` ordering, intervening group/annotation-aware fits, fixture explicit seeds, retained/dropped chains, fixed and learned modes, explicit update frequency/bounds, and disabled LD-swap. Only documented core-count metadata is normalized.

## 17. Public API and schema

Arguments, defaults, native signatures, public routing, `NAMESPACE`, generated wrappers, `stblr_raw_v1`, actual-`NULL` behavior, and formatted-fit schema remain byte- or fixture-exact.

## 18. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior BayesC, group BayesC, block-eigen, BED, scheduled CSR/BED, and multivariate files are unchanged from `1b4a0e4`. Permanent hashes and the full suite protect their schemas, seeds, core counts, LD-swap, selection, group, and annotation behavior.

## 19. Performance and memory baseline

Commands: `Rscript tools/benchmarks/blr_phase9f3_csr_learned_annotation_bayesc.R` and `Rscript tools/benchmarks/blr_phase9g_csr_learned_annotation_bayesc.R`. The canonical workload uses 40-marker smoke and 2,000-marker primary fits, one trait, three annotations/coefficients, 100 iterations after 25 burn-in, four learning modes, one/two chains and cores, retained chains, update frequency 5, probability bounds `[1e-8,0.5]`, and multiplier bounds `[1e-3,1e3]`. Five primary repetitions follow a warm-up.

Phase 9G means were 0.418 s fixed, 0.510 s probability-only, 0.492 s multiplier-only, and 0.562 s both-learning at one chain/core; both-learning was 1.120 s for two chains/one core and 1.406 s for retained two chains/two cores. The latter ranged 1.24--1.64 s and reflects noisy Windows scheduling. Phase 9F3 comparable means were 0.562, 1.122, and 1.102 s; no native code changed, the first two match closely, and the variable 2x2 result has no accompanying memory or correctness regression. Completed-fit RSS was 128.6--135.8 MiB, comparable with Phase 9F3. RSS is whole-process completed-fit sampling, not interval peak; interval sampling was unavailable. No unexplained material regression remains.

## 20. Test results

- Phase 9A and Phase 9F1--9F3 remain in the full suite.
- New Phase 9G: 35 passes, 0 failures, warnings, or skips.
- Baseline full suite: 4,327 passes.
- Final fresh-build full suite: 4,362 passes, 0 failures, warnings, or skips.
- Canonical-model and protected-backend checks pass through permanent hashes and focused/full tests.

## 21. Deviations and blockers

Interval peak-memory sampling was unavailable, so completed-fit process RSS is reported accurately. Windows timing, particularly retained two-chain/two-core execution, was noisy; no code or memory change accompanies it. No blocker remains.

## 22. Recommended next phase

> begin a scheduled-CSR contract, deterministic-reference, and RNG-ownership audit, including explicit investigation of stateful distribution lifetime and the previously identified `static thread_local` risk, without migrating scheduled execution yet.

## 23. Readiness marker

PHASE 9G COMPLETE — LEARNED-ANNOTATION BAYESC CANONICALIZED AND STABILIZED
