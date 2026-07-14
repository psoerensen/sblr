# Unified BLR Framework: Phase 8 Report

## 1. Executive summary

Ordinary-CSR SBayesRC is canonicalized and stabilized. The typed borrowed context, single operator-templated numerical core, typed result, and sole ordinary-CSR binding converter remain active, with no statistical code change.

## 2. Repository baseline

Branch `master`, starting commit `8d4d7c2` (Phase 7B3), initially clean. Toolchain: R 4.4.1, Rtools44/GCC 13.2, C++17, and OpenMP. Baseline compilation succeeded; the full suite passed 3,859 assertions and the focused suite passed 2,100, with no failures, warnings, or skips. The rerun Phase 7B3 benchmark medians were 0.24, 0.28, 0.54, 0.28, and 0.28 seconds; completed-fit RSS was 130.6–138.0 MiB.

## 3. Phase 7B3 structure inventory

| Item | Classification | Disposition |
|---|---|---|
| typed context/core/result/converter | retain permanently | unchanged |
| Phase-numbered test and fixture names | retain permanently | historical provenance and stable fixture lookup |
| fixture generator | retain permanently | manual maintenance tool; never run by tests |
| `legacy` wording in active converter comment | rename for canonical use | changed to model-neutral wording |
| conversion aliases | retain permanently | scoped inside the sole converter and express stable R orientation |
| context/result fields | retain permanently | all support execution, conversion, diagnostics, or documented typed vocabulary |
| repeated historical structural tests | retain | migration checkpoints remain useful; Phase 8 adds the canonical contract |
| commented pre-framework prototypes | defer | inactive archival source, not a route or fallback; deleting it would be unrelated broad cleanup |
| Phase 7B3 benchmark wording | consolidate | Phase 8 benchmark records the canonical baseline and names RSS accurately |

No active old/new selector, fallback, duplicate converter, or duplicate numerical loop was found.

## 4. Files changed

- `src/st_sbayesrc_omp_csr.cpp`: removes stale active-path “legacy backend” wording only.
- `tests/testthat/test-blr-framework-phase8.R`: permanent canonical structural, exact-reference, reproducibility, ownership, and protected-file tests.
- `tools/benchmarks/blr_phase8_csr_sbayesrc.R`: canonical five-repeat timing and completed-fit RSS benchmark.
- Framework plan and capability matrix: mark ordinary CSR SBayesRC canonical.
- This report records stabilization and validation.

## 5. Canonical execution path

Public R validation/alignment → native decoding and Armadillo/operator preparation → `CsrSBayesRCExecutionContext<Operator>` → `run_csr_sbayesrc()` → `CsrSBayesRCExecutionResult` → `stblr_csr_sbayesrc_result_to_raw()` → unchanged `stblr_raw_v1` → unchanged formatted fit.

## 6. Cleanup and naming

The canonical type/function names were already stable and retained. One stale “legacy backend” comment was renamed. Historical phase test/fixture names are retained as provenance; no migration alias exists in the active route. Inactive commented prototypes were deferred because they neither compile nor provide a fallback and removing thousands of archival lines would exceed this bounded phase.

## 7. Numerical core and implementation header

`blr_csr_sbayesrc_core_impl.h` has a conventional guard, is included once from `st_sbayesrc_omp_csr.cpp`, inherits that translation unit's Armadillo configuration, constructs no R objects, and contains the sole `run_csr_sbayesrc<Operator>()` and active MCMC/marker loop. Both `CsrOperator` and `BlockEigenOperator` use that template. No core statement changed.

## 8. Native adapter

The adapter decodes R inputs, aligns annotations/markers, validates, prepares Armadillo/operator state, builds the typed context, calls the core, and calls the converter. It contains no active MCMC, marker, alpha/probability update, assignment draw, variance update, accumulation, or aggregation loop.

## 9. Result converter

`stblr_csr_sbayesrc_result_to_raw()` remains the sole ordinary-CSR converter. It preserves names/order, storage, dimensions, dimnames, classes, actual `NULL`, marker/trait/annotation/component/chain ordering, schema metadata, alpha and component summaries, variances, VLE/VLD, diagnostics, timing, failure data, LD-swap, selection, optional chains, and input metadata.

## 10. Ordered-probit and alpha policy

Annotations remain marker-by-column in their supplied order with the established intercept convention. Alpha remains annotations-by-modeled-stick; component/null and stick order are unchanged. `Phi(A alpha_j)`, successive remaining mass, probability floor, row normalization, fixed/learned modes, priors, and `alpha_update_every` timing are unchanged.

## 11. Ownership

Operator storage, annotations, prepared statistics, component scales, priors, initialization, LD-swap friends, and metadata are borrowed immutable resources whose lifetime exceeds execution. Effects, assignments, residuals, variances, mutable probabilities, RNG engine/distributions, accumulators, diagnostics, and workspace are task/chain-owned. No per-chain CSR or annotation payload exists.

## 12. Logging

The three retained `std::cout` diagnostics are outside the marker loop, do not consume RNG or affect returned state/control flow, and preserve the prior messages. Configuration messages occur before the parallel region; per-task completion messages occur after task state is complete. No worker marker loop writes output.

## 13. Permanent regression fixtures

Six compact fixtures remain unchanged: fixed one-chain, fixed two-chain, learned two-core, learned explicit seeds with retained chains/update frequency, multiple traits, and explicit component scales. Together they cover fixed/learned alpha, one/two chains and cores, seeds, traits, retained/dropped chains, annotation order/intercept, scales, probabilities/counts/summaries, disabled LD-swap, and selection fields.

## 14. Exact reference results

Raw references: 6/6 exact. Formatted references: 6/6 exact. Comparisons cover values, storage, dimensions, names/dimnames, classes, actual `NULL`, ordering, schema, alpha summaries, probabilities, counts, and diagnostics.

## 15. Reproducibility

Repeated calls, one/two cores with only documented `input$ncores` normalization, reversed 1-2-2-1 ordering, an intervening canonical BayesR fit, explicit seeds, multiple traits, retained/dropped chains, fixed/learned alpha, and explicit update frequency remain exact.

## 16. Public API and schema

Arguments/defaults, native signatures, routing, `NAMESPACE`, generated wrappers, `stblr_raw_v1`, formatted fit, present `NULL`, and scalar matrix conventions are unchanged.

## 17. Protected backends

MD5 and focused tests protect canonical BayesC/BayesR, block-eigen sources, BED BayesRC, fixed-prior/group/learned-annotation BayesC, and `NAMESPACE`. Git comparison also confirms scheduled, BED, multivariate, and other non-SBayesRC backends are untouched. Block-eigen continues template instantiation but does not use the ordinary-CSR converter independently.

## 18. Performance and memory baseline

`Rscript tools/benchmarks/blr_phase8_csr_sbayesrc.R` used 2,000 markers, one trait, four annotations/components, 12 alpha parameters, 120 iterations, warm-up, and five repetitions. Median seconds for fixed 1×1, learned 1×1, learned 2×1, learned 2×2, and learned 2×2 retained-chain runs were 0.24, 0.28, 0.53, 0.27, and 0.28. Completed-fit RSS was 130.8–140.5 MiB. Phase 7B3 medians were 0.24, 0.28, 0.54, 0.28, 0.28 with RSS 130.6–138.0 MiB; differences are immaterial short-run/process variability. RSS is sampled after each completed fit, not an interval peak or sampler-only measure. No material runtime or memory regression is evident.

## 19. Test results

- Baseline full suite: 3,859 passed; 0 failed/warned/skipped.
- Baseline focused Phase 1–7/backend suite: 2,100 passed; 0 failed/warned/skipped.
- Phase 8 suite: 58 assertions after correction of one overly specific structural token; 0 final failures/warnings/skips.
- Final combined focused suite: 2,158 passed; 0 failed/warned/skipped.
- Final full suite: 3,917 passed; 0 failed/warned/skipped.

## 20. Deviations and blockers

No blocker exists. Compiler identity is known from the Rtools build (GCC 13.2), although `R.version$compiler` is unavailable. Memory is accurately reported as completed-fit whole-process RSS because interval peak sampling was unavailable. Minimal output is not supported by the public SBayesRC route. Archival commented prototypes were intentionally deferred as non-active historical material.

## 21. Recommended next phase

Perform a comparative contract-and-reference audit of fixed-prior, group, and learned-annotation BayesC backends to determine which probability, scale, annotation, and result infrastructure is genuinely shared before migrating any of them.

## 22. Readiness marker

PHASE 8 COMPLETE — CSR SBAYESRC CANONICALIZED AND STABILIZED
