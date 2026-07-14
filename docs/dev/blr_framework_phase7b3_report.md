# Unified BLR Framework: Phase 7B3 Report

## 1. Executive summary

Ordinary-CSR SBayesRC migration is complete with one typed operator-aware numerical core and one named binding-layer result converter. Exact statistical, annotation, alpha, RNG, public-schema, and protected-backend behavior is preserved.

## 2. Repository baseline

Branch `master`, starting commit `26c7a3c` (Phase 7B2), initially clean. Toolchain: R 4.4.1, Rtools44/GCC 13.2, C++17, and OpenMP. The baseline full suite passed 3,801 assertions with zero failures, warnings, or skips. The directly rerun Phase 7A moderate benchmark produced 0.03–0.06 second medians and approximately 128–138 MiB whole-process RSS; these short runs were explicitly treated as noisy.

## 3. Phase 7B2 structure inventory

| Item | Classification | Phase 7B3 disposition |
|---|---|---|
| typed context and callable core | retain permanently | unchanged |
| typed execution result | retain permanently | unchanged |
| result aliases | centralize now | moved inside the converter and reduced to consumed fields |
| marker/trace/diagonal lambdas | centralize now | owned by the named converter |
| component, selection, chain, and raw-list construction | centralize now | one converter body |
| `Rcpp::Rcout` migration wording | rename/document | `std::cout` decision recorded as permanent binding-neutral logging |
| byte-identical production assumptions | remove now | permanent structural and numeric protections used instead |
| fixture generation | retain outside tests only | permanent RDS references unchanged |
| binding-derived summary calculations | retain permanently | preserved operation-for-operation inside conversion |

No old lexical execution path, comparison route, fallback, or old/new selector existed.

## 4. Files changed

- `src/st_sbayesrc_omp_csr.cpp`: adds binding metadata and the single named result converter; the generic native adapter now ends in one converter call.
- Phase 7B1/7B2 tests: replace inline-converter wording with permanent named-converter protection.
- `tests/testthat/test-blr-framework-phase7b3.R`: adds closure regressions and structural protections.
- `tools/benchmarks/blr_phase7b3_csr_sbayesrc.R`: establishes the five-repeat post-migration baseline.
- Framework plan and capability matrix: mark ordinary CSR SBayesRC migrated and ready for canonicalization.
- This report records implementation and validation.

## 5. Final execution path

`stblr_csr()` validation/alignment → existing SBayesRC native decoding and Armadillo/operator preparation → `CsrSBayesRCExecutionContext<Operator>` → `run_csr_sbayesrc()` → `CsrSBayesRCExecutionResult` → `stblr_csr_sbayesrc_result_to_raw()` → unchanged `stblr_raw_v1` → unchanged formatter.

## 6. Centralized result converter

`stblr_csr_sbayesrc_result_to_raw()` accepts the typed result plus a small binding-only metadata view. It is the sole ordinary-CSR SBayesRC converter. Marker and trace orientation, diagonal covariance conversion, component and step names, chain ordering, selection summaries, diagnostics, optional `NULL`, schema metadata, class vector, and all field ordering remain identical. It invokes no sampler and consumes no RNG.

## 7. Native adapter boundary

The generic native adapter retains argument decoding, validation/alignment, Armadillo preparation, operator construction, Phase 7A contract/context construction, the single core call, and the single converter call. It contains no active MCMC loop, marker loop, alpha update, assignment draw, variance update, posterior accumulation, or aggregation formula.

## 8. Numerical core

One `run_csr_sbayesrc<Operator>()` and one active MCMC/marker loop remain in `blr_csr_sbayesrc_core_impl.h`. No numerical statement was changed in Phase 7B3.

## 9. Ordered-probit and alpha preservation

Annotation order, annotations-by-sticks alpha orientation, stick/component order, `Phi(A alpha_j)`, successive remaining mass, probability flooring, row normalization, fixed/learned alpha behavior, and `alpha_update_every` timing are unchanged and protected by exact references.

## 10. Operator sharing

The existing ordinary `CsrOperator` and `BlockEigenOperator` factories continue to instantiate the same templated core through the generic native implementation. Block-eigen source and its result conversion were not modified; it does not call the ordinary-CSR converter independently.

## 11. Logging decision

The three Phase 7B2 `std::cout` statements remain at the same pre-parallel/post-parallel diagnostic positions formerly occupied by `Rcpp::Rcout`. They are outside the marker hot path, do not touch statistical state or RNG, preserve message content and timing intent, and avoid Rcpp in the reusable core. Worker threads do not write console output.

## 12. Exact frozen references

Raw references: 6/6 exact. Formatted references: 6/6 exact. Values, storage types, dimensions, names, dimnames, classes, actual `NULL`, marker/trait/annotation/component/chain order, schema, alpha summaries, component probabilities, and diagnostics match.

## 13. Reproducibility

Repeated calls, one/two cores with only documented `input$ncores` normalization, reversed 1-2-2-1 order, an intervening canonical BayesR fit, explicit chain seeds, multiple traits, retained/dropped chains, fixed alpha, learned alpha, and explicit update frequency are exact.

## 14. Annotation and alpha behavior

Annotation column/intercept semantics, alpha orientation/dimensions and posterior summaries, ordered-probit marker probabilities, component counts and summaries, disabled LD-swap, and selection fields are unchanged. No alpha trace or acceptance counter was invented.

## 15. Public API and schema

Public arguments/defaults, routing, native signatures, `NAMESPACE`, generated wrappers, `stblr_raw_v1`, formatted fit, optional `NULL`, and one-marker/one-trait conventions remain unchanged.

## 16. Protected backends

Canonical BayesC and BayesR source hashes and permanent references pass. Block-eigen sources and focused tests pass. BED BayesRC, fixed-prior, group, learned-annotation, scheduled, BED, and multivariate sources are unchanged.

## 17. Performance and memory

Command: `Rscript tools/benchmarks/blr_phase7b3_csr_sbayesrc.R`. The stable workload used 2,000 markers, one trait, four annotations, four components, 12 alpha parameters, 120 iterations, and five timed repetitions after warm-up.

| Configuration | Times (s) | Mean | Median | Range | RSS MiB |
|---|---|---:|---:|---:|---:|
| fixed, 1 chain/1 core | 0.24, 0.24, 0.21, 0.24, 0.22 | 0.230 | 0.24 | 0.03 | 140.5 |
| learned, 1 chain/1 core | 0.28, 0.28, 0.26, 0.26, 0.28 | 0.272 | 0.28 | 0.02 | 131.0 |
| learned, 2 chains/1 core | 0.55, 0.52, 0.55, 0.55, 0.56 | 0.546 | 0.55 | 0.04 | 131.1 |
| learned, 2 chains/2 cores | 0.33, 0.27, 0.29, 0.30, 0.30 | 0.298 | 0.30 | 0.06 | 130.5 |
| learned, 2 chains/2 cores, retained | 0.26, 0.28, 0.27, 0.28, 0.30 | 0.278 | 0.28 | 0.04 | 131.8 |

RSS is whole-process memory sampled after each completed fit, not sampler-only interval peak. The current public route has no minimal-output mode, so ordinary and retained-chain modes were measured. Phase 7A's 500-marker timings are too short for direct percentage claims; operation order and exact references are the primary migration evidence. The larger baseline shows stable scaling, no unexplained runtime regression, and no meaningful memory increase.

## 18. Test results

Baseline full suite: 3,801 passed. Phase 7A: 82; Phase 7B1: 25; Phase 7B2: 49; Phase 7B3: 55. The final focused SBayesRC, annotation/component, canonical BayesC/BayesR, block-eigen, schema, backend, LD-swap, and routing suite passed 2,100 assertions. The final full suite passed 3,859 assertions. Failures, warnings, and skips: zero.

## 19. Deviations and blockers

The public SBayesRC route does not expose a minimal-output control, so that benchmark mode is documented as unsupported rather than invented. Memory is sampled at fit completion because sampler-interval peak instrumentation is unavailable. The converter retains existing binding-derived summaries to preserve the public schema exactly. No blocker remains.

## 20. Recommended Phase 8 boundary

> canonicalize and stabilize ordinary-CSR SBayesRC, remove remaining migration-only wording or aliases, retain permanent exact references, and establish it as the canonical annotation-aware mixture implementation before migrating adjacent annotation backends.

## 21. Readiness marker

PHASE 7B3 COMPLETE — CSR SBAYESRC MIGRATED WITH BEHAVIOR PRESERVED
