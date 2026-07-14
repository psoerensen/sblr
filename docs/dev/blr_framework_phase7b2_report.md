# Unified BLR Framework: Phase 7B2 Report

## 1. Executive summary

The lexically dependent shared SBayesRC execution block was converted into an explicit operator-templated callable core. Ordinary CSR and block-eigen continue to instantiate one numerical body; the existing inline Rcpp conversion remains in place.

## 2. Repository baseline

Branch `master`, starting commit `a7acf74` (Phase 7B1), initially clean. The toolchain is R 4.4.1 with Rtools44/GCC 13.2, C++17, and OpenMP. Phase 7B1 had 6/6 raw and 6/6 formatted references exact and a 3,752-assertion full suite.

## 3. Dependency inventory

| Dependency | Dimensions/ownership | Explicit representation |
|---|---|---|
| concrete operator and LD friends | shared borrowed, execution lifetime | `op`, `ld_swap_friends` references |
| marker/trait prepared statistics | traits × markers, shared borrowed | `wy_mat`, `ww_mat`, `yy_vec`, sample-size references |
| initial effects/residuals/components | traits × markers; prepared mutable aggregate plus borrowed optional inputs | matrix and nested-vector references |
| annotation design | markers × annotations, column-major, shared borrowed | `A` plus `SBayesRCAnnotationDesignView` |
| components | component ordered scales, null at zero | `gamma` plus `SBayesRCComponentSpec` |
| alpha | annotations × sticks and one variance per stick | native references plus `SBayesRCAlphaSpec` |
| ordered-probit policy | component stick order and probability floor | `SBayesRCProbabilityPolicy` |
| priors | trait covariance diagonals and scalar degrees of freedom | borrowed matrices and `SBayesRCPriors` |
| tasks/seeds/MCMC | scalar controls and borrowed explicit seeds | context scalars and `SBayesRCControls` |
| LD swap/selection | borrowed friend/ranking/scaling vectors and scalar controls | explicit context fields |
| task and aggregate outputs | task/trait × marker/trace/component layouts | `CsrSBayesRCExecutionResult` |

All borrowed resources are constructed above the seam and outlive the call. Effects, assignments, residuals, alpha state, probabilities, variances, RNG, accumulators, and scratch remain chain-local inside the task loop.

## 4. Files changed

- `src/blr_csr_sbayesrc_core_impl.h`: adds the templated context, typed execution result, callable function, aliases, and result population around the unchanged numerical block.
- `src/st_sbayesrc_omp_csr.cpp`: constructs Phase 7A contract views and the operator-aware context, calls the core, and exposes result fields to the unchanged inline converter.
- `tests/testthat/test-blr-framework-phase7b1.R`: extends the earlier structural protection to recognize the callable boundary.
- `tests/testthat/test-blr-framework-phase7b2.R`: adds permanent structural, binding-neutral, exact-reference, reproducibility, and protected-file checks.
- `docs/dev/blr_framework_phase7b2_report.md`: records this bounded migration.

## 5. Execution context

`CsrSBayesRCExecutionContext<Operator>` holds the concrete operator and all large prepared matrices/vectors by reference. Small dimensions, flags, iteration controls, update frequencies, probability floor, and scalar priors are copied. The Phase 7A annotation, component, alpha, probability, prior, control, and output contracts are constructed at the approved seam and referenced by the context.

## 6. Typed result

`CsrSBayesRCExecutionResult` owns task-level marker, trace, component, alpha, timing, failure, LD-swap, and selection outputs plus trait-level aggregates, chain summaries, final variances/probabilities, posterior alpha/component summaries, and diagnostics. Layout and ordering match the former lexical locals consumed by the Rcpp converter.

## 7. Callable core

The callable signature is `template <class Operator> CsrSBayesRCExecutionResult run_csr_sbayesrc(CsrSBayesRCExecutionContext<Operator>&)`. The single generic binding implementation constructs the context, so the existing `CsrOperator` and `BlockEigenOperator` factories instantiate the same function at compile time.

## 8. Numerical preservation

The ordered-probit `Phi(A alpha_j)` calculation, alpha and stick order, remaining-mass calculation, flooring, normalization, fixed/learned alpha branches, update frequency, marker traversal, component sampling, RNG construction and calls, task mapping, static OpenMP scheduling, accumulation, aggregation, and diagnostics are unchanged. Three binding stream expressions were mechanically changed from `Rcpp::Rcout` to `std::cout` at the same positions; they do not affect numerical state.

## 9. Inline R conversion

The existing converter remains inline in `src/st_sbayesrc_omp_csr.cpp`. It reads references exposed from the typed result; field construction, dimensions, names, classes, optionality, schema, annotation/component/chain order, and diagnostics are unchanged. Converter centralization is deferred to Phase 7B3.

## 10. Exact references

Raw references: 6/6 exact. Formatted references: 6/6 exact.

## 11. Reproducibility

Repeated calls, one/two cores with only declared `input$ncores` normalization, reversed 1-2-2-1 order, explicit seeds, multiple traits, retained/dropped chains, fixed alpha, learned alpha, and explicit update frequency remain exact. The Phase 7A intervening canonical BayesR check also remains active.

## 12. Annotation and alpha behavior

Annotation order and intercept handling, annotations × sticks alpha orientation, fixed and learned alpha, posterior alpha summaries, update timing, marker-specific component probabilities, component counts/summaries, and disabled LD-swap behavior match the frozen references exactly. No unsupported alpha traces or acceptance diagnostics were invented.

## 13. Protected backends

Canonical BayesC/BayesR, block-eigen source files, BED BayesRC, fixed-prior/group/learned-annotation sources, generated wrappers, public signatures, schemas, and `NAMESPACE` remain unchanged.

## 14. Tests

Phase 7A passed 82 assertions, Phase 7B1 passed 24, and Phase 7B2 passed 47. The broader focused suite passed 1,910 assertions. The full suite passed 3,801 assertions. All runs had zero failures, warnings, or skips.

## 15. Deviations and blockers

The shared concrete-operator requirement is represented by a templated context rather than type erasure. Phase 7A's raw CSR view is not imposed on block-eigen; the concrete operator is the binding-neutral borrowed execution resource. No numerical blocker remains.

## 16. Recommended Phase 7B3 task

> centralize the typed SBayesRC execution-result-to-R converter, replace obsolete migration assertions with permanent structural tests, benchmark pre/post execution, update framework documentation, and complete final Phase 7B validation.

## 17. Readiness marker

PHASE 7B2 COMPLETE — SBAYESRC TYPED EXECUTION BOUNDARY ACTIVE
