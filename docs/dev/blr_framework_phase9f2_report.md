# Unified BLR Framework Phase 9F2 Report

## 1. Executive summary

The lexically dependent learned-annotation CSR BayesC execution block was converted into an explicit, binding-neutral typed context, callable core, and typed result. The existing numerical statements, inline R conversion, and wrapper-level multichain aggregation were retained.

## 2. Repository baseline

- Branch: `master`
- Starting and Phase 9F1 commit: `97da778` (`Extract learned-annotation BayesC execution block`)
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1 UCRT and GCC 13.2.0 (Rtools44).
- Phase 9F1 checkpoint: full suite 4,254 passes, no failures/errors; focused Phase 9A/9F1 suite 103 passes. A fresh Phase 9F2 native build and the same 103 checks passed before the new Phase 9F2 test was run.

## 3. Dependency inventory

| Names/category | Original type/dimensions | Mutability/ownership | Scope and lifetime | Contract/context/result |
|---|---|---|---|---|
| `ld`, CSR row pointers/indices/values | native `SpMat2`, `m` markers | immutable borrowed | shared; adapter outlives call | context `ld_storage`, `ld_row_ptr_count`; operator reconstructed as a reference alias |
| `ww`, `wy`, `yy`, `n` | marker/trait matrices and trait vectors | `ww`, `yy`, `n` borrowed immutable; prepared `wy` mutable | shared prepared storage | explicit matrix/vector pointers; no result field |
| `b`, `r`, `d`, initial states | `m x nt` Armadillo matrices | adapter-owned prepared state; chain-local copies remain in core | call and chain lifetime | explicit pointers; posterior raw slots in result |
| marker/trait dimensions and order | scalars and `vector<int>(m)` | copied controls; order borrowed immutable | call lifetime | explicit counts/order pointer |
| `A` and annotation metadata | column-major `m x K` `arma::mat` | borrowed immutable | shared; adapter outlives call | policy `annotation_values`, layout/order/intercept plus context matrix pointer |
| `eta_pi`, `eta_vb` initial values | `K x nt` matrices | borrowed immutable inputs; coefficient state copied per trait/chain | shared input/chain-owned state | policy coefficient pointers/counts and explicit context pointers; summaries in raw result |
| coefficient priors/proposals | four scalars | copied immutable controls | call lifetime | policy and context fields |
| learning flags/update frequency | booleans and integer | copied immutable controls | call lifetime | policy and context fields |
| probability/multiplier links, centering, bounds | policy tags and four bounds | copied metadata | call lifetime | active learned-annotation policy; centered-logistic and exponential tags |
| effective probabilities/multipliers | marker vectors | chain-owned mutable | recomputed at existing proposal timing | accumulated in typed raw result |
| global `pi`, `B`, update flags | trait vector/matrix and booleans | borrowed initial values; chain-local mutable state | call/chain lifetime | explicit context fields; raw result outputs |
| variance/residual priors and controls | matrices/scalars | borrowed immutable/copy controls | call lifetime | explicit context fields; traces in result |
| MCMC, seed, thread controls | scalar values | copied immutable controls | call lifetime | explicit context fields |
| LD-swap friends/ranking/control | native friend object, order, scalars | borrowed immutable/copy controls | call lifetime | explicit opaque native pointer and controls; raw output slots |
| accumulators, proposal/acceptance counts, timing, diagnostics | vectors/matrices/scalars | chain-owned mutable | created inside callable core | aggregated into the typed raw result |

Every former lexical dependency is now either a context field or local/chain-owned state created by the callable core. The result owns all 24 raw numerical slots consumed below the seam.

## 4. Files changed

- `src/blr_csr_learned_annotation_bayesc_types.h`: added context, result, validation, and callable declaration.
- `src/blr_csr_learned_annotation_bayesc_core_impl.h`: made the extracted block a callable implementation and added original-name aliases.
- `src/st_cpg_omp_csr_annot.cpp`: constructs the active policy/context, calls the core, and passes its raw result to the unchanged inline converter.
- `tests/testthat/test-blr-framework-phase9f2.R`: permanent architecture, exact-reference, reproducibility, and protected-file tests.
- Historical framework tests: removed only the obsolete learned-annotation production-source MD5 expectation; permanent structural checks replace it.
- This report records the boundary and validation.

## 5. Execution context

`CsrLearnedAnnotationBayesCExecutionContext` explicitly carries prepared operator storage, marker/trait data, mutable adapter-owned initial state, marker order, annotation matrix, initial coefficients, priors, proposal/bound controls, MCMC/seeding/thread controls, and LD-swap inputs. Large resources are borrowed; scalar controls are copied. All borrowed resources must outlive `run_csr_learned_annotation_bayesc()`.

## 6. Learned-annotation contract activation

`LearnedAnnotationBayesCPolicyView` is constructed at the native seam and validated by the core. It explicitly records column-major marker-by-annotation storage, ordinal annotation-column order, no implicit intercept, `K x nt` coefficient orientation, priors, proposal scales, learning flags, update frequency, centered-logistic probability link, exponential multiplier link, and clipping bounds. No annotation reordering, intercept insertion, or dimension broadening was made.

## 7. Callable core

```cpp
CsrLearnedAnnotationBayesCExecutionResult
run_csr_learned_annotation_bayesc(
    const CsrLearnedAnnotationBayesCExecutionContext& context);
```

The implementation header is included by only `st_cpg_omp_csr_annot.cpp` and contains the sole numerical implementation and marker loop.

## 8. Typed result

`CsrLearnedAnnotationBayesCExecutionResult` owns the existing 24-slot native raw payload: marker posterior/state summaries, variance and global summaries, VLE/VLD, LD-swap/selection data, effective annotation probabilities and multipliers, coefficient summaries, proposal/acceptance diagnostics, and metadata already consumed by conversion and aggregation. No public field or new trace was added.

## 9. Numerical preservation

Coefficient initialization, priors, proposal construction/timing, acceptance calculations and counters, centered-logistic centering, exponential transformation, clipping order, update frequency, global interactions, RNG construction/calls, task ordering, OpenMP scheduling, posterior accumulation, LD-swap, selection, and diagnostics are operation-for-operation unchanged. Only dependency aliases and the typed return envelope were added.

## 10. Wrapper-level aggregation

The callable still returns the native result for one wrapper invocation. Existing per-chain conversion, chain-order handling, coefficient/effective-marker summary aggregation, proposal/acceptance aggregation, and `keep_chains` behavior remain in the wrapper. There is one aggregation path; no result is converted or counted twice.

## 11. Inline R conversion

Rcpp conversion remains inline in `src/st_cpg_omp_csr_annot.cpp`. It consumes `execution_result.raw` with the existing field order, types, dimensions, names, `NULL` behavior, and schema.

## 12. Exact references

All Phase 9A learned-annotation fixtures passed exactly: raw 3/3 and formatted 3/3, including types, dimensions, names, annotation/coefficient order, effective values, proposal/acceptance fields, and diagnostics.

## 13. Reproducibility

Exact checks passed for repeated calls, one/two cores and reversed `1,2,2,1` execution order, intervening group fits, explicit seeds represented by the fixtures, retained/dropped chains, fixed/learned coefficient modes, update frequencies, probability/multiplier bounds, and disabled LD-swap. Only already-declared execution metadata is normalized.

## 14. Unsupported cases

Existing validation and limitations remain: annotation rows/columns and coefficient dimensions must align; annotations, proposal scales, update frequency, and bounds must be valid; current shared-`ww`/sample-size construction and supported trait dimensions were not broadened.

## 15. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior BayesC, group BayesC, block-eigen, BED, scheduled, and multivariate production files remain unchanged from `97da778`. Generated wrappers, `NAMESPACE`, public signatures, and schemas remain unchanged.

## 16. Tests

- Phase 9A plus Phase 9F1: 103 passes, 0 failures, 0 warnings, 0 skips.
- Phase 9F2: 39 passes, 0 failures, 0 warnings, 0 skips.
- Full suite after a fresh compile: 4,294 passes, 0 failures, 0 warnings, 0 skips.
- The focused Phase 9F2 suite is included in that full result and independently passed 39/39.

## 17. Deviations and blockers

The pre-edit baseline is the clean committed Phase 9F1 validation recorded in its report; the fresh native build was performed after the mechanical typed-boundary edit. `std::cout` replaces the former binding-specific diagnostic stream inside the binding-neutral core; messages remain outside the marker loop and do not affect RNG or returned state. No blocker remains.

## 18. Recommended Phase 9F3 task

> centralize the typed learned-annotation BayesC result-to-R converter, remove obsolete migration assertions and aliases, establish post-migration runtime and memory baselines, update framework documentation, and close the learned-annotation migration.

## 19. Readiness marker

PHASE 9F2 COMPLETE — LEARNED-ANNOTATION BAYESC TYPED EXECUTION BOUNDARY ACTIVE
