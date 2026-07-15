# Unified BLR Framework: Phase 9B1 Report

## 1. Executive summary

The fixed-prior CSR BayesC execution block was mechanically extracted into one guarded implementation header without changing its numerical statements, production route, or inline R result conversion.

## 2. Repository baseline

Phase 9B1 started on branch `master` at commit `e89c920` (`Add annotation BayesC contracts and migration references`) with a clean working tree. This is the committed Phase 9A baseline. The environment was R 4.4.1 with Rtools44/GCC 13.2, C++17, and OpenMP. Native compile/load completed successfully. The pre-edit focused suite passed 298 assertions with zero failures, warnings, or skips; the committed Phase 9A full-suite baseline was 3,984 assertions.

## 3. Extraction seam

The relocated span begins after input validation, Armadillo preparation, fixed `pi_marker` and `vb_multiplier` preparation, shared CSR construction, LD-swap friend preparation, and marker-order construction. Its first statement allocates the native output matrices under `// Output storage`. It ends with `return result;` in `stblr_cpg_omp_csr_prior_single()`, immediately before `cpg_raw_marker_matrix()` starts the binding conversion helpers.

The production source now includes `blr_csr_prior_bayesc_core_impl.h` at exactly that lexical position.

## 4. Files changed

- `src/st_cpg_omp_csr_prior.cpp`: replaced the 542-line execution span with one guarded implementation-header include.
- `src/blr_csr_prior_bayesc_core_impl.h`: contains the mechanically relocated execution span and an include guard.
- `tests/testthat/test-blr-framework-phase9a.R`: replaced the intentionally obsolete fixed-prior production-source MD5 assertion while retaining byte protection for the still-unmigrated group and learned-annotation sources.
- `tests/testthat/test-blr-framework-phase7a.R`, `test-blr-framework-phase7b2.R`, `test-blr-framework-phase7b3.R`, and `test-blr-framework-phase8.R`: removed only the now-obsolete fixed-prior source entry from historical protected-file MD5 lists; every other historical hash remains protected.
- `tests/testthat/test-blr-framework-phase9b1.R`: adds permanent extraction, single-loop, inclusion-site, fixed-policy, exact-reference, and reproducibility protections.
- `docs/dev/blr_framework_phase9b1_report.md`: records the extraction and validation.

## 5. Extracted execution content

The header contains native output and failure storage, OpenMP trait execution, chain-local RNG and sampler state, residual/effect initialization, the MCMC loop, marker traversal, fixed marker probability and variance-multiplier access, optional global and variance updates, LD-swap attempts, posterior accumulation, timing and failure diagnostics, and construction of the native 24-slot execution result.

The existing public wrapper's chain-seed resolution, repeated invocation of the single-chain execution function, cross-chain summary construction, and optional chain payload handling remain byte-for-byte in the production translation unit. They were not part of the contiguous numerical span between operator preparation and the first conversion helper.

## 6. Fixed-prior preservation

`pi_marker_t(i)` and `vb_multiplier_t(i)` are accessed in the same marker order and at the same statements. Their fallbacks to global `pi_t[1]` and multiplier `1.0`, respectively, are unchanged. `updatePi`, `samplePi_ST_prior`, `updateB`, `sampleB_ST_csr_prior`, global `pi`, global `B`, residual and variance updates, RNG construction and calls, OpenMP static scheduling, retained-sample timing, LD-swap, selection output, accumulation, and result population were relocated without reordering or rewriting.

## 7. Dimension limitation

The existing fixed-prior production trait-dimension behavior is unchanged. In particular, the unsupported multiple-trait construction documented in Phase 9A was not broadened and no validation was relaxed.

## 8. Existing R conversion

`cpg_prior_raw_v1()` and its matrix/trace/chain conversion helpers remain inline in `src/st_cpg_omp_csr_prior.cpp`. Field construction, schema, actual `NULL` values, classes, ordering, diagnostics, LD-swap fields, selection fields, and optional chain payloads were not modified.

## 9. Frozen references

All three fixed-prior raw fixtures matched exactly (3/3), and all three formatted fixtures matched exactly (3/3). Comparisons use `identical()` after only the permanent deterministic timing/path normalization defined by Phase 9A; references were not regenerated.

## 10. Reproducibility

Repeated fixed-seed calls were exact. One-core and two-core results were exact after normalizing only documented `input$ncores` metadata. Reversed `1, 2, 2, 1` execution order was exact. An intervening learned-annotation fit did not change the result. The frozen configurations retain exact explicit chain seeds, retained and dropped chain behavior, fixed `pi_marker`, fixed `vb_multiplier`, both inputs together, disabled global update flags, and disabled LD-swap behavior.

## 11. Protected backends

Git comparison against `e89c920` confirmed no changes to canonical CSR BayesC, BayesR, or SBayesRC; group or learned-annotation CSR BayesC; block-eigen; BED; scheduled; or multivariate production sources. `NAMESPACE`, generated wrappers, public signatures, and public R files are unchanged.

## 12. Tests

The pre-edit focused baseline passed 298 assertions. Phase 9A passed 67 assertions after replacing only the obsolete fixed-prior MD5 assumption. New Phase 9B1 tests passed 28 assertions. Native compilation succeeded. The final focused and full-suite results are recorded after final validation below:

- Phase 9A: 67 passed, 0 failed, 0 warned, 0 skipped.
- Phase 9B1: 28 passed, 0 failed, 0 warned, 0 skipped.
- Focused fixed-prior/protection suite: 1,177 passed, 0 failed, 0 warned, 0 skipped.
- Historical framework guard rerun: 339 passed, 0 failed, 0 warned, 0 skipped.
- Full suite: 4,012 passed, 0 failed, 0 warned, 0 skipped.

## 13. Deviations and blockers

The numerical span is 542 lines, rather than the broader wrapper-level chain orchestration. This follows the actual contiguous source seam: wrapper-level multichain aggregation is separated from the sampler by binding conversion helper definitions and remains unchanged for Phase 9B2. The first full run identified four historical MD5 assertions that still named the intentionally changed fixed-prior source; those obsolete entries were removed while all other hashes and the new structural protection were retained. The subsequent full suite passed. No numerical, compiler, reference, or ownership blocker was found. The implementation header remains a deliberately lexically dependent block for this subphase.

## 14. Recommended Phase 9B2 task

> replace the lexically dependent fixed-prior execution include with an explicit typed execution context and callable core using the Phase 9A contracts, while retaining the existing inline R result conversion until all frozen references pass again.

## 15. Readiness marker

PHASE 9B1 COMPLETE — FIXED-PRIOR BAYESC EXECUTION BLOCK EXTRACTED
