# Unified BLR Framework Phase 9F1 Report

## 1. Executive summary

The learned-annotation CSR BayesC execution block was mechanically extracted
into one guarded implementation header without changing behavior. The 646
moved lines are byte-for-byte equivalent apart from their new guarded location.

## 2. Repository baseline

- Branch: `master`.
- Starting and Phase 9E commit: `0cf12d1` (`Canonicalize group CSR BayesC`).
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1 UCRT, Rtools44 GCC, C++17 and OpenMP on x86_64 Windows.
- Baseline native compile/load passed.
- Baseline full suite: 4,218 passed, 0 failed, 0 warned, 0 skipped.
- Baseline learned-annotation/framework focused tests passed.

## 3. Extraction seam

The moved block begins immediately after CSR construction, optional LD-swap
friend construction, marker ranking and `std::sort(order...)`. Its first line is
`arma::mat bm_mat(nt, m, arma::fill::zeros);`. It ends with `return result;`
immediately before the native single-run function closes and before
`cpg_raw_marker_matrix()` starts the binding conversion helpers.

## 4. Files changed

- `src/st_cpg_omp_csr_annot.cpp`: replaced the contiguous 646-line execution
  body with the guarded include at the same lexical seam.
- `src/blr_csr_learned_annotation_bayesc_core_impl.h`: contains the exact moved
  body plus include/translation-unit guards and an implementation-detail comment.
- `tests/testthat/test-blr-framework-phase9a.R`: replaced the intentionally
  obsolete learned-annotation source MD5 with structural route protections.
- `tests/testthat/test-blr-framework-phase9f1.R`: permanent extraction,
  policy, exact-reference, reproducibility and protected-backend checks.
- Historical Phase 7A, 7B2, 7B3, 8, 9B3, 9C, 9D1 and 9E tests: refreshed only
  the intentionally changed learned-annotation production-source hash; all
  other protected hashes remain unchanged. The new structural tests supersede
  that byte-level assumption for the migration boundary.
- This report.

## 5. Extracted execution content

The header now contains output/task allocation, OpenMP trait execution,
chain-local RNG and sampler state, coefficient state and proposals,
probability/multiplier updates, MCMC and marker loops, posterior accumulation,
acceptance/proposal diagnostics, trait summaries, native result assembly and
single-run aggregation.

## 6. Learned-annotation policy preservation

Annotation dimensions/layout/order, intercept convention, coefficient
dimensions/order/initialization, centered-logistic inclusion probabilities,
exponential variance multipliers, priors, proposal distributions/scales,
update frequency, probability and multiplier bounds/clipping, acceptance and
proposal counters, effective marker values, global `pi`/`B` interactions, RNG
construction/invocation order, static OpenMP scheduling, retained-sample timing,
LD-swap and selection fields are unchanged.

## 7. Wrapper-level aggregation boundary

The extracted single-run block returns one native result. Public multichain
orchestration and its single `for (int chain...)` aggregation path remain below
the native function in `st_cpg_omp_csr_annot.cpp`. Chain order and
`keep_chains` handling are unchanged; no aggregation was moved or duplicated.

## 8. Existing R conversion

All `cpg_raw_*`, `cpg_annot_chains_raw_v1()` and `cpg_annot_raw_v1()` binding
conversion remains inline in `st_cpg_omp_csr_annot.cpp` and unchanged. Field
names/order, R types, dimensions, dimnames, actual `NULL`, schema and optional
chain payload behavior are preserved.

## 9. Frozen references

Learned-annotation raw references: 3/3 exact. Learned-annotation formatted
references: 3/3 exact. Fixtures and tolerances were not changed or regenerated.

## 10. Reproducibility

Exact checks pass for repeated calls, one/two cores, reversed core execution
order through the permanent Phase 9A suite, intervening canonical group fits,
explicit chain seeds, retained/dropped chains, fixed coefficients,
inclusion-only, multiplier-only and joint learning configurations, explicit
update frequency/bounds and disabled LD-swap.

## 11. Unsupported cases

Current annotation/marker alignment, matrix and coefficient dimensions,
intercept representation, probability/multiplier bound validation,
equal-sample-size/shared-`ww` construction and unsupported trait-dimension
behavior remain unchanged. No support was broadened.

## 12. Protected backends

Hashes and focused/full tests confirm unchanged canonical BayesC, BayesR,
SBayesRC, fixed-prior BayesC and group BayesC, plus block-eigen, BED, scheduled
and multivariate implementations. Generated wrappers, `NAMESPACE`, public
signatures and raw schema remain unchanged.

## 13. Tests

- Phase 9A plus Phase 9F1 focused suite: 103 passed, 0 failed, 0 warned, 0 skipped.
- New Phase 9F1 tests: 33 passed, 0 failed, 0 warned, 0 skipped.
- Broad framework/annotation/schema/backend focused suite: 2,332 passed,
  0 failed, 0 warned, 0 skipped.
- Final full suite: 4,254 passed, 0 failed, 0 warned, 0 skipped.

## 14. Deviations and blockers

No implementation deviation or blocker occurred. The extraction compiled on
the first attempt. The first broad/full run exposed eight historical tests that
still froze the intentionally moved source by its Phase 9A hash; only that hash
was refreshed, after which the full suite passed. The header continues to
resolve surrounding function names lexically by design for Phase 9F1. Existing
Rcpp diagnostic output is part of the moved numerical block and remains unchanged.

## 15. Recommended Phase 9F2 task

> replace the lexically dependent learned-annotation execution include with an explicit typed execution context and callable core using the Phase 9A contracts, while retaining the existing R result conversion and wrapper-level aggregation until all frozen references pass again.

## 16. Readiness marker

PHASE 9F1 COMPLETE — LEARNED-ANNOTATION BAYESC EXECUTION BLOCK EXTRACTED
