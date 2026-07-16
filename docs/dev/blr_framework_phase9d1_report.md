# Unified BLR Framework Phase 9D1 Report

## 1. Executive summary

The group CSR BayesC execution block was mechanically extracted into one guarded implementation header without changing numerical behavior, stochastic trajectory, public behavior, inline R conversion, or wrapper-level multichain aggregation.

## 2. Repository baseline

Phase 9D1 started clean on branch `master` at commit `8691427` (`Canonicalize fixed-prior CSR BayesC`), the committed Phase 9C baseline. `git diff --check` passed. R 4.4.1 UCRT on x86_64 Windows used Rtools44/GCC, C++17, and OpenMP. `Rcpp::compileAttributes()` left generated wrappers unchanged, native compile/load succeeded, and the focused baseline passed 1,914 assertions with zero failures, warnings, or skips.

## 3. Extraction seam

The moved block began immediately after group mapping/validation, group-order and prior preparation, CSR and LD-swap-friend construction, and deterministic marker ordering. Its first statement was allocation of `bm_mat`; its final statement was `return result;` from `stblr_cpg_omp_csr_group_annot_single()`. The 424 source lines formerly at lines 922--1345 of `st_cpg_omp_csr_group.cpp` were moved verbatim. The source now defines the intended-translation-unit guard, includes `blr_csr_group_bayesc_core_impl.h`, and undefines the guard at that same lexical seam.

## 4. Files changed

- `src/st_cpg_omp_csr_group.cpp`: replaced the 424-line execution block with the guarded include.
- `src/blr_csr_group_bayesc_core_impl.h`: contains the verbatim execution block plus a conventional include guard, intended-translation-unit check, and implementation-detail comment.
- `tests/testthat/test-blr-framework-phase9a.R`: replaced only the intentionally obsolete group production-source MD5.
- `tests/testthat/test-blr-framework-phase9d1.R`: adds singularity, inclusion-site, policy-access, exact-reference, reproducibility, coverage, and protected-file checks.
- Historical Phase 7A/7B2/7B3/8/9B3/9C protection tests: replaced only the intentionally obsolete group production-source MD5; all adjacent hashes remain unchanged.
- This report records the extraction and validation.

## 5. Extracted execution content

The header contains trait-task allocation and OpenMP execution, chain-local RNG construction, effect/residual and group-state initialization, group updates, the MCMC and marker loops, group lookup and policy access, optional LD swaps, global and residual variance updates, group probability and multiplier updates, normalization invocation, posterior accumulation, group summaries, VLE/VLD calculations, diagnostics, native result allocation, and native result population.

## 6. Group-policy preservation

Marker-to-group mapping, group count/order, zero-based execution indices, group-specific inclusion probabilities and variance multipliers, Beta priors, update flags and timing, `normalize_group_vb`, global `pi`/`B` interactions, marker traversal, arithmetic and branch order, residual/variance updates, RNG engine/distribution construction and invocation order, seed mapping, OpenMP static scheduling, retained-sample timing, accumulation order, LD-swap behavior, selection fields, and diagnostics are unchanged. No group was recoded, reordered, converted to annotations, or removed.

## 7. Wrapper-level aggregation boundary

Wrapper-level multichain aggregation was already outside the contiguous single-chain execution function and remains unchanged in `stblr_cpg_omp_csr_group_annot()`. It invokes the extracted single-chain block once per chain, preserves chain order and `keep_chains`, aggregates once, and constructs chain summaries once. It was not moved because it is a separate established boundary rather than part of the contiguous numerical block.

## 8. Existing R conversion

The `cpg_raw_*` helpers and `cpg_group_raw_v1()` remain inline in `st_cpg_omp_csr_group.cpp`, after the extracted single-chain function, and are unchanged. Public signatures, generated wrappers, `NAMESPACE`, raw schema, classes, field names/order/types, and formatted output remain unchanged.

## 9. Frozen references

All permanent Phase 9A group references matched exactly: raw 3/3 and formatted 3/3. Comparisons include values, storage, dimensions/dimnames, names/classes, actual `NULL`, marker/trait/group/chain order, zero-based mapping, schema, group probabilities/multipliers/counts, normalization metadata, variance and VLE/VLD traces, diagnostics, LD-swap and selection fields, and optional chain payloads. No fixture was regenerated and no tolerance was weakened.

## 10. Reproducibility

Repeated identical calls, one versus two cores, reversed `1,2,2,1` core order, an intervening learned-annotation fit, explicit chain seeds, retained and dropped chains, the nontrivial two-group marker mapping, updated group variance multipliers, normalization enabled and disabled, and disabled LD swap remained exact. Existing focused tests cover fixed and updated group probability/multiplier policies, group chains, group LD swap, group variance components, and intervening canonical backend execution. Only declared execution metadata (`input$ncores`, timing, and fixture LD prefix) was normalized where already documented.

## 11. Dimension and unsupported cases

Current production restrictions are unchanged: marker groups use zero-based native indices, every declared group must be represented, initial group vectors and priors must match `ngroup`, shared-LD multiple-trait inputs require equal sample sizes and matching `ww`, and no broader trait construction was introduced. Both supported normalization modes remain available; no unsupported case was forced.

## 12. Protected backends

MD5 and focused/full tests confirm canonical BayesC, BayesR, SBayesRC, canonical fixed-prior BayesC, learned-annotation BayesC, block-eigen, BED, scheduled, and multivariate production sources are unchanged. Generated wrappers, `NAMESPACE`, public R/native signatures, and schema artifacts are unchanged. The only protected production hash intentionally replaced is the extracted group source; the new tests protect its single include, sole active marker loop, inline converter, and sole aggregation path structurally.

## 13. Tests

- Phase 9A: 67 assertions passed.
- Phase 9D1: 33 assertions passed.
- Final focused framework/group/canonical/schema/interface suite: 2,764 passed.
- Final full suite: 4,111 passed.
- Failures: zero.
- Test warnings: zero.
- Skips: zero.

Native compile/load succeeded before and after extraction. R only reported that the installed `testthat` package was built under R 4.4.3; existing compiler unused-function warnings were unchanged and were not test warnings.

## 14. Deviations and blockers

The first combined Phase 9A/9D1 run found one new test-only helper-name error; `phase9a_fixture()` was corrected to the existing `phase9a_inputs()` helper without changing production code or fixtures. The first broad focused run then identified six historical MD5 assertions that encoded the intentionally obsolete group-source layout. Only those group-source hash entries were updated; all other hashes remain protected, and the rerun passed. No numerical, reference, compiler, API, schema, or unsupported-case blocker remains.

## 15. Recommended Phase 9D2 task

> replace the lexically dependent group execution include with an explicit typed group execution context and callable core using the Phase 9A contracts, while retaining the existing R result conversion and wrapper-level aggregation until all frozen references pass again.

## 16. Readiness marker

PHASE 9D1 COMPLETE — GROUP BAYESC EXECUTION BLOCK EXTRACTED
