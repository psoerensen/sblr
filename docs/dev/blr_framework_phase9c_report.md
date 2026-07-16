# Unified BLR Framework Phase 9C Report

## 1. Executive summary

Fixed-prior CSR BayesC is canonicalized and stabilized. Its public route uses one typed numerical core, one typed result, one binding-layer converter, and one wrapper-level multichain aggregation path. No numerical statement or public contract changed.

## 2. Repository baseline

The baseline was clean on branch `master` at starting commit `7766551` (`Migrate fixed-prior CSR BayesC to typed core`), which is also the Phase 9B3 commit. `git diff --check` was clean. R 4.4.1 UCRT on x86_64 Windows with Rtools44 compiled and loaded the package. The baseline full suite passed 4,052 tests with zero failures, warnings, or skips; the Phase 9B3 focused suite passed 2,916 tests, including Phase 9A/9B1/9B2/9B3 counts of 67/28/16/24.

The original Phase 9B3 benchmark used 2,000 markers, 100 retained iterations, 25 burn-in iterations, and five repetitions. Representative medians were 0.09--0.11 seconds for one chain/one core, 0.21 seconds for two chains/one core, and 0.19 seconds for two chains/two cores; completed-fit whole-process RSS was approximately 128.7--129.6 MiB. These timings are retained as historical measurements, not performance claims.

## 3. Phase 9B3 structure inventory

Phase 9B3 already had the binding-neutral `CsrPriorBayesCExecutionContext`, sole `run_csr_prior_bayesc()` core, typed `CsrPriorBayesCExecutionResult`, sole `stblr_csr_prior_bayesc_result_to_raw()` converter, sole wrapper aggregation loop, and three permanent raw/formatted fixtures.

Classification:

- Retain permanently: all four stable architecture names, the typed context/result contracts, the converter, wrapper aggregation, helper functions, and permanent Phase 9A fixtures.
- Remove now: active documentation saying the route was merely ready for canonicalization, and migration-comparison wording in the canonical benchmark.
- Rename for canonical use: the new benchmark and report use Phase 9C canonical-baseline terminology.
- Consolidate: permanent Phase 9C structural assertions protect singularity and binding neutrality; historical phase tests remain as regression provenance.
- Defer with justification: `CsrPriorBayesCLdFriendsView`/`LDLDFriendsPrior` and the raw payload member remain because they are active internal contracts used by preserved numerical helpers and conversion. Renaming would add churn without improving ownership. The fixture helper retains its Phase 9A name because that name records permanent fixture provenance.

No unused context field, unused result field, duplicate converter, duplicate aggregation alias, active selector, or fallback was found.

## 4. Files changed

- `src/blr_csr_prior_bayesc_core_impl.h`: clarified canonical implementation-detail and inclusion safety; numerical code is unchanged.
- `tests/testthat/test-blr-framework-phase9c.R`: added permanent architecture, reference, reproducibility, limitation, and protected-file checks.
- `tools/benchmarks/blr_phase9c_csr_prior_bayesc.R`: added the manual canonical timing/RSS baseline.
- `docs/dev/blr_framework_implementation_plan.md`: marked the route canonical and Phase 9C complete.
- `docs/dev/blr_model_capability_matrix.md`: marked fixed-prior CSR BayesC canonical while retaining adjacent-policy status.
- This report records implementation and validation.

## 5. Canonical execution path

The unchanged public fixed-prior CSR route performs R validation/alignment, native CSR and prepared-state construction, typed context construction, `run_csr_prior_bayesc()`, typed-result conversion, wrapper aggregation, `stblr_raw_v1` validation, and canonical raw-to-fit formatting.

## 6. Cleanup and naming

Stable architecture names were retained. Migration-only active documentation and benchmark wording were replaced by canonical terminology. No code alias was removed because the inventory found none that was both migration-only and unused. Historical phase names in tests and reports remain intentionally as provenance.

## 7. Numerical core and implementation header

Exactly one `run_csr_prior_bayesc()` implementation and one active fixed-prior marker loop remain in `blr_csr_prior_bayesc_core_impl.h`. The conventional include guard and inline definition prevent duplicate definitions, and the header is included only by `st_cpg_omp_csr_prior.cpp` after required native helpers. It neither redefines compiler/Armadillo configuration nor contains Rcpp/Python result construction. No numerical statement changed.

## 8. Native adapter

The adapter is limited to argument decoding, alignment, prior validation/preparation, CSR/native-state preparation, typed context construction, core invocation, result conversion, wrapper aggregation, and exception translation. It contains no MCMC/marker loop, RNG draw, effect/residual/variance update, global-parameter update, or posterior formula.

## 9. Result converter

`stblr_csr_prior_bayesc_result_to_raw()` remains the sole binding converter. It preserves field names/order/types, dimensions/dimnames, classes, actual `NULL`, schema metadata, prior metadata, global outputs, variance/VLE/VLD outputs, diagnostics/timing/failures, LD-swap and selection fields, chain payloads, and input metadata.

## 10. Wrapper-level multichain aggregation

Each core invocation returns one chain's typed execution result, which is converted once. The single wrapper loop preserves chain order and forms existing means and standard deviations while honoring `keep_chains`. It remains outside the core because aggregation is schema-aware and moving it could change chain order or raw layout. No duplicate or fallback aggregation exists.

## 11. Fixed-prior policy

`pi_marker` and `vb_multiplier` values, dimensions, marker order, and interactions with global `pi` and `B` are unchanged. `updatePi`, `updateB`, all formulas, and the public validation path are unchanged. The unsupported shared-LD multiple-trait construction remains rejected.

## 12. Ownership

CSR/operator storage, summary statistics, prepared priors, `pi_marker`, and `vb_multiplier` are borrowed immutable inputs. Effects, indicators, residuals, variances, global parameters, RNG/distributions, accumulators, diagnostics, and workspace remain chain-owned mutable state.

## 13. Permanent regression fixtures

The three unmodified Phase 9A configurations cover marker probability only, variance multiplier only, and both inputs, including single/multiple chains, core counts, explicit seeds, retained/dropped chains, supported update flags, disabled LD-swap, and selection fields. Fixture generation remains manual and is not run by tests.

## 14. Exact reference results

Phase 9C matched 3/3 raw references and 3/3 formatted references exactly, including values, storage, dimensions, names/classes, actual `NULL`, ordering, schema, prior metadata, variance/diagnostic, LD-swap, and selection fields.

## 15. Reproducibility

Permanent and inherited tests passed for repeated calls, one/two cores, reversed `1,2,2,1` core order, an intervening annotation-aware fit, explicit chain seeds, retained/dropped chains, all three fixed-prior policies, supported update flags, and disabled LD-swap. Only documented execution metadata is normalized.

## 16. Public API and schema

Public arguments/defaults, native public signatures, routing, `NAMESPACE`, generated wrappers, `stblr_raw_v1`, and formatted-fit fields are unchanged. Actual present-but-`NULL` fields and supported one-marker/one-trait dimensions remain intact.

## 17. Protected backends

Starting-commit hashes and focused/full tests confirm canonical BayesC, BayesR, SBayesRC, group BayesC, learned-annotation BayesC, block-eigen, BED, scheduled, and multivariate production files are unchanged. Group and learned-annotation remain production-unchanged, contract/reference-ready, and not migrated.

## 18. Performance and memory baseline

`Rscript tools/benchmarks/blr_phase9c_csr_prior_bayesc.R` ran a tiny workload and seven 2,000-marker primary configurations with one warm-up and five timed repetitions. Phase 9C medians were 0.46 (`pi_marker`), 0.47 (`vb_multiplier`), 0.49 (both, 1x1), 0.88 (both, 2x1), 0.97 (both, 2x2 retained), 0.48 (`updateB`), and 0.47 seconds (`updatePi`). IQRs were 0.01--0.05 seconds except the tiny workload. Completed-fit whole-process RSS was 127.6--130.0 MiB for moderate workloads and 134.8 MiB for the tiny first row.

The measurements use `system.time()` and `ps::ps_memory_info()` after each completed fit; RSS is not an interval peak. R did not report a compiler string, while compilation used Rtools44. A contemporaneous rerun of the unchanged Phase 9B3 script produced matching 1x1 medians of 0.41--0.48 seconds and RSS of 127.6--130.1 MiB, demonstrating that the earlier 0.09--0.21-second historical run is not directly comparable under the current Windows process conditions. Multichain timings were order-sensitive in the contemporaneous runs. Because Phase 9C changed only a comment in native code and identical scripts showed the same current regime, there is no evidence of a code-induced runtime or memory regression; no speed improvement is claimed.

## 19. Test results

Baseline: Phase 9A 67, Phase 9B1 28, Phase 9B2 16, Phase 9B3 24, focused suite 2,916, and full suite 4,052 passed. The final Phase 9C file passed 26/26, and the final full suite passed 4,078/4,078. Failures: zero. Warnings/skips: zero test warnings/skips; R reports only that `testthat` was built under R 4.4.3.

## 20. Deviations and blockers

No blocker remains. Interval peak RSS was unavailable, so completed-fit RSS is reported accurately. No optional larger-than-2,000-marker workload was run. Historical Windows timings were noisy and process-dependent, so the contemporaneous comparison is used and no speed claim is made. No native field or alias was removed because none was proven obsolete and unused.

## 21. Recommended next phase

> begin Phase 9D1 by mechanically extracting the group BayesC execution block while preserving group mapping, group order, probability and variance-multiplier updates, normalization, RNG ordering, public schema, and all Phase 9A group references.

## 22. Readiness marker

PHASE 9C COMPLETE — FIXED-PRIOR BAYESC CANONICALIZED AND STABILIZED
