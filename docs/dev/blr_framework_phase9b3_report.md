# Unified BLR Framework: Phase 9B3 Report

## 1. Executive summary

Fixed-prior CSR BayesC migration is complete with one typed numerical core and one named binding-layer result converter. Numerical behavior, wrapper-level multichain aggregation, public routing, raw schema, and formatted fits remain unchanged.

## 2. Repository baseline

Phase 9B3 started clean on branch `master` at commit `cde7330` (`Activate typed fixed-prior BayesC execution boundary`), the committed Phase 9B2 baseline. Phase 9B1 is `c7ff7cb`; Phase 9A is `e89c920`. R 4.4.1, Rtools44/GCC 13.2, C++17, and OpenMP were used. `Rcpp::compileAttributes()`, native compile/load, and the baseline full suite completed successfully: 4,028 passed, zero failures, warnings, or skips.

The Phase 9A directly comparable fixed-prior baseline used 2,000 markers and 100 iterations. Median seconds at one chain/one core, two chains/one core, and two chains/two cores retained were `0.11`, `0.68`, and `0.68`; completed-fit whole-process RSS ranged from 129.2 to 138.0 MiB. RSS was measured after fits and was not an interval peak.

## 3. Phase 9B2 structure inventory

| Item | Phase 9B3 classification | Resolution |
|---|---|---|
| `CsrPriorBayesCExecutionContext` | retain permanently | Explicit borrowed execution dependencies remain active. |
| `run_csr_prior_bayesc()` and its one marker loop | retain permanently | Numerical header unchanged. |
| `CsrPriorBayesCExecutionResult` | retain permanently | Used through execution, aggregation handoff, and conversion. |
| `cpg_raw_marker_matrix`, trace, trait, selection, and chain helpers | retain permanently | Private implementation details of the sole converter. |
| `cpg_prior_raw_v1()` field builder | centralize/rename now | Renamed to `stblr_csr_prior_bayesc_result_to_raw()` and made explicitly typed-result based. |
| Direct `.raw` return from the native single execution adapter | rename for stable use | The adapter now returns the typed result; wrapper aggregation borrows its raw payload. |
| Wrapper chain seed resolution, aggregation, summaries, and retention | retain permanently | Kept single and operation-for-operation unchanged. |
| Phase 9B1 inline-converter assertion | remove now | Replaced with permanent one-named-converter protection. |
| Group/learned source hashes | retain permanently | Adjacent unmigrated policies remain byte protected. |
| Generic raw vector storage inside the typed result | defer | Canonical typed field decomposition belongs to Phase 9C; changing it here could affect conversion or aggregation. |

No unused context or result field was proven removable. No old/new route selector, fallback, duplicate converter, or duplicate aggregation path existed.

## 4. Files changed

- `src/st_cpg_omp_csr_prior.cpp`: typed single-execution return, one named converter, and mechanical typed aggregate-to-converter handoff.
- Phase 9B1/9B2 tests: replaced obsolete inline-converter wording and assertions.
- `tests/testthat/test-blr-framework-phase9b3.R`: permanent architecture, exact-reference, reproducibility, trait-restriction, namespace, and adjacent-policy guards.
- `tools/benchmarks/blr_phase9b3_csr_prior_bayesc.R`: post-migration correctness, runtime, and completed-fit RSS benchmark.
- Implementation plan and capability matrix: migrated/ready-for-Phase-9C status.
- This report: migration inventory, validation, and baseline.

The numerical core/type headers, public R code, generated wrappers, schemas, and protected production backends were not modified.

## 5. Final execution path

The unchanged public fixed-prior CSR route performs R decoding, validation, and alignment; prepares fixed marker priors, CSR data, and native state; constructs `CsrPriorBayesCExecutionContext`; invokes `run_csr_prior_bayesc()`; preserves wrapper chain aggregation; passes one typed aggregate result to `stblr_csr_prior_bayesc_result_to_raw()`; and returns unchanged `stblr_raw_v1` for unchanged canonical formatting.

## 6. Centralized result converter

`stblr_csr_prior_bayesc_result_to_raw()` is the sole named fixed-prior binding converter. It consumes `CsrPriorBayesCExecutionResult` plus binding metadata and retains the existing private matrix/trace/trait/chain helpers. Field names and order, storage modes, dimensions, dimnames, marker/trait/chain ordering, classes, actual `NULL`, schema/version, fixed-prior metadata, global probability and variance outputs, traces, VLE/VLD, timing, diagnostics, LD-swap, selection, optional chains, and input metadata remain exact.

The converter performs no sampling, updates, aggregation, or numerical inference and the binding-neutral core constructs no R objects.

## 7. Wrapper-level multichain aggregation

Each `run_csr_prior_bayesc()` invocation produces one chain's native 24-slot typed payload. The wrapper resolves seeds in existing order, invokes each chain once, accumulates sums/squares/minima/maxima and LD-swap counts once, and optionally appends retained chain payloads once. Only the final aggregate typed result is converted. Chain order and `keep_chains` behavior are unchanged.

Moving aggregation into the core would alter the proven seam and risk raw slot shapes, chain ordering, and optional-chain schema; it remains intentionally outside.

## 8. Native adapter boundary

The adapter retains argument validation/alignment, fixed-prior and CSR preparation, context construction, typed core invocation, wrapper aggregation, one converter call, and exception propagation. It contains no active MCMC loop, marker traversal, RNG draw, probability/variance update, residual update, or posterior accumulation formula.

## 9. Numerical core

Exactly one `run_csr_prior_bayesc()` implementation and one active fixed-prior marker loop remain in `src/blr_csr_prior_bayesc_core_impl.h`. Its numerical statements, RNG construction/calls, arithmetic and branch order, OpenMP directives/scheduling, accumulation, and diagnostics were not changed in Phase 9B3.

## 10. Fixed-prior policy preservation

`pi_marker` remains borrowed immutable marker-by-trait input and replaces global `pi` only in marker draws. `vb_multiplier` remains borrowed immutable input and scales global marker variance. Marker order, trait order, global `pi` and `B` interactions, `updatePi`, `updateB`, validation, and fixed-policy metadata are unchanged. The documented multiple-trait construction with incompatible shared-LD scaling remains rejected.

## 11. Exact frozen references

All three fixed-prior raw references match exactly (3/3), and all three formatted references match exactly (3/3). Fixtures were not regenerated.

## 12. Reproducibility

Repeated calls, one/two cores, reversed `1,2,2,1` order, intervening annotation-aware execution, explicit chain seeds, retained/dropped chains, fixed probability input, fixed multiplier input, both inputs, supported update flags, and disabled LD-swap remain exact under the existing documented timing/path/core-metadata normalization only.

## 13. Public API and schema

Public arguments/defaults, native exported signature, routing, `NAMESPACE`, generated wrappers, `stblr_raw_v1`, actual `NULL` fields, marker/trait dimensions, and formatted fit are unchanged.

## 14. Protected backends

Git comparison against `cde7330` confirms no Phase 9B3 changes to canonical CSR BayesC, BayesR, or SBayesRC; group or learned-annotation BayesC; block-eigen; BED/BayesRC; scheduled CSR/BED; or multivariate implementations. Existing exact-reference, reproducibility, LD-swap, selection, annotation/alpha, schema, routing, and backend inventory suites protect these routes.

## 15. Performance and memory

Command: `Rscript tools/benchmarks/blr_phase9b3_csr_prior_bayesc.R`.

The benchmark warms each row once and runs five repetitions. Moderate workloads use 2,000 markers, one trait, 100 iterations plus 25 burn-in iterations. Median seconds were: probability-only `0.14`, multiplier-only `0.14`, both inputs 1×1 `0.14`, both 2×1 `0.27`, both 2×2 retained `0.27`, both with `updateB` `0.16`, and both with `updatePi` `0.14`. Primary both-input 1×1 times were `0.14,0.14,0.15,0.15,0.13` (range 0.13–0.15; IQR 0.01). Two-chain/one-core times were 0.26–0.30; two-chain/two-core retained times were 0.25–0.28.

Moderate completed-fit whole-process RSS ranged from about 128.1 to 128.9 MiB. Measurement used `system.time()` and `ps::ps_memory_info()` after completed fits; it is not interval sampling or peak RSS. Row order and a noisy tiny-workload outlier are recorded by the script.

Against Phase 9A's directly comparable medians `0.11/0.68/0.68` and RSS 129.2–138.0 MiB, the 1×1 difference is small in absolute terms and the multichain rows are lower. These timings are noisy and do not support a speed-improvement claim. No unexplained material runtime or completed-fit-memory regression is present.

## 16. Test results

- Baseline full suite: 4,028 passed, zero failures/warnings/skips.
- Phase 9A: 67 assertions.
- Phase 9B1: 28 assertions.
- Phase 9B2: 16 assertions.
- Phase 9B3: 24 assertions.
- Migration tests collectively: 135 assertions.
- Focused framework/fixed-prior/schema/canonical/protection suite: 2,916 passed, zero failures/warnings/skips.
- Final full suite: 4,052 passed, zero failures/warnings/skips.

## 17. Deviations and blockers

Memory sampling is completed-fit RSS rather than interval peak RSS and is labeled accordingly. The typed result intentionally retains its proven raw native payload until Phase 9C canonicalization. The trait restriction test protects the documented incompatible shared-LD construction rather than adding broader multiple-trait support. No blocker remains.

## 18. Recommended Phase 9C task

> canonicalize and stabilize fixed-prior CSR BayesC, remove remaining migration-only wording or aliases, retain permanent exact fixtures, establish the post-migration baseline as canonical, and then proceed to the group backend.

## 19. Readiness marker

PHASE 9B3 COMPLETE — FIXED-PRIOR BAYESC MIGRATED WITH BEHAVIOR PRESERVED
