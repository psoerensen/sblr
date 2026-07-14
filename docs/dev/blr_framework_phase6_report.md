# Unified BLR Framework: Phase 6 Report

## 1. Executive summary
Ordinary-CSR BayesR is canonical and stabilized with exact behavior preserved.
## 2. Repository baseline
`master` at `05805eb`, initially clean; R 4.4.1, Rtools44/GCC 13.2, C++17/OpenMP; baseline 3,620 tests and the Phase 5B compact benchmark passed.
## 3. Phase 5B structure inventory
Retained: typed context/result, one templated core/converter, and six permanent fixtures. Renamed: migration-era core/converter names. Removed: Phase 5-only comments and one unused context control. Deferred: the broader Phase 5A result vocabulary, which remains a validated contract.
## 4. Files changed
Canonical names/comments, structural tests, benchmark, plan/matrix, and this report changed; no backend outside BayesR changed.
## 5. Canonical call path
Public R validation -> unchanged native entry/operator construction -> borrowed typed context -> `run_csr_bayesr()` -> typed result -> one binding converter -> unchanged raw/fit.
## 6. Cleanup and renaming
`run_bayesr_execution` became `run_csr_bayesr`; the converter received the stable `stblr_csr_bayesr_result_to_raw` name; extraction-only wording and an unused field were removed.
## 7. Numerical core and implementation header
One guarded implementation header contains one marker loop. Both operator types instantiate it in the established Armadillo translation unit; no configuration macro is redefined.
## 8. Native adapter
It validates/prepares inputs, constructs the operator/context, invokes the core, and converts the result. No MCMC loop remains in the binding body.
## 9. Result conversion
One binding-layer converter preserves all raw fields, ordering, types, classes, `NULL`, chains, diagnostics, selection, and LD-swap payloads.
## 10. Ownership
Operator/prepared data, scales, and priors are borrowed immutable resources. Effects, components, residuals, variances, probabilities, RNG/distributions, accumulators, diagnostics, and workspace are chain-owned.
## 11. Permanent regression fixtures
Six configurations cover 1/2 chains, 1/2 cores, explicit seeds, multiple traits, retained/dropped chains, default/explicit scales, fixed/updated probabilities, summaries, disabled LD-swap, and stable selection fields.
## 12. Exact regression results
Raw 6/6 exact; formatted 6/6 exact.
## 13. Reproducibility
Repeated calls, 1/2 cores, 1-2-2-1 order, intervening BayesC, explicit seeds, multiple traits, and chain retention are exact.
## 14. Component behavior
Scales, null component, probabilities/traces/finals, counts, assignments, and summaries are exact.
## 15. Public API and schema
Arguments, defaults, signatures, routing, generated wrappers, `NAMESPACE`, `stblr_raw_v1`, and formatted fit are unchanged.
## 16. BayesC protection
Protected files are unchanged and permanent exact/reproducibility tests pass.
## 17. Block-eigen protection
Protected files and route are unchanged; it continues to instantiate the shared core and is not marked independently migrated.
## 18. Performance and memory baseline
Command: `Rscript tools/benchmarks/blr_phase6_csr_bayesr.R`. The moderate fixture has 2,000 markers, two traits, four components, 80 iterations/20 burn-in. Five-run medians were 1.72s (1/1), 3.25s (2/1; high order variability), 0.97s (2/2), and 0.79s (2/2 retained chains). Sampled whole-process RSS was 128.1--129.3 MiB. RSS is post-run process sampling, not instantaneous sampler-only peak. No material code-path regression is indicated and no speed claim is made.
## 19. Test results
Phase 1--5 and new Phase 6 exact/structural suites, CSR BayesR/component, BayesC, block-eigen, schema, consistency, and interface tests are included in the final result: 3,648 passed, zero failures, zero warnings, and zero skips.
## 20. Deviations and blockers
Process interval peak sampling was unavailable; sampled whole-process RSS and its limitation are recorded. No correctness blocker remains.
## 21. Recommended next phase
> begin a bounded contract and reference phase for CSR BayesRC/SBayesRC, preserving its current mathematics and annotation-specific behavior before attempting execution migration.
## 22. Readiness marker
PHASE 6 COMPLETE — CSR BAYESR CANONICALIZED AND STABILIZED
