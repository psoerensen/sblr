# BLR test-suite consolidation report

## 1. Executive summary

Phase 17F2 removes executable migration history that repeatedly ran identical fixtures and replaces it with permanent contract ownership. Canonical scalar, annotation, scheduled, packed-BED, and corrected multivariate references remain active once in the ordinary suite. Historical Phase 17B evidence remains immutable without executing the current sampler. Fresh-process and thread-environment coverage now has one extended owner.

The ordinary suite fell from 6,462 to 4,286 passing expectations (33.67%) and measured test time fell from 74.49 to 31.86 seconds (57.23%), with no production change or loss of a canonical reference family.

## 2. Baseline

Baseline commit was `69a83b6` on clean `master`. Compilation succeeded. The baseline contained 86 executed test files, 617 `test_that` blocks, 6,462 passes, 0 failures, 0 warnings, 13 skips, and 74.49 seconds of test-recorded elapsed time. The prior fast filter selected 22 files, 191 blocks, 1,737 passes, and 15.93 seconds. The workflow declared nine `SBLR_RUN_*` variables including peak RSS.

The slowest baseline files were individual BayesRC (6.30 s), selection-s (5.18 s), block eigen (4.22 s), annotation backends (3.61 s), CSR LD swap (3.34 s), CSR BayesR (2.16 s), CSR interface (2.10 s), backend field inventory (2.06 s), annotation interface (1.63 s), credible sets (1.61 s), raw schema (1.60 s), annotation chains (1.52 s), Phase 9G (1.24 s), Phase 9F3 (1.20 s), SBayesRC helpers (1.16 s), Phase 17C (1.11 s), BED interface (1.08 s), Phase 17F (1.02 s), Phase 3 (1.01 s), and Phase 9A (1.01 s).

## 3. Inventory method

`tools/audit/blr_test_suite_inventory.R` scans every test file and writes a machine-readable CSV covering blocks, expectation calls, fixture usage, process/thread calls, source assertions, hashes, schema/scientific checks, and protected-backend checks. Runtime and actual expectation counts were measured from `testthat` result records, aggregated by file. The static inventory guides review but does not replace semantic ownership decisions.

## 4. Contract ownership matrix

`docs/dev/blr_test_contract_ownership.md` is authoritative. Permanent numerical owners are Phase 3, 6, 8, 9C, 9E, 9G, 10D, 11D, 13E, 14E, and 17C. Public schemas remain owned by schema/interface suites. Phase 17E owns typed-core architecture; Phase 17F owns typed finalization; Phase 17B owns historical MT evidence only.

## 5. Tier model

- Tier 1 fast: canonical references, focused units, public/schema contracts, scientific identities, and essential architecture assertions.
- Tier 2 complete ordinary: all inexpensive supported and experimental integration tests, without enabled child-process or memory work.
- Tier 3 extended/manual: representative fresh/thread processes, peak RSS, experimental RNG risk, and benchmarks.

## 6. Phase 17B–17F consolidation

Phase 17B now reads and hashes historical fixtures and proves the historical fixed-B and denominator defects; it never executes current MT production. Phase 17C is the sole corrected raw/formatted, update-control, retained-count, normalization, scientific-identity, and same-process owner. Phase 17D retains one minimal assertion that the old lexical extraction has one typed successor. Phase 17E retains typed binding-neutral contracts, explicit ownership, one core/Gibbs loop, and fit-local RNG boundaries. Phase 17F retains the finalized type/finalizer, arithmetic separation, positional mapping, and future naming/operator documentation. Repeated Phase 17C executions and active-source hash matrices were removed from 17D–17F.

## 7. Scalar and packed-BED consolidation

Superseded migrations were removed while their final canonical owners remain: Phases 1–2 into Phase 3; 5A–5B into Phase 6; 7A–7B3 into Phase 8; 9B1–9B3 into 9C; 9D1–9D3 into 9E; 9F1–9F3 into 9G; 10A–10C3 into 10D; 11A–11C3 into 11D; 13A–13D into 13E; 14A–14D into 14E; and 15A into 15B. Their fixtures and historical reports were not deleted.

## 8. Source-assertion policy

Retained source assertions cover boundaries not efficiently observable at runtime: binding-neutral headers, singular public/core/finalizer calls, absence of adapter RNG/samplers, absence of worker identity, and positional compatibility. Exact incidental occurrence counts and repeated whole-source scans were removed from superseded migrations. Update controls, retention, normalization, covariance identities, and reproducibility are behavioral contracts.

## 9. Hash policy

Hashes remain for frozen RDS fixtures, historical binary evidence, and deliberate generated-wrapper immutability. Repeated whole-file hashes for actively maintained numerical sources were removed from Phase 17B–17F. Remaining legacy function-region hashes are temporary protection where active and legacy functions share a translation unit; the exit is disposition or separation of that legacy route.

## 10. Shared helpers

`helper-blr-test-contracts.R` supplies only repository-root text loading, frozen-fixture MD5 assertion, and extraction of the public MT function region. Existing `helper-source-architecture.R` remains the zero-safe source assertion owner. Helpers do not hide scientific comparisons.

## 11. Extended reproducibility

`test-blr-extended-reproducibility.R` is gated by `SBLR_RUN_EXTENDED_REPRODUCIBILITY=true` and covers one ordinary scalar CSR route, one scheduled CSR route, one packed-BED route, one corrected MT route, plus MT thread-environment equivalence. Model-specific seed/task unit tests remain with their models.

## 12. CI changes

The fast workflow selects the permanent canonical phase owners plus stable raw-schema, backend-consistency, CSR-interface, and BED-interface categories. This is a documented transition toward wholly subsystem-named files. Extended CI now uses only `SBLR_RUN_EXTENDED_REPRODUCIBILITY` and `SBLR_RUN_PEAK_RSS`; seven superseded phase-specific fresh variables plus the unused generic fresh variable were removed.

## 13. Defect sensitivity

`tools/audit/blr_phase17f2_mutation_sensitivity.R` applies controlled in-memory mutations and requires the owning contract predicate to reject each mutation. No working-tree file is written. Results were 8/8:

| Mutation | Permanent owner | Detected |
|---|---|---|
| remove `updateB` guard | Phase 17C behavior / Phase 17E boundary | yes |
| restore strict `it > nburn` | Phase 17C retained-count contract | yes |
| restore absolute thinning | Phase 17C retained-count contract | yes |
| wrong pi denominator | Phase 17C normalization / Phase 17F finalizer | yes |
| change logical-chain seed constant | Phase 15B seed unit contract | yes |
| change legacy output position | Phase 17F positional contract | yes |
| swap MT trait orientation | Phase 17C exact formatted reference | yes |
| introduce Rcpp into MT types | Phase 17E binding-neutral contract | yes |

## 14. Before/after metrics

| Metric | Before | After | Change |
|---|---:|---:|---:|
| executed test files | 86 | 52 | -34 (-39.53%) |
| `test_that` blocks | 617 | 363 | -254 (-41.17%) |
| passing expectations | 6,462 | 4,286 | -2,176 (-33.67%) |
| test-recorded elapsed | 74.49 s | 31.86 s | -42.63 s (-57.23%) |
| fast passes | 1,737 | 1,404 | -333 (-19.17%) |
| fast elapsed | 15.93 s | 14.02 s | -1.91 s (-11.99%) |
| ordinary skips | 13 | 2 | -11 |
| workflow `SBLR_RUN_*` variables | 9 | 2 | -7 |

Extended reproducibility adds five expectations and passed in 7.8 seconds locally when enabled.

## 15. Canonical reference preservation

Ordinary CSR BayesC, CSR BayesR, CSR SBayesRC, fixed-prior, group, learned-annotation, scheduled CSR, packed-BED BayesC, packed-BED BayesR, packed-BED BayesRC, and corrected dense MT reference owners all passed. Raw/formatted schema and public interface suites remain active.

## 16. Historical fixture preservation

All Phase 17B and Phase 17C RDS hashes are unchanged. Historical defect evidence is read directly and no current sampler execution is used as its expectation. No canonical fixture file changed.

## 17. Production protection

No numerical C++, public R, native signature, wrapper, schema, fixture, RNG, scheduling, or `NAMESPACE` file changed. No CSR, block-eigen, operator, or statistical implementation began.

## 18. Deviations and blockers

The fast filter still names final phase-owner files as a safe transition; permanent subsystem filenames should replace that regex only when those owners are deliberately renamed. Mutation checks use controlled in-memory transformations rather than recompiling eight deliberately defective package variants; each transformed contract was rejected and no mutation entered the working tree. `R CMD check` reached DESCRIPTION validation and stopped on the repository's pre-existing missing `Author` and `Maintainer` fields; compilation and the full test suite succeeded independently. No Phase 17F2 test-consolidation blocker remains.

## 19. Recommended next phase

Audit and formalize reuse of the canonical scalar CSR representation, marker-alignment contract, and validation vocabulary for a trait-specific multivariate LD-operator bundle, now using the consolidated permanent test architecture.

## 20. Readiness marker

PHASE 17F2 COMPLETE — BLR TEST CONTRACTS CONSOLIDATED
