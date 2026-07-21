# Unified BLR Framework Phase 17J2 report

## 1. Executive summary

Phase 17J2 separates portable installed-package validation from source-checkout architecture audits. Scientific, public, schema, and fixture contracts remain active in both contexts; repository assertions skip narrowly when their assets are unavailable.

## 2. Repository baseline

Baseline: `master` at `944ce603122e3e4f9fe09d4e08b5beb9929fda41`, initially clean. R 4.4.1 UCRT uses GCC/GFortran 13.2.0. Both workflows parsed. Hosted CI was not locally visible.

## 3. Initial check reproduction

The initial built tarball check reproduced 3,536 passes, 157 failures, 75 test warnings, and 2 skips: status 1 ERROR, 2 WARNINGs, 6 NOTEs.

## 4. Failure inventory

Failures comprised repository-root/source reads, fixture paths coupled to that root, unqualified namespace internals, toolchain-specific runtime hashes/exact continuous comparisons, cross-test parsing, and documentation defects. The audit script reports occurrences by category.

## 5. Root-cause analysis

One `blr_test_root` incorrectly represented both repository and portable assets. `load_all()` also exposed internals that installed namespaces do not place in the test environment. GCC builds differed from frozen MSVC/Rtools representations only at roughly 1e-16 in continuous values; discrete trajectories did not change.

## 6. Test-context architecture

The final helpers define `blr_test_dir`, `blr_fixture_path()`, structurally validated optional `blr_source_root`, `blr_repo_path()`, and `blr_source_text()`.

## 7. Fixture migration

Fixture scripts and RDS/MD5 reads resolve beneath `testthat::test_path("fixtures", ...)`. Frozen fixture contents were not regenerated.

## 8. Source-only architecture tests

The complete block inventory and reasons are maintained in `blr_source_only_test_inventory.md`.

## 9. Installed scientific tests

Permanent scalar CSR, scheduled, packed-BED, dense MT, internal MT CSR, public MT CSR, raw-schema, metadata, alignment, backend consistency, and public-interface owners remain installed-required.

## 10. Namespace access

Production-R sourcing fallbacks were removed. Tests use the installed namespace for internals.

## 11. Shared helper refactor

MT CSR data builders now live in `helper-mtblr-csr-fixtures.R`; Phase 17J no longer parses Phase 17I or evaluates definitions in `.GlobalEnv`.

## 12. Hash and numerical-reference policy

Frozen files retain exact MD5. Active source hashes are source-tier. Same-build calls remain exact. Cross-toolchain continuous fixture comparisons use `1e-12`, except the learned-annotation `dm_sd` zero-boundary difference (`1.756119e-9`) uses a field-owner tolerance of `1e-8`; structures and discrete values remain exact.

## 13. Documentation fixes

The credible-set Rd percent operator was escaped correctly and `finemap_stblr_csr(verbose)` is documented.

## 14. Metadata fixes

R >= 3.5.0 is explicit for version-3 RDS fixtures; unused `RcppArmadillo` Imports was removed while LinkingTo remains; missing stats/utils imports were declared.

## 15. Source-tree validation

The final source-tree run covered 55 test files and 387 `test_that()` blocks, with zero failures and zero test warnings. Source architecture blocks executed; only the two established opt-in extended tests skipped. The final compiler-enabled validation took 599.8 seconds wall time (including compilation and verbose native output).

## 16. Built-package validation

The final fresh-tarball check ran 3,627 passing expectations with zero failures and zero test warnings; 94 blocks skipped (93 explicitly inventoried source-architecture blocks and the extended-reproducibility opt-in block). `R CMD check` exited zero with zero ERRORs and zero WARNINGs. Two NOTES remain: the frozen Phase 17C fixture helper has a long portable path, and legacy numerical sources use `std::cout`. Removing either would require renaming frozen evidence or editing protected numerical sources, so both are documented non-failing legacy limitations.

## 17. CI workflow

Fast CI includes Phase 17J2 and invokes the shared built-tarball check driver, which fails on ERROR or WARNING.

## 18. Numerical protection

Phase 17C, 17I, 17J, and scalar permanent owners remain unchanged numerically and retain their existing tolerances.

## 19. Public API protection

No Phase 17J signature, default, schema, marker/allele policy, or output changed.

## 20. Diff hygiene

Only test infrastructure, tests, documentation/metadata, workflow, and audit/check tooling changed. Numerical sources and fixtures are unchanged.

## 21. Deviations and blockers

The check status contains two NOTES rather than the literal text `Status: OK`: one long frozen-fixture path and one pre-existing compiled-code `std::cout` portability note. There are no ERRORs, WARNINGs, test failures, or documentation warnings. Phase 17I continues to define its local reusable builders as well as the new helper copy, but no test parses another test file or mutates `.GlobalEnv`; removing those local definitions would be maintenance-only follow-up.

## 22. Recommended next phase

> audit and formalize the canonical scalar block-eigen storage, filtering, ownership, marker-block, and operator contracts for trait-specific MT reuse before implementing a multivariate block-eigen route.

## 23. Readiness marker

PHASE 17J2 COMPLETE — R CMD CHECK RELIABLE ACROSS SOURCE AND INSTALLED CONTEXTS
