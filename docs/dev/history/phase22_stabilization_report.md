# BLR stabilization report

## 1. Executive summary

Phase 22 converts the development-history-oriented tree into permanent
scientific and software ownership without changing numerical behavior.

## 2. Starting repository state

Starting commit: `03db04a8eb6044b78a44f49c7366bf30cc03cd7f`.

## 3. Cleanup authority and compatibility policy

Backward compatibility was explicitly unnecessary; canonical scientific and
numerical contracts were protected.

## 4. Baseline validation

Baseline: 83 test files, 567 blocks, 6,133 expectations, zero failures/errors/
warnings, two opt-in skips; package check had zero errors/warnings and three
notes.

## 5. Cleanup inventory

See `../blr_cleanup_manifest.md`.

## 6. Neutrality harness

`tools/validation/blr_phase22_neutrality_matrix.R` serializes normalized fits
for 21 supported route combinations against a detached local baseline. Only
timing, process/path, and explicitly retired metadata are normalized.

## 7. Permanent test ownership

See `../blr_test_ownership.md`.

## 8. Removed phase-numbered tests

Thirty-eight phase-numbered test files were removed after their enduring
expectations were transferred. The active suite has 44 scientifically named
test owners and no phase-numbered test file.

## 9. Removed migration/compatibility tests

The dedicated historical source-hash owner and migration/alias assertions
embedded in phase tests were removed. Canonical invalid inputs remain tested.

## 10. Consolidated helpers

Sixteen helpers were reduced to four focused helpers for fixtures, operators,
convergence, and repository/source checks. `codetools::checkUsagePackage()`
reports no findings.

## 11. Fixture consolidation

Fifty-nine fixture artifacts were reduced to seven scientifically named files
(six deterministic reference builders plus their ownership README). No
scientific fixture values were regenerated.

## 12. Source-only test consolidation

Source-only checks are limited to permanent architecture contracts.

## 13. Audit consolidation

Four permanent audits replace 43 phase/history-oriented scripts:
architecture, generated interfaces, documentation, and mutation sensitivity.

## 14. Mutation consolidation

Twenty permanent mutations protect the current contracts. Redundant
phase-transition mutations were removed; the retained set covers API/model
semantics, schemas, RNG/trace ownership, preparation, generated interfaces,
memory guards, annotations, selected markers, convergence, and CI checking.

## 15. CI simplification

Fast and extended CI use permanent owners; benchmarks are manual.

## 16. Benchmark consolidation

Model, operator, diagnostics, and peak-RSS benchmarks replace 62 one-off
scripts. They remain manual or extended-tier tools.

## 17. Developer-documentation consolidation

Eleven active developer documents replace 120 historical files. Architecture,
models/capabilities, output, convergence, backends, testing, development, and
history each have one owner.

## 18. Phase-report disposition

Detailed phase reports were removed; Git history and `history.md` retain the
milestones. This completion record is the sole detailed stabilization record.

## 19. User-documentation consolidation

Nine maintained Quarto sources passed link/navigation inspection and rendered
successfully to a task-specific temporary output tree. Generated HTML/site
libraries and obsolete drafts were removed from the source tree.

## 20. Workflow/example consolidation

Eighteen workflow files were reduced to ten (README, one shared helper, and
eight canonical workflows). Each uses a small reproducible fixture or states
its external-data requirements.

## 21. Removed R compatibility code

Legacy model spellings and annotation aliases, the retired architecture
argument, and migration-only resolver code were deleted rather than deprecated.

## 22. R-layer simplification

Duplicate model resolution was removed and canonical resolvers remain narrow.
The installed usage audit reports no unused-variable or undefined-global
findings.

## 23. Final public fit schema

See `../blr_output_schema.md`.

## 24. Final raw schema

Raw schema version 1 and model semantics version 2 remain authoritative.

## 25. Route inventory and dispositions

See `../blr_backend_inventory.md`.

## 26. Removed native routes

Twelve dead native registrations/reference bridges were removed (55 to 43).
No canonical sampler or independently owned scientific oracle was removed.

## 27. Native-code simplification

Dead audit-only headers and phase binding code were deleted. Live kernel loop,
update, RNG, and accumulation order were intentionally retained.

## 28. Generated-interface cleanup

`Rcpp::compileAttributes()` reports 43 wrappers and 43 registrations; the
permanent generated-interface audit passes all 43 owners. There are 31 public
exports and no stale or duplicate registration.

## 29. Native-output cleanup

All active `std::cout`, `std::cerr`, `printf`, `fprintf`, and R/Rcpp console
output was removed from production native paths. Worker code remains R-API free.

## 30. Package-check note cleanup

The long fixture-path and legacy native-output notes were eliminated. The sole
remaining note is installed size (5.5 MB, including a 4.4 MB native library),
classified as the expected cost of the retained scientific kernels.

## 31. Numerical-neutrality evidence

All 21 normalized route results are bitwise identical (`identical()`) to the
detached Phase 21 installation. Posterior summaries, final states,
probabilities, covariance and variance quantities all match exactly.

## 32. RNG-neutrality evidence

Resolved logical seeds, stochastic states, and outputs are bitwise identical
for single/multichain and serial/OpenMP representatives.

## 33. Trace-neutrality evidence

Core/extended convergence bundles and explicitly selected marker traces are
bitwise identical to baseline. A compiler-layout-induced last-bit mismatch was
detected during cleanup and eliminated by retaining harmless local structure;
no numerical tolerance was accepted for the final same-route comparison.

## 34. Test-suite before/after

Before: 83 files, 567 blocks, 6,133 expectations, two opt-in skips, 173.16 s.
After: 44 files, 309 blocks, 3,687 passing expectations, one opt-in skip,
113.3 s; zero failures, errors, or test warnings in both validated endpoints.

## 35. Repository-size before/after

Tracked active files fell from 655 (11,770,073 bytes) to 278 (4,760,656
bytes). The temporary installed package fell from 6,543,358 to 5,716,884
bytes; its DLL fell from 4,723,712 to 4,543,488 bytes.

## 36. Fast CI

The exact permanent filter passed 977/977 expectations with no skips or
warnings in 42.5 s.

## 37. Full source suite

The final full suite passed 3,687/3,687 expectations with one justified
fresh-process opt-in skip and no failures, errors, or warnings in 113.3 s.

## 38. Installed tests

Built-package tests passed.

## 39. Package check

`R CMD check` completed with 0 errors, 0 warnings, and the one classified
installed-size note.

## 40. Architecture audit

All 19 permanent architecture guards pass. The generated-interface audit also
passes 43/43 registrations and permanent mutation sensitivity passes 20/20.

## 41. Documentation audit

The documentation audit passes: 17 active documents and eight candidate links
were checked; YAML, Rd, and Quarto validation also pass.

## 42. Diff and artifact hygiene

`git diff --check` passes after generated Rd normalization. Generated compiler,
DLL, tarball, check, Quarto, baseline-clone, and serialized comparison artifacts
were removed after validation.

## 43. Remaining historical/phase references

Only this completion record, `history.md`, the cleanup manifest's historical
paths, and the neutrality driver's baseline name retain milestone references.
No active test, audit, CI filter, public documentation, or workflow is
phase-organized.

## 44. Deviations and blockers

No scientific, RNG, convergence, or trace blockers remain. The installed-size
note is retained and classified; reducing it further would require unsafe or
scientifically unjustified native refactoring.

## 45. Long-term stable contracts

The canonical APIs, model semantics, schemas, RNG topology, operator
mathematics, and convergence engine are frozen by permanent owners.

## 46. Recommended next work

> prepare a release candidate, expand simulation and empirical validation of the MT BayesRC/SBayesRC methodology, and manage future changes through focused issues, releases, and named scientific capabilities rather than additional framework phases.

## 47. Readiness marker

PHASE 22 COMPLETE — BLR FRAMEWORK STABILIZED AND HISTORICAL SCAFFOLDING REMOVED
