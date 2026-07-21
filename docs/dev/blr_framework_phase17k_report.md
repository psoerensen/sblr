# Phase 17K report

## 1. Executive summary

One binding-neutral block-filtered storage/view contract is active for all three retained scalar block-eigen routes. Numerical construction and sampler policies are unchanged; no MT sampler or public API was added.

## 2. Repository baseline

Baseline: clean `master` at `7e8b5170e4697d33e667a5b94cdfa522ef44dbe6`. R 4.4.1 UCRT, Rtools44 GCC/GFortran 13.2.0; BLAS identity was not reported by `extSoftVersion()`. Workflows parsed. Existing block-eigen tests passed 193/193. Phase 17J2 source/installed test architecture was retained.
The committed Phase 17J2 ordinary suite contained 55 files, 387 blocks, 4,437 passing expectations, zero failures/warnings, and two opt-in skips. Its built check had zero errors/warnings and the same two legacy notes retained below.

## 3. Previous architecture

`EigenBlock` and `BlockEigenOperator` lived in `st_ld_operator.h`; filter types/diagnostics lived separately; three translation units duplicated parsing and diagnostics conversion. The operator already stored reconstructed float-packed dense blocks.

## 4. Complete implementation inventory

| Route | Model/status | Builder/storage | Filter / wy | Numerical owner |
|---|---|---|---|---|
| scalar block eigen | BayesC, internal experimental | canonical builder/storage/view | all three; hard projects | BayesC operator template |
| scalar block eigen | BayesR, internal experimental | same | all three; hard projects | BayesR operator template |
| scalar block eigen | SBayesRC, internal experimental | same | all three; hard projects | SBayesRC operator template |
| scalar CSR counterparts | supported public | SparseLdCsrStorage/View | none | canonical CSR cores |
| `mtblr_eigen()` | unsupported research | XX row lists, unrelated | unrelated | legacy independent loop |

## 5. BED and standardization contract

SNP-major, supplied file/column/row order; codes 00/01/10/11 map 2/missing/1/0; missing maps to standardized zero; finite `p` in `(0,1)`; float Z, double accumulation/reconstruction, float runtime packing, double wy.

## 6. Block-domain contract

Zero-based starts begin at zero, strictly ascend below `m`, imply contiguous nonempty blocks, and cover `[0,m)`. Native code inserts no chromosome/file boundary.

## 7. Filter contract

Hard uses `max(tau,0.01)` and largest-eigenvalue fallback. Fixed ridge uses clamped `eta/(1+eta)`. Ledoit–Wolf retains its existing data-derived formula and clamp. Inactive parameters remain ignored.

## 8. Hard-truncate transformation

The existing reconstructed matrix and `wy D^-1/2 V_k V_k' D^1/2` projection are unchanged and protected by an independent float-stage oracle.

## 9. Ridge transformations

Both ridge modes reconstruct `(1-a)A+a diag(diag(A))`; neither transforms `wy`.

## 10. Canonical storage

`sblr::core::BlockEigenStorage` owns marker count, block floats, mappings, and runtime diagonal. `EigenBlock`/`BlockEigenOperator` are compatibility aliases, not duplicate definitions.

## 11. Canonical view

`sblr::core::BlockEigenView` borrows immutable storage and implements identical Armadillo/vector apply/rebuild traversal.

## 12. Central validation

Checks marker/block presence, contiguity/coverage, overflow-safe packed lengths, finiteness, mapping identity, positive diagonal, and exact packed/runtime diagonal identity.

## 13. Builder ownership

BED and quadratic workspaces are transient. One completed owner outlives all task borrowers. No per-chain copy, MCMC-time BED I/O, or eigendecomposition occurs.

## 14. Scalar active integration

BayesC, BayesR, and SBayesRC factories call `build_block_eigen()`, move one owner into their operator context, and all operator methods delegate through `view()`.

The preserved audit result is:

```text
BLOCK_EIGEN_OPERATION_ORDER_PRESERVED=TRUE
BLOCK_EIGEN_RNG_ORDER_PRESERVED=TRUE
BLOCK_EIGEN_RETENTION_PRESERVED=TRUE
BLOCK_EIGEN_WY_TRANSFORM_PRESERVED=TRUE
BLOCK_EIGEN_FLOAT_PACKING_PRESERVED=TRUE
```

## 15. Diagonal audit

Input `ww`, unfiltered BED diagonal, and filtered float-derived runtime diagonal are distinct. Same-BED fixture differences are reported by tests/audit; equality of filtered diagonal and input ww is not required.

## 16. Marker-ranking audit

Ranking remains before builder consumption and uses current `wy` plus original input `ww`. Sampling/rebuild/LE use the operator diagonal. No switch occurred.

## 17. Diagnostics

Names and meanings remain `start,size,n_kept,mu_min,shrink`; one shared Rcpp converter now serves all models.

## 18. Numerical oracle

Portable multi-block tests cover a one-marker block, missing genotype, two wy rows, nonzero effects, hard floor/equality/subset/fallback, fixed ridge zero/positive/large eta, and Ledoit–Wolf.
All filter reconstructions, transformed/non-transformed `wy`, packed diagonals, diagnostics, and Armadillo/vector rebuilds passed (`1e-12` where float-rounded inputs coincide and `1e-10` for the independent eigensystem). The installed compiler required only the documented `1e-12` Armadillo/vector continuous comparison; mapping/rank/fallback results remained exact.

## 19. Scalar numerical protection

Existing BayesC/BayesR/SBayesRC block-eigen tests and permanent scalar CSR owners remain the model-level references; no fixture was regenerated.
All three modes passed for all three scalar model families. Existing tests retain multiple-trait/shared-operator and chain/seed coverage, `keep_chains`, fixed and sampled `selection_s` where implemented, model-specific B/E/pi or annotation/component updates, and unchanged raw metadata/diagnostics. LD-swap remains rejected; scheduled and public block-eigen execution remain unsupported.

## 20. Future MT compatibility

The view supports fully shared storage, shared boundaries with independent values/projections, and fully independent operators. Numerical sharing requires complete construction-input or verified content identity.

## 21. External-summary limitation

The route builds reference cross-products from BED but cannot prove external-study/reference sample or allele identity. Hard-projected external summaries require a future explicit provenance contract.

## 22. Complexity and memory

Packed storage and `O(sum s_b^2)` update/rebuild formulas are reported by the audit. Hard truncation remains dense-block runtime storage.

## 23. Benchmark

The analytical small/moderate benchmark reports storage, largest transient estimate, and sweep/rebuild visits; it is a regression signal, not peak RSS or a speed claim.
For 100 markers in four blocks of 25, retained bytes were 6,800 and a sweep visited 2,500 within-block values. For 2,000 markers in twenty blocks of 100, retained bytes were 436,000 and a sweep visited 200,000 values. Largest transient estimates were 40,000 and 1,120,000 bytes respectively under the script's stated workspace model.

## 24. Retained MT research comparison

The legacy route differs in input, storage, construction, transformation, diagonal, conditional, G/E, RNG, retention, result, and scientific status. It remains noncanonical, unmodified, internal, and unusable as an oracle.

## 25. Public-interface protection

Public ST routes remain CSR; no block-eigen selector/export was added. LD-swap rejection is retained.

## 26. Installed-check behavior

Portable inspection/oracle/validation blocks run installed. Only source-text architecture/disposition blocks use the Phase 17J2 narrow source-root skip.

## 27. Generated wrappers

One internal registered inspection wrapper was added; no namespace export or public signature changed.

## 28. Tests and CI

Focused Phase 17K passed 69/69. The focused legacy-hash/block-eigen suite passed 457/457 after actively maintained translation units were transferred from obsolete whole-file hash ownership to Phase 17K behavioral ownership. The final source suite passed 4,506 expectations across 56 files and 395 blocks, with zero failures/warnings and only two established opt-in skips (101.8 seconds test-recorded, 106.7 seconds wall). Fast CI includes `17k`; its owned selection passed through the full source run.

## 29. Package check

The fresh built-tarball check exited zero: zero ERRORs, zero WARNINGs, and two unchanged Phase 17J2 NOTES (the long frozen Phase 17C helper path and legacy `std::cout` use). Installed tests passed, including the portable Phase 17K filter/validation/oracle blocks; 97 skips comprise the established installed source-only inventory plus opt-in coverage. No new NOTE was introduced.

## 30. Mutation sensitivity

All fourteen required contract mutations were detected.

## 31. Diff hygiene

Only the single internal inspection wrapper/registration changed generated Rcpp files. `NAMESPACE` and fixtures are unchanged. `git diff --check` passes; final cleanup removed compiled objects/DLLs and no tarball or `.Rcheck` directory is in the repository.

## 32. Deviations and blockers

The inspection seam returns runtime matrices for small maintenance fixtures only; it is registered internally because installed checks cannot compile ad-hoc C++ harnesses. No public exposure results.

## 33. Recommended next phase

> implement an internal trait-specific MT block-eigen data view and execution route that reuses the Phase 17I shared MT BayesC Gibbs core, uses one canonical BlockEigenView and matching transformed `wy` per trait/study, and retains the corrected dense and CSR routes as numerical oracles without yet adding a public block-eigen interface.

## 34. Readiness marker

PHASE 17K COMPLETE — CANONICAL BLOCK-EIGEN CONTRACT FORMALIZED FOR TRAIT-SPECIFIC MT REUSE
