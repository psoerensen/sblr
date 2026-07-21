# Phase 17L report

## 1. Executive summary

Internal trait-specific MT block-eigen execution is active through the single Phase 17I Gibbs implementation.

## 2. Repository baseline

Baseline: `763319b0910cba02089730bf2182fae4355022fb`, branch `master`, initially clean. R 4.4.1 UCRT with Rtools44 GCC/GFortran 13.2.0 was used. Hosted CI state was not queried; both workflow files were locally parsed.

## 3. Phase 17K verification

`build_block_eigen()` produces the sole `BlockEigenStorage`; `BlockEigenView` borrows packed blocks, mappings, and its packed-float-derived diagonal.

## 4. Existing MT seam

Representation access remains trait count, marker count, diagonal, residual rebuild, and marker-difference application.

## 5. MT block-eigen types

`MtBlockEigenBundleView` contains `marker_count` and `vector<BlockEigenView> trait_operator`. `MtBlockEigenDataView` borrows `wy`, `yy`, and `n` and contains the bundle.

## 6. Central validation

Validation requires positive trait/marker counts, consistent summary and operator counts, positive `n`, valid canonical views, model row width, and in-domain set indices.

## 7. Representation access

Each trait uses its own view. Rebuild delegates to that view. Difference application subtracts the runtime diagonal once, then traverses the block off-diagonal in canonical local order.

## 8. Shared-core integration

`mtblr_block_eigen_internal()` → `run_mt_block_eigen_core()` → `run_mt_bayesc_core_impl()` → `MtDefaultCoreResult` → shared finalizer → shared legacy adapter.

## 9. Native internal signature

The signature is `(wy, yy, b, operator_descriptors, sets, B, E, ssb_prior, sse_prior, models, pi, nub, nue, updateB, updateE, updatePi, n, nit, nburn, nthin, seed, method=4)`.

## 10. Descriptor parsing

Descriptors follow the exact schema in the internal contract. Parsing and all basic validation precede construction; `method != 4` is rejected.

## 11. Shared-owner mode

One descriptor yields one packed-BED build, one owner, all rows transformed together, and `nt` views.

## 12. Trait-specific mode

`nt` descriptors yield `nt` independent owners, transformations, and views. No descriptor deduplication occurs.

## 13. Mixed filters

Hard-truncated rows are projected; ridge rows are unchanged, independently by trait.

## 14. Transformed-wy contract

The transformed container is used for ranking, the core, covariance updates, and the legacy output.

## 15. Runtime diagonal

The packed operator diagonal is used; the route deliberately accepts no competing `ww`.

## 16. Dense numerical reductions

Permanent tests cover shared hard truncation with fixed/all updates, fixed ridge, Ledoit–Wolf ridge, nonzero initialization, multiple sets, three traits, trait-specific values, independent boundaries, mixed filters, differing `n`, and one trait.

## 17. CSR numerical reductions

One-marker diagonal blocks and a nontrivial CSR-representable block are compared with the canonical CSR route.

## 18. Scientific identities

Tests protect transformed residual identities and legacy equality, which jointly protect finite/discrete states, retained counts, probabilities, covariance summaries, marker ordering, and trace dimensions.

## 19. Statistical-order audit

`PHASE17I_MT_STATISTICAL_STATEMENTS` and `PHASE17L_MT_STATISTICAL_STATEMENTS` are identical because the template body is unchanged. `STATISTICAL_ORDER_EQUIVALENT=TRUE`, `RNG_CALL_ORDER_EQUIVALENT=TRUE`, `UPDATE_CALL_ORDER_EQUIVALENT=TRUE`, `RETENTION_POLICY_EQUIVALENT=TRUE`, and `MARKER_CONDITIONAL_EQUIVALENT=TRUE`. Construction precedes RNG; transformed `wy` is used for ranking, core, and output.

## 20. Ownership, copies, and memory

Shared mode owns one operator; trait-specific mode owns `nt`. Both hold `nt` views, zero per-chain operator bytes, and perform zero MCMC-time BED reads/eigendecompositions. Packed bytes are `4 sum s(s+1)/2` per owner, mappings `8m`, diagonal `8m`, transformed data `8ntm`.

## 21. Performance

The small shared and mixed trait-specific cases took 0.02 and below 0.01 seconds for build-plus-MCMC; dense execution rounded below 0.01 seconds. The moderate 2-trait, 500-marker, 40-reference-sample case with twenty 25-marker blocks took 0.08 seconds for each route, with 26,000 packed bytes versus a 4,000,000-byte dense estimate. Timings are regression signals, not peak-memory claims.

## 22. External-summary limitation

The numerical route cannot prove marker/allele/reference provenance or validity of projecting external summary statistics.

## 23. Sample-overlap semantics

Sample overlap and cross-study sampling-error likelihoods remain unsupported.

## 24. Dense-route protection

The dense wrapper and shared core body are unchanged.

## 25. CSR-route protection

CSR types, adapter, and public route are unchanged.

## 26. Scalar block-eigen protection

The Phase 17K builder, view, filters, inspection, and scalar routes are unchanged.

## 27. Research-route protection

`mtblr_eigen()` remains unchanged, unsupported, and unused.

## 28. Public-interface protection

No public function, argument, selector, schema, or export changed.

## 29. Generated wrappers

One internal Rcpp wrapper and registration were generated; `NAMESPACE` is unchanged.

## 30. Installed-check behavior

Numerical and validation tests are portable. Only source-text architecture assertions use the Phase 17J2 narrow source-tree skip contract.

## 31. Mutation sensitivity

The audit covers wrong trait view, lost/projected `wy`, wrong output, omitted diagonal, duplicate loops/helpers, late construction, descriptor count, forced boundaries, copies, research/public routing, wrapper drift, and a competing `ww`.

## 32. Tests and CI

Phase 17L focused tests passed 44 expectations. The exact fast filter passed with its one established peak-RSS opt-in skip. The complete 58-file, 403-block source suite passed approximately 4,550 expectations with zero failures/warnings and the two established opt-in skips. The final source rerun also completed all tests successfully; its post-run summary calculation failed after testing because testthat result columns are list-valued.

## 33. Package check

The built-tarball driver completed with 0 errors and 0 warnings. Installed tests reported `OK`. Two unchanged notes remain: the existing long Phase 17C fixture path and existing `std::cout` references in legacy scalar translation units.

## 34. Diff hygiene

`git diff --check` passes. Fixtures and `NAMESPACE` are unchanged. Attribute generation added only the internal wrapper and one registration. No tracked binary, tarball, `.Rcheck`, or temporary artifact was created; Git's normal Windows LF-to-CRLF checkout notices do not indicate content-only EOL churn.

## 35. Deviations and blockers

The check driver removes its temporary logs after success, so it reports installed tests as `OK` rather than retaining a testthat expectation total. This does not affect numerical or package-check acceptance.

## 36. Recommended next phase

> design and expose a canonical public `mtblr_block_eigen()` interface that reuses Phase 17J statistics normalization, marker/allele alignment, trait metadata, overlap policy, named `mtblr_raw` version 1, and `mtblr_fit` formatting, while adding explicit BED-reference descriptors, selected rows/columns, allele frequencies, block boundaries, filter parameters, operator-sharing policy, transformed-summary provenance, and block diagnostics, using the validated Phase 17L internal route without changing its numerical implementation.

## 37. Readiness marker

PHASE 17L COMPLETE — INTERNAL TRAIT-SPECIFIC MT BLOCK-EIGEN CORE ACTIVE
