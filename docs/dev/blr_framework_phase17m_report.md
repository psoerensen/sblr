# Phase 17M report

## 1. Executive summary

The public same-BED multivariate block-eigen route is active in the working
tree, subject to the final validation recorded below.

## 2. Repository baseline

Baseline commit: `e30ed7e42b4712eb302ec3df4e6ea82933bd7948`; branch `master`;
initial tree clean. R 4.4.1 uses Rtools44 GCC 13.2.0 and R's `Rblas`/`Rlapack`.
Hosted CI was not visible from this offline workspace. The first baseline
compile/test attempt lacked Rtools on `PATH`, and its package installation
failed with an MSYS signal-pipe error. Consequently no valid pre-change source
or installed total was recorded; this is an explicit baseline deviation.

## 3. Phase 17L verification

`mtblr_block_eigen_internal()` → descriptor parsing → BED read →
`build_block_eigen()` → owner/view construction → transformed `wy` →
`run_mt_block_eigen_core()` → shared default core → finalizer → legacy adapter.

## 4. Phase 17J reuse

Statistics normalization, marker metadata, scale normalization, alignment,
models, sets, covariance preparation, raw validation, and fit formatting are
reused through the existing Phase 17J helpers.

## 5. Public API

The exact signature is recorded in
`docs/dev/blr_mt_block_eigen_public_contract.md`.

## 6. Support boundary

Only same-selected-BED, by-construction statistics are supported.

## 7. Stats normalization

Multi-trait objects and named lists of single-trait objects retain normalized
per-trait genotype provenance.

## 8. Genotype provenance

Normalized paths, `n_bed`, rows/order, columns/file order, marker IDs/order,
allele frequencies, marker metadata, source, standardized scale, and analysis
sample size are compared exactly.

## 9. Glist resolution

One Glist is replicated as a candidate common reference; `nt` Glists remain
trait-specific. Conflicting `sparseLD` carriers are rejected.

## 10. Marker alignment

Strict and reorder-only policies reuse `.mtblr_align()`; intersections fail.

## 11. Allele orientation

Explicit alleles match exactly, or BED file/column provenance establishes
same-Glist construction.

## 12. Block normalization

Mandatory one-based public starts are validated and converted to zero-based
native starts.

## 13. Filter normalization

Hard, fixed-ridge, and Ledoit-Wolf filters accept scalar or per-trait
configuration.

## 14. Operator sharing

Auto, shared, trait-specific, and trait-specific-shared-boundary metadata are
implemented without native-owner deduplication.

## 15. Adapter-result refactor

One adapter execution owns construction, diagnostics, transformed `wy`, views,
core execution, finalization, and legacy conversion.

## 16. Shared raw conversion

CSR and block-eigen named raw use the single `mtblr_legacy_to_raw()` converter.

## 17. Named raw

The block backend is `mt_block_eigen_bayesc`, data level `summary`, schema
version 1.

## 18. Block diagnostics

Diagnostics come from the actual build and include complete native block
fields, owner metadata, and trait-owner mapping.

## 19. Transformed wy

Hard traits report projected `wy`; ridge traits report unchanged `wy`.

## 20. Summary ww policy

`ww` is validated by construction and is not passed as the runtime diagonal.

## 21. Fit formatting

`.as_mtblr_fit()` remains the sole general formatter; block diagnostics and
operator metadata are appended.

## 22. Public/internal reductions

Focused shared, ridge, mixed-filter, differing-row-count, one-trait raw, and
trait-specific reductions pass in the Phase 17M owner.

## 23. Provenance rejection tests

External source, missing construction fields, rows, frequencies, blocks,
filters, overlap, and reference-policy failures are covered.

## 24. Scientific identities

Marker summaries, residuals, states, traces, covariances, probabilities, and
retention continue to originate in the shared Phase 17L result.

## 25. CSR protection

Phase 17J and 17J2 tests protect public/raw CSR behavior.

## 26. Internal-route protection

Phase 17L reductions protect legacy numerical output.

## 27. Scalar protection

Phase 17K remains the scalar canonical block-eigen owner.

## 28. Research-route protection

`mtblr_eigen()` remains unsupported and is not called.

## 29. Installed-check behavior

Execution, provenance, validation, raw, fit, diagnostic, and reduction tests
are portable. Only source architecture assertions may skip without a source
tree.

## 30. Generated wrappers and namespace

One internal raw registration and one public export are generated.

## 31. Mutation sensitivity

All 17 required public and architecture mutation probes are detected.

## 32. Benchmark

The small regression benchmark reports shared-hard, shared-ridge, and
trait-specific mixed-filter timing, owner count, operator size, fit size, and
separable build/internal estimates.

## 33. Tests and CI

Focused and exact fast tiers pass. The full source suite passes 4,605
assertions with zero failures and zero warnings; only the two established
opt-in skips remain. Installed tests pass during `R CMD check`.

## 34. Package check

Built-tarball check has zero errors and zero warnings. The two unchanged notes
are the long Phase 17C fixture path and pre-existing compiled `std::cout`
usage.

## 35. Diff hygiene

Generated changes are one raw wrapper/registration and one public export. No
fixture regeneration, tarball, `.Rcheck` directory, or compiled artifact is
retained.

## 36. Deviations and blockers

The only deviation is the unavailable valid pre-change baseline total
described above. Final source and installed validation is complete.

## 37. Recommended next phase

audit and formalize the existing packed-BED owner/view, standardized marker decoding, phenotype/covariate preparation, sample-space residual, full or diagonal residual-covariance, missing-phenotype, ownership, memory, and output contracts required for an internal individual-level multivariate BayesC route before implementing that sampler.

## 38. Readiness marker

PHASE 17M COMPLETE — PUBLIC SAME-BED MT BLOCK-EIGEN INTERFACE ACTIVE
