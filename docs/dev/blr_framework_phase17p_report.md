# Unified BLR Framework Phase 17P report

## 1. Executive summary

The supported public `mtblr_bed()` joint individual-level multivariate BayesC
adapter is active over the unchanged Phase 17O numerical route.

## 2. Repository baseline

Work began on clean `master` at
`7f22a237f32c1c1a016f70f26ad3ad26e37e4461`. The baseline used R 4.4.1 UCRT
and Rtools44 GCC/G++/Fortran 13.2; local BLAS identity was not reported by
`extSoftVersion()`. Hosted CI was not visible locally. Baseline source tests
reported 4803 passes and two established skips. Baseline built checking was
environmentally blocked while resolving dependency repositories and also
reported the established long fixture-path note.

## 3. Phase 17O verification

The unchanged numerical path is BED validation, one `PackedBedMatrix`, one
`BedPackedGenotypeView`, marker maps/`X'Y`/order, sample residual, one fit-local
RNG and serial MT BED core, `MtDefaultCoreResult`, shared finalizer, legacy
adapter, and `mtblr_legacy_to_raw()`.

## 4. BED preparation reuse

The adapter calls `.make_bed_marker_data()` once. That existing path remains
the owner of Glist/BED/file/column/sample/ID/order/frequency/name/default-set
alignment.

## 5. Shared MT helper reuse

Models use `.mtblr_models()`, sets `.mtblr_sets()`, covariances `.mtblr_cov()`,
raw validation `.validate_mtblr_raw()`, and formatting `.as_mtblr_fit()`.

## 6. Public API

The exact signature is recorded in `blr_mt_bed_public_contract.md`; it exposes
public data, covariance, initialization, MCMC, memory-warning and metadata
controls but no paths, native indices, residuals, workers, or chains.

## 7. Glist contract

Exactly one BED-backed Glist supplies one genotype matrix for every trait.
Lists of per-trait Glists and incomplete BED metadata fail closed.

## 8. Phenotype forms

Numeric vectors, matrices, and all-numeric data frames are accepted. Names are
preserved after alignment or generated as `T1`, `T2`, and so on.

## 9. Individual alignment

Established explicit-row and phenotype-ID behavior, unmatched-ID warning,
duplicate rejection, selected order, multiple files, default markers, and
explicit columns are preserved and described in `fit$alignment`.

## 10. Complete-data validation

Aligned phenotypes require `n>1`, `nt>0`, finite complete values, unique trait
names, and positive finite column variances. No rows are dropped for phenotype
missingness.

## 11. Centering

`center=TRUE` subtracts aligned column means in R. `center=FALSE` verifies the
Phase 17O scale-aware tolerance without transforming data. Units and variances
are retained and centering evidence is returned.

## 12. Covariates

Non-NULL `covar` is rejected with a pre-adjustment instruction. No design,
projection, residualization, or fixed-effect computation occurs.

## 13. Genotype scale and frequencies

Only standardized `scale=TRUE` is accepted. Selected Glist frequencies must be
finite and strictly inside `(0,1)` and are never recomputed or clipped.

## 14. Marker metadata

Metadata contains marker ID, chromosome/file, BED column, and allele frequency
in final order, plus explicit alleles only when exact Glist alignment exists.

## 15. Trait metadata

Default trait IDs or exactly ordered user rows are enriched with sample size,
pre/post means, variance, centering, covariance mode, and individual data level.

## 16. Models and probabilities

Shared MT pattern expansion/restrictive modes, validation, null pattern,
probability normalization, names, and size limit are unchanged.

## 17. Sets

Explicit sets use the shared disjoint complete partition contract. Defaults are
converted from stable BED preparation labels and do not reorder markers.

## 18. Covariance modes

Full residual covariance is the public default and accepts full SPD E/prior.
Diagonal mode is an explicit reduction and requires exact diagonal inputs.

## 19. Prior defaults

Defaults use centered sample variances, `h2`, non-null probability, marker
count, and the documented Phase 17P `vg`, `ve`, `vb`, `ssb_prior`, and
`sse_prior` formulas. Degrees of freedom exceed `max(2, nt-1)`.

## 20. Initialization

Marker-by-trait matrices and trait lists normalize to latent/effective/state
matrices. Zero defaults, b-only inference, explicit inactive latent effects,
model-pattern membership, binary masks, and `b=state*beta` are enforced.

## 21. Memory estimate

Component byte formulas and total GiB are returned. The configurable warning
states explicitly that this is analytical working memory, not measured RSS or
peak RSS, and never stops execution by itself.

## 22. Native call

The adapter calls `mtblr_bed_internal()` exactly once with aligned descriptors,
centered Y, trait-major initialization/prior data, native sets/models, method 4,
and the selected residual-covariance mode.

## 23. Raw enrichment

Named raw version 1 gains public model/pi names, data metadata, and alignment;
it does not duplicate Y, packed genotype bytes, sample residuals, or genetic
values.

## 24. Fit formatting

The shared MT formatter is called once. Phase 17P then adds BED diagnostics,
phenotype preprocessing, and the memory estimate without a second formatter.

## 25. Public/internal equality

Portable tests own full/diagonal fixed and all-update, restrictive/full,
one/three-trait, multiple-set, initialization, row/ID/file/column, and centering
equality at exact or `1e-12` tolerance.

## 26. Alignment tests

Portable tests cover vectors, matrices, data frames, all/ID/explicit rows,
multiple files, default/explicit columns and sets, generated names, metadata,
and failure-closed Glist/ID/row/frequency/marker contracts.

## 27. Centering tests

Tests verify aligned means are subtracted, units/variance retained,
precentered input accepted, uncentered opt-out rejected, and native Y centered.

## 28. Covariance tests

Tests cover full and diagonal defaults/inputs, diagonal rejection of
off-diagonal matrices, E preservation when disabled, and trait-dimension degree
of freedom errors.

## 29. Initialization tests

Tests cover zero, b-only, state-only, explicit matrices/lists, inactive latent
values, and dimension/finiteness/binary/pattern/masking failures.

## 30. Memory-warning tests

Exact small-formula tests cover below/above threshold, `Inf`, invalid controls,
and the required non-peak-RSS wording.

## 31. Raw-schema tests

Both internal and public calls retain version 1, backend `mt_bed_bayesc`, data
level `individual`, standard namespaces, and validated MT BED diagnostics.

## 32. Phase 17O protection

Phase 17O production numerical files and generated native registrations are
unchanged; its permanent marker, execution, reduction, initialization,
reproducibility, ownership and diagnostic owners remain authoritative.

## 33. Scalar BED protection

Existing BayesC, BayesR, BayesRC and marker-interface behavior remains owned by
its permanent tests and is not routed through MT BED.

## 34. Summary-MT protection

Dense, CSR, and block-eigen routes retain their formals, numerical backends,
named raw outputs and shared general formatter semantics.

## 35. Public API protection

Only `mtblr_bed()` is newly exported. Existing public formals/defaults remain
unchanged and no `ncores` or `nchains` control is introduced.

## 36. Mutation sensitivity

The Phase 17P audit owns all required alignment, singularity, centering,
missingness, scale, initialization, covariance, memory, output, schema,
protected-route, forbidden-feature and export mutations.

## 37. Benchmark

The maintenance benchmark reports small/moderate full/diagonal public and
prepared internal total times, analytical memory, object bytes, and unavailable
separated/RSS quantities without performance claims.

## 38. Installed-check behavior

All execution and scientific/public contracts are portable. Only source-text
singularity and forbidden-route assertions narrow-skip outside a source tree.

## 39. Generated documentation and namespace

Roxygen generates one `mtblr_bed.Rd` page and one `export(mtblr_bed)` line.
No Rcpp wrapper or native registration is generated or changed.

## 40. Tests and CI

The fast framework filter includes `17p`. Phase 17O/17P focused tests and the
exact fast filter passed. The final source suite reported 4,917 passes, zero
failures, zero warnings, and two established opt-in skips. Built-package tests
completed successfully.

## 41. Package check

The built-tarball check completed with zero errors and zero warnings. It
reported three classified notes: the established long Phase 17C fixture path,
the established compiled `std::cout` finding, and installed size 5.1 MB
(`libs` 4.2 MB). The size note is a package-size threshold signal; Phase 17P
made no native or DLL change.

## 42. Diff hygiene

Final checks cover fixtures, generated files, whitespace/EOL, build/check
artifacts, and the absence of commits or pushes.

## 43. Deviations and blockers

The first baseline package check was blocked by transient dependency-repository
access. The final authorized check succeeded. One sandboxed post-check source
rerun could not create a `processx` child pipe; the same fresh-process test
passed in the earlier source run, built check, and final unsandboxed source run.
No production or numerical workaround was made.

## 44. Recommended next phase

> audit and formalize the multichain, OpenMP dispatch, chain-level state ownership, deterministic chain seeds, aggregation, diagnostics, and optional retained-chain output contracts for public individual-level MT BayesC before implementing parallel or multichain execution.

## 45. Readiness marker

PHASE 17P COMPLETE — PUBLIC INDIVIDUAL-LEVEL MT PACKED-BED INTERFACE ACTIVE
