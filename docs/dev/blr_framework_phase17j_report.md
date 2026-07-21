# Phase 17J: public trait-specific MT CSR interface

## 1. Executive summary

The public `mtblr_csr()` interface is active and routes validated R inputs
through the single Phase 17I CSR numerical path, named raw schema, validator,
and stable formatter.

## 2. Repository baseline

Work began clean on `master` at `18b2177`. Both workflows parsed. Hosted CI
status was not visible locally. R 4.4.1 and Rtools44 were used.

## 3. Phase 17I architecture verification

`mtblr_csr_raw_internal()` calls `mtblr_csr_internal()`, which owns all CSR
resource preparation and calls the one representation-neutral Gibbs core and
existing finalizer. No numerical core statement changed.

## 4. Previous public data architecture

`sblr()` remains the dense reference route. `make_summary_stats()` previously
returned numerical sufficient statistics and basic marker fields;
`make_sparse_ld()` recorded CSR construction fields without the canonical
provenance data frame now required at the public MT boundary.

## 5. Public API

The exported signature is the documented `mtblr_csr(stats, Glist,
ld_prefix, ld_metadata, trait_metadata, marker_policy, sample_overlap, method,
n, sets, b, h2, pi, models, pimodels, vg, vb, ve, ssb_prior, sse_prior,
updateB, updateE, updatePi, nub, nue, nit, nburn, nthin, seed, verbose)`.

## 6. Accepted stats forms

One multi-trait stats object or one named list of single-trait objects is
normalized to trait-major native vectors. Ambiguous structures fail.

## 7. Glist and LD resolution

One Glist, one per trait, one explicit prefix, or `nt` prefixes are accepted.
Explicit prefixes without Glist require metadata; simultaneous sources must
agree.

## 8. LD sharing modes

The modes are `fully_shared_operator`, `shared_correlation_reference`, and
`trait_specific_reference`. Differing `ww` with one raw reference expands the
prefix to `nt` before native scaling.

## 9. Marker-domain policy

Strict mode requires exact IDs and order. Reorder mode changes only R-side
stats and named initial effects; LD descriptors must already share one order.

## 10. Marker-intersection policy

Automatic intersection, union, insertion, and CSR subsetting are absent.
Mismatch is an error and `marker_intersection_policy` is recorded as `error`.

## 11. Allele-orientation policy

Explicit alleles match exactly. Swaps and complement-only matches stop. BED
inputs can be accepted by construction when stats and LD provenance and marker
domain agree. External unidentified orientation stops.

## 12. Scaling

Only `standardized_genotype` cross-products are accepted.

## 13. Trait metadata

Output columns are trait ID, study ID, ancestry, population, LD reference,
sample size, resolved prefix, and sharing mode.

## 14. Sample-overlap semantics

Only `not_modeled` is accepted. The likelihood consumes marginal `yy` and does
not reinterpret residual covariance as sampling-error overlap.

## 15. Residual covariance

The public policy requires diagonal positive-definite `ve` and `sse_prior`.

## 16. Models, priors, and sets

Binary unique patterns include the null state and preserve ordering; model
probabilities are normalized. Marker sets are a complete disjoint partition.
Existing dense prior formulas and full B/G prior covariance support are kept.

## 17. Named MT raw schema

Version 1 namespaces are `schema`, `meta`, `marker`, `trace`, `variance`, `pi`,
`model`, `diagnostics`, `data`, and `alignment`, with dimensions documented in
`mtblr_raw_schema.md`.

## 18. Native raw converter

The converter only names and shapes the finalized legacy-compatible values. It
contains no RNG, updates, posterior division, or biological alignment.

## 19. Formatted fit

`mtblr_fit` contains the established marker, trace, covariance, probability,
and correlation fields plus marker order, model patterns, metadata, alignment,
input contract, and raw schema version. No CSR payload is returned.

## 20. Numerical reduction

Shared, trait-scaled repeated-prefix, and trait-specific/independent-pattern
public configurations agree with the internal Phase 17I route at `1e-12`.

## 21. Dense-route protection

The Phase 17C raw and formatted fixture owners remain unchanged.

## 22. Internal-route protection

The Phase 17I reductions remain the owner of dense/CSR numerical equivalence.

## 23. Scalar protection

The permanent Phase 3, 6, 8, 9C, 9E, 9G, and 10D scalar owners remain active.

## 24. Research-route protection

The public code contains no call to `mtblr_eigen()` or
`mtblr_cpg_omp_csr()`.

## 25. Metadata and output size

Only compact descriptors and reports are returned; CSR row pointers, indices,
values, and diagonal buffers are not retained in the fit.

## 26. Generated wrappers and namespace

Rcpp generation adds one internal named-raw registration. `NAMESPACE` adds
only `export(mtblr_csr)`.

## 27. Mutation sensitivity

The audit covers marker IDs/order, allele swaps, scale, off-diagonal `yy`,
overlap, raw validation, research routing, duplicate Gibbs loops, pattern names,
and duplicate raw wrappers.

## 28. Tests and CI

The Phase 17J owner passed 37 expectations. The fast permanent tier passed.
The complete ordinary suite passed 4,417 expectations in 55 files and 380
blocks with zero failures/warnings, two opt-in skips, and 152 seconds wall time.

## 29. Package check

Compilation, loading, and documentation generation succeeded. `R CMD check`
reached package tests but failed on pre-existing installed-check-context source
path/hash assumptions, compiler-level fixture differences, and existing Rd
warnings (157 failures, 75 warnings, two skips, 3,536 passes in that check
context). The source-tree ordinary suite is green.

## 30. Deviations and blockers

No numerical blocker is known. Explicit allele metadata is supported only when
provided; the repository has no established universal Glist allele-label field.

## 31. Recommended next phase

Audit and formalize the canonical scalar block-eigen storage, filtering,
ownership, marker-block, and operator contracts for trait-specific MT reuse
before implementing a multivariate block-eigen route.

## 32. Readiness marker

PHASE 17J COMPLETE — PUBLIC TRAIT-SPECIFIC MT CSR INTERFACE ACTIVE
