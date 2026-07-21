# Public MT-BLR CSR contract

## 1. Purpose

`mtblr_csr()` is the canonical scalable public multivariate BayesC summary-
statistics interface.

## 2. Statistical target

It uses the corrected Phase 17C model and the single Phase 17I dense/CSR Gibbs
implementation, finalizer, and result contract.

## 3. Accepted stats forms

Inputs are either one multi-trait object (`wy`/`ww` lists, `yy`, `n`, marker
and trait IDs) or a named list of single-trait objects. Marginal beta, SE, and
z-score inputs are not converted.

## 4. LD resources

One Glist, one Glist per trait, one explicit prefix, or one prefix per trait is
accepted. Prefix-only use requires explicit LD metadata.

## 5. Sharing modes

Identical diagonals plus one prefix use a fully shared operator. One raw
correlation prefix with differing `ww` is expanded per trait for trait-specific
scaling. Multiple prefixes use trait-specific references.

## 6. Marker domain

`strict` requires identical order. `reorder_stats` reorders only R-side stats
to the first LD descriptor order. All LD descriptors must already agree.

## 7. Marker intersection

Missing or extra markers are errors. There is no automatic intersection,
union, zero insertion, or CSR subsetting.

## 8. Allele orientation

Explicit effect/other alleles must match exactly. Swaps and strand complements
are rejected; no automatic flip occurs. Missing labels are accepted only for
same-construction BED statistics and LD with matching marker provenance.

## 9. Scaling

Only standardized-genotype cross-products are supported.

## 10. Trait metadata

Trait, study, ancestry, population, LD reference, sample size, prefix, and
sharing-mode metadata accompany the fit. Unknown scientific metadata remains
missing rather than being inferred.

## 11. Sample overlap

Only `sample_overlap = "not_modeled"` is accepted. Marginal `yy` is used and
residual covariance is not interpreted as GWAS sampling-error correlation.

## 12. Residual covariance

The public policy is diagonal `ve` and diagonal `sse_prior`.

## 13. Models and sets

Patterns are unique binary vectors including the null pattern; probabilities
are normalized without reordering. Sets are a disjoint complete 1-based marker
partition.

## 14. Reproducibility

The supplied integer seed is passed unchanged to the fit-local native engine.

## 15. Raw schema

Native output is named `mtblr_raw` version 1 and validated before formatting.

## 16. Formatted fit

`mtblr_fit` preserves the established numerical fields and adds marker order,
patterns, marker/trait metadata, alignment, input metadata, and schema version.

## 17. Limitations

No automatic harmonization, block eigen, multichain, OpenMP MT execution,
marker union, or overlap likelihood is implemented.
