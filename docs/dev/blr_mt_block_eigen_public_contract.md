# Public multivariate block-eigen contract

## 1. Purpose

`mtblr_block_eigen()` is the supported public same-BED multivariate
block-filtered BayesC route.

## 2. Public function

```r
mtblr_block_eigen(
  stats, Glist, block_start,
  operator_sharing = c("auto", "shared", "trait_specific"),
  eigen_filter = "hard_truncate", eigen_tau = 0.01, eigen_eta = 0,
  summary_reference = "same_bed_by_construction", trait_metadata = NULL,
  marker_policy = c("strict", "reorder_stats"),
  sample_overlap = "not_modeled", method = "bayesC", n = NULL,
  sets = NULL, b = NULL, h2 = 0.5, pi = 0.001, models = NULL,
  pimodels = NULL, vg = NULL, vb = NULL, ve = NULL,
  ssb_prior = NULL, sse_prior = NULL, updateB = TRUE, updateE = TRUE,
  updatePi = TRUE, nub = 4, nue = 4, nit = 1000, nburn = 500,
  nthin = 1, seed = 1, verbose = FALSE
)
```

## 3. Statistical target

The target is the corrected Phase 17C/17I serial multivariate BayesC model.

## 4. Internal execution

Phase 17L is the sole numerical route. The public adapter calls
`mtblr_block_eigen_raw_internal()`, which builds owners once, transforms `wy`,
runs `run_mt_block_eigen_core()` once, finalizes once, and converts once.

## 5. Supported provenance

Statistics and operator must have the same normalized BED paths, original BED
sample count, selected rows and order, selected marker columns and file order,
allele frequencies, marker identities and orientation, and standardized scale.

## 6. Unsupported provenance

External GWAS statistics plus a separate reference panel are rejected even
when marker metadata agrees.

## 7. Statistics forms

Both the Phase 17J multi-trait statistics object and a named list of
single-trait statistics objects are accepted.

## 8. Glist forms

Supply one Glist or one Glist per trait. `Glist` is mandatory.

## 9. Marker alignment

`strict` requires identical order. `reorder_stats` permits R-side reordering
only when marker sets and all reordered BED provenance agree. Intersection is
an error.

## 10. Allele orientation

Explicit effect/other alleles must match exactly. Otherwise package BED
construction requires identical marker, file, and BED-column provenance.
Swaps and strand-complement-only matches are rejected.

## 11. Blocks

Blocks are mandatory. Public starts are one-based, strictly ascending, begin
at one, and are converted to zero-based starts only in native descriptors.
No boundary is inserted automatically; explicit blocks may cross file or
chromosome boundaries.

## 12. Filters

`hard_truncate`, `ridge_fixed`, and `ridge_lw` are supported. Filter, `tau`,
and `eta` may be scalar or per trait. The canonical hard floor is 0.01.

## 13. `wy` transformation

Hard truncation projects `wy`. Fixed and Ledoit-Wolf ridge leave `wy`
unchanged. The returned `wy` is the exact value consumed by the core.

## 14. `ww` policy

Summary `ww` remains finite and positive and is validated by construction, but
is not the runtime diagonal. The runtime diagonal belongs to the canonical
filtered block operator.

## 15. Operator sharing

`auto` shares only unambiguously common public inputs and provenance.
`shared` enforces that condition. `trait_specific` creates exactly one owner
per trait without deduplication.

## 16. Trait-specific references

Trait-specific owners may use independent Glist references, selected rows,
blocks, numerical values, and filters while retaining a common marker domain.

## 17. Sample size

Each analysis size must equal its selected reference-row count. `rows = NULL`
means all `n_bed` rows.

## 18. Sample overlap

Sample overlap is not modeled.

## 19. Residual covariance

Residual covariance and its prior are diagonal.

## 20. Models and sets

The Phase 17J model-pattern, probability, complete disjoint marker-set, prior,
initialization, and MCMC contracts are reused unchanged.

## 21. Named raw

The native boundary returns `mtblr_raw` version 1 with backend
`mt_block_eigen_bayesc`.

## 22. Formatted fit

The shared `.as_mtblr_fit()` formatter owns general fields. The public adapter
adds `block_diagnostics` and `operator_metadata`.

## 23. Diagnostics

Diagnostics record owner count, one-based trait-owner mapping, sharing mode,
owner filter parameters, sample counts, marker count, and native block
`start`, `size`, `n_kept`, `mu_min`, and `shrink`. Native starts stay
zero-based; formatted diagnostics additionally contain `start_1based`.

## 24. Memory

Runtime storage is reconstructed packed dense blocks, not retained low-rank
eigen factors.

## 25. Limitations

The route is serial, single-chain, BayesC-only, summary-level, and does not
support external-reference transport.

## 26. Evolution

External summaries require a separate, explicit contract for operator
rescaling, allele/scale transport, transformed-summary projection, reference
sample provenance, and runtime-diagonal semantics.
