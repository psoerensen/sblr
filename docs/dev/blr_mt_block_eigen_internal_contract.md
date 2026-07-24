# Internal MT block-eigen contract

## 1. Purpose

This contract defines internal trait-specific multivariate BayesC execution over canonical block-filtered operators.

## 2. Public status

The route is internal only. No public function or algorithm selector exposes it.

## 3. Statistical target

It targets the corrected Phase 17C/17I MT BayesC model, not the scalar block-eigen sampler or legacy `mtblr_eigen()`.

## 4. Shared statistical core

Dense, CSR, and block-eigen wrappers call `run_mt_bayesc_core_impl<DataView>()`. There is one Gibbs loop and one active `sampleBetaCPG_Mt_latent<DataView>()` conditional.

## 5. Data contract

`MtBlockEigenDataView` borrows transformed `wy`, marginal `yy`, analysis sample sizes `n`, and one immutable `BlockEigenView` per trait.

## 6. Operator descriptor

Each internal descriptor contains `bed_files`, positive `n_bed`, 1-based positive `cls` per BED, optional 1-based `rows`, finite `af` in `(0,1)`, zero-based strictly increasing `block_start`, `eigen_filter`, and finite nonnegative `eigen_tau` and `eigen_eta`.

## 7. Fully shared mode

One descriptor builds one owner. All trait `wy` rows enter the builder together; all traits borrow the same owner. Analysis `n` may differ.

## 8. Trait-specific mode

With `nt` descriptors, each trait builds and owns its operator independently and retains its matching transformed `wy` row. Boundaries and numerical inputs may differ.

## 9. Mixed filters

Filter choice is per descriptor. Hard truncation projects only its matching rows; fixed and Ledoit–Wolf ridge preserve their input rows.

## 10. Marker domain

Descriptor position, summary position, initial effect position, output row, marker identity, and effect-allele orientation must already coincide. Native validation is dimensional, not biological.

## 11. Runtime diagonal

The only numerical diagonal is the canonical operator diagonal recovered from packed filtered floats. The adapter has no `ww` argument.

## 12. Hard-truncate data transformation

Projected `wy` is mandatory and is used for ranking, residual construction, covariance calculations, and output.

## 13. Ridge data transformation

Both ridge modes retain the original `wy`.

## 14. Ownership

Packed BED data and builder workspaces are transient. Adapter-owned `BlockEigenStorage` objects outlive all borrowed views and core/finalizer execution. There are no per-chain copies.

## 15. Disk I/O

BED reads, filtering, and eigendecomposition finish before the core creates its RNG. MCMC performs no disk I/O or filtering.

## 16. Result path

The route reuses `MtDefaultCoreResult`, `finalize_mt_default_result()`, and `make_mt_default_legacy_result()`; the latter receives transformed `wy`.

## 17. Dense oracle

Tests expand inspected packed blocks into dense rows: diagonal first, then increasing global off-diagonal positions. A union with explicit zero padding handles differing trait boundaries without entering production code.

## 18. CSR oracle

Diagonal blocks and CSR-representable filtered blocks reduce to the canonical Phase 17I route at float-aware tolerances.

## 19. Study/reference semantics

`n` is the analysis sample size; `n_bed` is the reference BED sample count. They may differ. This numerical route does not establish transportability.

## 20. Sample overlap

Only marginal `yy` is consumed. Sample overlap, `SSY`, and sampling-error covariance are not modeled.

## 21. External summaries

Hard-truncation projection is provenance-sensitive. Marker, allele, study, ancestry, reference, and overlap validation are outside this internal route.

## 22. Research-route exclusion

Legacy `mtblr_eigen()` has different storage, sampler, retention, and output logic. It is neither called nor used as an oracle.

## 23. Public Phase 17M requirements

A future public adapter must define biological alignment, reference/BED metadata, selected rows and columns, frequencies, blocks, filters, sharing policy, projection provenance, named raw output, fit formatting, and diagnostics.

## 24. Phase 17M public adapter

Phase 17M supplies those contracts through `mtblr_block_eigen()`. The internal
legacy symbol remains unchanged; `mtblr_block_eigen_raw_internal()` consumes
the same single adapter result and adds actual-build diagnostics.
