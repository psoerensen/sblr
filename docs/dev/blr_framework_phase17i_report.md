# Unified BLR Framework Phase 17I report

## 1. Executive summary

Internal trait-specific MT CSR execution is active and shares the corrected
dense MT BayesC Gibbs implementation, result, finalizer, and positional adapter.

## 2. Repository baseline

Started clean on `master` at Phase 17H commit `cfbcb49820afba878f988276726229f2606deeb7`.
The fast workflow had the reviewed indentation defect; it was corrected first,
and both workflow files parsed. Hosted CI was not visible locally. R 4.4.1 and
Rtools44 compiled the package. Phase 17H baseline was 53 files, 369 blocks,
4,346 expectations, zero failures/warnings, and two opt-in skips.

## 3. Previous dense architecture

Dense-specific work comprised trait/marker dimensions, `ww[t][i]`, initial
`XXindices`/`XXvalues` residual reconstruction, and the same row traversal for
each marker difference. Statistics, RNG, accumulation, finalization, and schema
mapping were representation-independent.

## 4. Shared MT core architecture

`mtblr()` -> dense view -> `run_mt_default_core()` and internal adapter -> CSR
bundle -> `run_mt_csr_core()` both delegate to
`run_mt_bayesc_core_impl<DataView>()` -> `MtDefaultCoreResult` -> shared
finalizer -> shared legacy adapter.

## 5. Dense representation access

Dense overloads retain exact `ww`, diagonal-first row representation, stored
index/value order, initial rebuild nesting, and difference-subtraction order.

## 6. CSR representation access

CSR subtracts the double diagonal first, then calls the matching trait view's
float-to-double off-diagonal traversal. Vector overloads allocate nothing.

## 7. MT CSR contracts

`MtSparseLdBundleView { marker_count, trait_ld }` and `MtCsrDataView { wy, yy,
n, ld }` are borrowed, binding-neutral, and contain no priors, state, RNG,
files, outputs, or sharing assumption.

## 8. Validation

Validation requires positive dimensions, matching trait vectors and views,
positive sample sizes, valid canonical CSR views, model rows of trait width,
and in-domain set indices. It never reorders, flips, fills, or aligns biology.

## 9. Active marker conditional

The active `sampleBetaCPG_Mt_latent<DataView>()` is singular. Only diagonal and
residual-update access moved behind overloads; likelihood, Cholesky/jitter,
uniform/normal draws, state order, masking, and assignments are unchanged.

## 10. Statistical-order audit

```text
PHASE17H_MT_STATISTICAL_STATEMENTS=240
PHASE17I_MT_STATISTICAL_STATEMENTS=240
STATISTICAL_ORDER_EQUIVALENT=TRUE
RNG_CALL_ORDER_EQUIVALENT=TRUE
UPDATE_CALL_ORDER_EQUIVALENT=TRUE
RETENTION_POLICY_EQUIVALENT=TRUE
```

Differences are typed representation qualification, shared wrapper delegation,
and relocation of positional copying. No draw or update statement changed.

## 11. Internal native adapter

`mtblr_csr_internal(wy, ww, yy, b, ld_prefixes, sets, B, E, ssb_prior,
sse_prior, models, pi, nub, nue, updateB, updateE, updatePi, n, nit, nburn,
nthin, seed, method)` is Rcpp-registered but not namespace-exported or routed.

## 12. Shared mode

Length-one prefixes require exact diagonal equality, build one owner, and form
`nt` views. Unequal diagonals fail before disk build/MCMC.

## 13. Trait-specific mode

Length-`nt` prefixes build reserved independent owners and permit distinct
values and sparsity. Every resource must match the common marker count.

## 14. Common marker domain

All trait rows are already-aligned positions. There is no union likelihood,
missing-marker insertion, allele flip, strand inference, or reordering.

## 15. Study, ancestry, and overlap semantics

Traits may represent different studies/cohorts/ancestries and references.
Sample overlap remains an explicit future metadata/likelihood concern and is
not inferred from `E`.

## 16. Result, finalization, and adapter reuse

There is one core result, one final result, one finalizer, and one binding-neutral
20-position mapping. R formatting remains unchanged.

## 17. Dense reduction results

Shared fixed, shared all-updates, nonzero initialization, multiple sets,
three traits, trait-specific values/diagonals, and independent patterns all
matched exact structure and numerical tolerance `1e-12`. Invalid shared and
marker-count cases failed before sampling.

## 18. Scientific identities

Internal CSR output was finite, states binary, inclusion means bounded, fixed
B/E/pi preserved, and final residuals reconstructed trait by trait. Existing
Phase 17C owns covariance, probability, and retention identities.

## 19. Scalar reference preservation

Phase 3, 6, 8, 9C, 9E, 9G, and 10D permanent scalar owners are required in the
final validation; no scalar numerical source or fixture changed.

## 20. Dense MT reference preservation

Phase 17C corrected raw 3/3 and formatted 3/3 references passed after the core
and adapter changes.

## 21. Research-route protection

`mtblr_cpg_omp_csr()` and `mtblr_eigen()` are unchanged and unused by the new
route. No research storage, sampler, covariance policy, RNG, or OpenMP path is reused.

## 22. Public-interface protection

`sblr()` remains dense/default-only. Arguments, defaults, rejection behavior,
formatting, correlations, and the native 20 positions are unchanged.

## 23. Ownership, copies, and memory

At 10,000 markers and 1,000,000 symmetric entries, shared four-trait storage is
8,160,008 bytes; four independent owners are 32,640,032 bytes versus a 3.2 GB
dense test oracle. Per-chain CSR bytes are zero.

## 24. Performance

Five-repetition tiny signals completed in 0–0.02 seconds per fit; medians were
at timer resolution. These are regression signals, not speed claims. Moderate
memory scaling is audited analytically; completed-fit RSS is not peak RSS.

## 25. Generated wrappers

`R/RcppExports.R` and `src/RcppExports.cpp` gained exactly the internal wrapper
and registration. `NAMESPACE` did not change.

## 26. Mutation sensitivity

All seven in-memory structural mutations were detected: trait-zero reuse,
missing diagonal update, duplicate loop, unequal shared diagonals, traversal
change, public reroute, and duplicate positional adapter.

## 27. Tests and CI

Focused Phase 17I (34 expectations), Phase 17C/17E/17F/17H owners, all scalar
owners, and the fast tier pass. The complete ordinary suite passed: 54 files,
374 blocks, 4,380 expectations, zero failures/warnings, and two
expected opt-in skips (38.9 seconds wall time). Fast CI includes Phase 17I and
no environment variable was added. Extended reproducibility could not create
its `processx` child pipe on this Windows host (system error 5).

## 28. Deviations and blockers

The benchmark uses tiny executable workloads plus a moderate analytical memory
case because dense `O(Tm^2)` oracle construction is deliberately excluded from
routine maintenance runs. `devtools::check()` stopped before build at the same
pre-existing Windows `processx` pipe error (system error 5); compile/load and
ordinary tests passed independently. No numerical blocker is known.

## 29. Recommended next phase

Design and expose a canonical public `mtblr_csr()` interface aligned with
`stblr_csr()`, including R/Glist marker and allele validation,
trait/study/ancestry and LD-reference metadata, explicit marker-intersection
policy, and conversion to a stable named MT raw/fit schema, using the validated
internal CSR core without changing its numerical implementation.

## 30. Readiness marker

PHASE 17I COMPLETE — INTERNAL TRAIT-SPECIFIC MT CSR CORE ACTIVE
