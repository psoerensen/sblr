# Unified BLR Framework Phase 17N Report

## 1. Executive summary

The future internal individual-level MT BayesC contract is implementation-ready.
Phase 17N specifies storage, likelihood, covariance, ownership, memory, output,
and permanent dense/BED oracles without implementing a sampler.

## 2. Repository baseline

Baseline commit `ffea4c0bf0ded29dd5c0c7267ecf0ffbbfa2835f` on `master`
had a clean tree. R 4.4.1 UCRT used Rtools44 GCC 13.2.0 and LAPACK
3.12.0; the runtime did not report a named external BLAS. Hosted CI was not
visible. The baseline source suite passed 4,605 assertions with zero failures
and warnings and two established opt-in skips. Built-tarball check passed with
the two existing NOTES.

## 3. Phase 17M verification

The current public summary family remains dense `sblr()`, CSR `mtblr_csr()`,
and same-BED block-filtered `mtblr_block_eigen()`. No public individual MT
route exists.

## 4. Packed-family history

Phase 15A rejected broad numerical-core consolidation. Phase 15B introduced
the conservative representation-identical `BedPackedGenotypeView` for BayesC
and BayesR while retaining the ownership-equivalent BayesRC view. Phase 17N
does not reverse those decisions.

## 5. Complete owner inventory

`PackedBedMatrix` is a move-only raw-pointer owner with rounded 64-byte marker
strides. `FastPackedBedMatrix` and `FastPackedBedMatrixBR` are vector owners
with compact `ceil(n/4)` strides. Allocation, padding, readers, consumers,
copy/move behavior, lifetime, threading, and owner tests are tabulated in the
internal contract.

## 6. Complete view inventory

BayesC/BayesR use `BedPackedGenotypeView`; BayesRC uses
`BedBayesRCPackedGenotypeView`. Both borrow immutable fit-owned storage. The
common view is retained unchanged for Phase 17O.

## 7. Physical BED contract

Only SNP-major PLINK BED with `6c 1b 01` is supported. Codes are
`00->2`, `01->missing`, `10->1`, `11->0`; four samples occupy a byte.
Partial-byte padding is never traversed as samples.

## 8. Standardization contract

Standardized values are `(g-2p)/sqrt(2p(1-p))`, with missing mapped to zero.
Raw dosage maps missing to `2p`. Frequencies count the allele represented by
dosage two; the repository does not justify labeling it minor/effect/alternate.
Phase 17O supports standardized genotypes and requires finite `p` in `(0,1)`.

## 9. Owner equivalence

The independent two-file, seven-sample R decoder matches native
`PackedBedMatrix` decoding exactly for nonidentity rows and marker columns.
Source audit proves the two vector readers use the same code extraction,
selection order, maps, and sample guards; existing BayesC/BayesR/BayesRC
references protect their decoded consequences. Direct vector-owner buffers are
private, so raw cross-owner buffer inspection would require a production/test
export and was deliberately not added.

## 10. Phase 17O owner decision

Use the existing move-only `PackedBedMatrix` and its existing reader as one
fit-owned owner. It avoids modifying or extracting protected scalar readers.
The nonblocked read and MinGW address-alignment caveats are explicit.

## 11. Sample and phenotype alignment

Current `.make_bed_marker_data()` behavior is fully inventoried. Phase 17O
requires one complete finite `n x nt` matrix, shared selected rows/order,
unique trait names, pre-adjustment, and centered columns. Scaling is not
performed.

## 12. Covariate decision

Phase 17O requires pre-adjusted phenotypes and has no covariate argument,
matching current `stblr_bed()` behavior. A public covariate policy is deferred.

## 13. Missing-phenotype decision

Missing phenotypes are unsupported. Masks and latent imputation would alter
scores, cross-trait likelihoods, sample sizes, `B/E` updates, and prediction.

## 14. Individual-level MT model

`Y=XB_eff+R`, `R_i~N(0,E)`, `beta_j~N(0,B)`,
`b_j=D_j beta_j`. Models, probabilities, latent/effective semantics, sets, and
covariance update order are the corrected MT contracts, not scalar BED logic.

## 15. Full-E marker conditional

For `s=x'R+w b`, `Omega=E^-1`, `P=B^-1`:

```text
C_k=P+w D_k Omega D_k
rhs_k=D_k Omega s
log q_k=log pi_k-.5 log|C_k|+.5 rhs_k'C_k^-1 rhs_k
beta_new~N(C_k^-1 rhs_k,C_k^-1), b_new=D_k beta_new.
```

The omitted `+.5 log|P|` is common. It makes the null-pattern weight exactly
`pi_null`. Inactive latent coordinates remain sampled but never enter the
sample residual directly.

## 16. Diagonal-E reduction

With diagonal `E`, the equations reduce exactly to the corrected summary-MT
conditional. Exact equality additionally requires identical `X'X`, `X'Y`,
marker residuals, state, models, sets, marker order, safeguards, and RNG.

## 17. Sample-residual operations

Use column-major `arma::mat R=Y-XB_eff`; always rebuild initially, calculate
`x_j'R` across traits, add back `w_j b_j`, and update by
`R <- R-x_j delta_j'`.

## 18. Marker decoding strategy

Decode one marker per visit into one reusable length-`n` double workspace.
This costs `8n` bytes and avoids repeated packed lookup across traits.

## 19. Marker order and sets

Use explicit disjoint complete sets, deterministic full sweeps, and the
corrected summary-MT marginal-score order with marker-index tie breaking.
There is no adaptive scheduling.

## 20. B update

Preserve guarded Phase 17C set-local `sampleBset`, latent inverse-Wishart
`sampleB_latent`, and final global `sampleB` order. Inactive latent effects enter
the latent update; only active effective effects enter the heuristic active
updates. Fixed `B` remains exact when disabled.

## 21. E update

Diagonal mode uses independent scaled inverse-chi-square draws from trait SSE.
Full mode uses `IW(nue+n, nue*sse_prior+R'R)`, symmetric SPD inputs, and the
binding-neutral form of the existing Bartlett/inverse-Wishart utility.

## 22. Genetic covariance

`U=XB_eff` and `G=U'U/n`, symmetrized. This matches scalar BED diagonal and
same-design summary identities.

## 23. LE/LD decision

Phase 17O reports `vgs/ves/vbs` and covariance matrices only. No cross-trait
LE/LD quantity is introduced without derivation.

## 24. CPO decision

CPO is unsupported and absent. Scalar marginal log-CPO is not a joint
multivariate predictive density.

## 25. RNG and execution

Serial, one chain, one fit-local `std::mt19937`, explicit seed, created after
validation/reading/maps/workspace and owned by the core. No OpenMP or static
mutable RNG state.

## 26. Ownership and lifetime

The adapter owns conversion and one genotype owner; the view borrows; core
state owns phenotype residuals, marker workspace, effects/states/covariances,
maps/order, and RNG; result owns summaries/final states; binding owns metadata.
There are no MCMC-time reads, handles, owner copies, or global caches.

## 27. Memory

The audit reports packed `m*64*ceil(ceil(n/4)/64)`, `Y/R=8nnt` each,
latent/effective effects `8mnt` each, state `4mnt`, workspace `8n`, marker maps,
model/covariance work, and traces. Analytical, transient, completed-fit RSS,
and peak RSS are explicitly distinct.

## 28. Complexity

Read/maps are `O(mn)`; a score/update is `O(nnt)`; a dense-pattern sweep is
`O(m(nnt+K nt^3))`; full `E` is `O(nnt^2+nt^3)`; reconstruction is
`O(mnnt+nnt^2)`.

## 29. Output schema

Reuse `mtblr_raw` version 1 with backend `mt_bed_bayesc` and data level
`individual`. Common marker meanings remain intact; no schema change is needed.

## 30. Marker-space outputs

Compute `wy=X'Y` once before MCMC and retain it. Do not retain marker-space
residual state during sampling; compute `r=X'R_final` once for output.

## 31. Sample-space outputs

Sample residuals, genetic values, fitted values, and individual predictions
remain internal and are not added to the raw schema in Phase 17O.

## 32. Dense-X oracle

Portable R tests independently verify BED codes/padding/order, standardized and
raw missing handling, `X'X`, `X'Y`, residual rebuild/update, genetic values,
`G`, diagonal SSE, full `R'R`, and marker-conditional matrices/weights.

## 33. Owner audit harness

The harness uses an independent BED writer/decoder and existing packed
inspection route. Private vector owners cannot be called directly without a
new native maintenance export; source equivalence plus permanent scalar-route
oracles is the conservative evidence boundary.

## 34. Existing-route protection

Phase 15A/15B, packed BayesC/BayesR/BayesRC, Phase 17C/I/J/K/L/M, scalar CSR,
raw schemas, BED interface, and backend consistency remain unchanged and pass.

## 35. Mutation sensitivity

All 18 required mutations are detected: code/missing/order/padding, ownership
and I/O, scalar/MT labeling, covariates/missingness, residual spaces, full-E,
schema decision, forbidden sampler/export/eigen base, and Phase 15B view.

## 36. Benchmark

The audit benchmark reports packed/dense/map/workspace bytes and tiny synthetic
decode, score, and update timings for small/moderate `nt`. These are regression
signals only and support no performance claim.

## 37. Installed-check behavior

R BED decoding, dense algebra, selection, residual, conditional, and contract
tests are portable. Only source architecture/repository assertions skip
narrowly without a source checkout.

## 38. Tests and CI

Focused Phase 17N and protected-owner tiers pass. The exact fast filter includes
17N and passed 1,812 assertions with zero failures or warnings and one
established opt-in skip. The final source suite passed 4,669 assertions with
zero failures or warnings and two established opt-in skips, both before and
after package checking. Built-tarball installed tests passed; `R CMD check`
reports installed tests as `OK` without an assertion total.

## 39. Package check

Built-tarball check has zero errors and warnings. The two unchanged NOTES are
the long Phase 17C fixture pathname and pre-existing compiled `std::cout`.

## 40. Diff hygiene

Changes are audit documentation, pure-R oracles/tests, audit/benchmark tooling,
and the fast CI filter only. No fixtures, generated wrappers, namespace,
schemas, production sources, or public R files change.

## 41. Deviations and blockers

Direct logical buffer comparison of the private vector owners is unavailable
without a new compiled inspection export. Independent decoding, identical
reader/map source contracts, and existing permanent numerical route oracles
provide the non-invasive equivalence evidence. This is not a Phase 17O blocker.

## 42. Recommended next phase

> implement the internal serial one-chain individual-level MT BayesC route defined by the Phase 17N contract, using one shared packed genotype owner/view, complete pre-adjusted phenotypes, sample-space residuals, corrected MT inclusion-pattern and B semantics, the selected residual-covariance policy, and the chosen common MT raw/finalization path, without adding a public interface.

## 43. Readiness marker

PHASE 17N COMPLETE — INDIVIDUAL-LEVEL MT PACKED-BED CONTRACT FORMALIZED
