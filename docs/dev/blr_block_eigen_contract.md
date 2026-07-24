# Canonical block-eigen contract

## 1. Purpose

This contract defines one block-filtered LD representation for the internal scalar family and future trait-specific MT reuse.

## 2. Terminology

“Block-eigen” names the construction/filtering technique. Runtime storage is not low rank: it is a float-packed upper triangle of each reconstructed dense filtered cross-product block.

## 3. Input provenance

BED files must be SNP-major. Files and their `cls` selections are concatenated in supplied order; `cls` is positive and one-based at the R/native boundary. Selected sample rows preserve supplied order. BED marker position, summary position, operator position, and posterior position share one marker/orientation domain. Native code does not reorder markers.

## 4. Standardization

BED codes map `00 -> 2`, `01 -> missing`, `10 -> 1`, `11 -> 0`. With allele frequency `p`, `z=(dosage-2p)/sqrt(2p(1-p))`; missing maps to standardized zero (mean imputation after centering). Valid `p` is finite and strictly inside `(0,1)`. Decoded `Z` is float, `A=ZZ'` accumulates in double, reconstructed entries narrow to float, the runtime diagonal is the corresponding float promoted to double, and `wy` remains double.

## 5. Block domain

Native starts are zero-based, nonempty, begin at zero, are strictly increasing, lie in `[0,m-1]`, and imply nonempty contiguous blocks ending at the next start or `m`. Blocks do not overlap and cover `[0,m)`. File/chromosome boundaries are not inserted automatically.

## 6. Unfiltered matrix

For each block, float `Z` produces double `A=ZZ'`, `D=diag(A)`, and `C=D^-1/2 A D^-1/2`; `C` is replaced by `0.5(C+C')` before filtering.

## 7. Hard truncation

The threshold is `max(tau,0.01)`. Eigenvalues at least the threshold are retained; if none qualify, the largest is retained. The reconstruction is `D^1/2 V_k diag(mu_k) V_k' D^1/2`. The matching summary vector is projected as `wy D^-1/2 V_k V_k' D^1/2`. Eigenvector signs/order are not contractual; the reconstruction and projection are.

## 8. Fixed ridge

`a=eta/(1+eta)`, clamped to `[0,1]`, and `tildeA=(1-a)A+a diag(diag(A))`. `wy` is unchanged. `tau` is ignored.

## 9. Ledoit–Wolf ridge

For standardized block columns, `C2=sum(C*C)`, `d2=C2-p`, each sample has `q_k=sum_i z_ik^2/d_i`, `bbar=sum_k q_k^2-C2/n`, and `a=min(bbar,d2)/d2` when `d2>0` (otherwise zero), clamped to `[0,1]`. The ridge reconstruction is used and `wy` is unchanged. `tau` and `eta` are ignored.

## 10. Runtime storage

`BlockEigenStorage` owns marker count, `BlockEigenBlockStorage` blocks, `block_of`, `local_of`, and the double runtime diagonal. Each block owns row-major upper-triangle floats. No eigensystem or retained-rank factor survives construction.

## 11. Borrowed view

`BlockEigenView` borrows immutable blocks, mappings, and diagonal. Armadillo and `std::vector<double>` overloads traverse increasing local positions identically. `apply_offdiag` excludes the diagonal; sampler code subtracts it once. `rebuild` applies the full within-block matrix and never crosses a block.

## 12. Validation

Builder-time checks cover dimensions, AF, starts, positive decoded diagonals, filter validity, and completed canonical storage. View-time validation checks non-null fields, coverage, packed lengths, finite values, exact mappings, positive diagonal, and equality to the float-packed block diagonal. Validation is outside marker loops. More expensive biological/content identity checks remain adapter responsibilities.

## 13. Diagonal semantics

Input summary `ww`, unfiltered reference-BED `diag(A)`, and filtered runtime diagonal are distinct. The marker ranking preceding construction retains current `wy`/input-`ww` behavior. Marker conditionals, residual rebuild, and LE calculations use the runtime operator diagonal. Filtered diagonal is not required to equal `ww`.

## 14. Marker ranking

All three scalar routes preserve their existing pre-sampling ranking and do not switch ranking to the filtered diagonal in Phase 17K.

## 15. Ownership

`PackedBedMatrix`, `Z`, `A`, `C`, reconstruction, and eigensystem are construction-only. One storage owner is completed before task creation. Trait/chain tasks borrow its immutable operations; they copy no operator. There is no MCMC-time BED I/O or eigendecomposition.

## 16. Diagnostics

Fields remain `start`, `size`, `n_kept`, `mu_min`, `shrink`. Hard truncation reports zero-based start, block size, retained count, minimum unfiltered correlation eigenvalue, and `1-n_kept/size`. Ridge reports `n_kept=size`, preserves the legacy `mu_min` value, and reports ridge weight `a` as shrink. Filter, tau, and eta remain separate input metadata; the fixed floor is 0.01.

## 17. Scalar model use

Internal experimental BayesC, BayesR, and SBayesRC factories call the same builder and store the canonical owner. Its public-compatible methods delegate to the borrowed view. Model-specific marker samplers remain unchanged.

## 18. Current limitations

Routes are internal/experimental, share one operator across scalar traits/chains, require existing equal-sample/input-ww agreements, reject LD-swap, have no scheduled execution or persisted operator format, and are not public selectors.

## 19. External summaries

The operator is built from a BED reference while `wy`/`ww` arrive separately. Current native code proves positional dimensions and shared scalar inputs but not study/reference sample identity or allele provenance. Hard projection of external `wy` therefore needs an explicit future provenance contract; external-summary support must not be inferred.

## 20. Future MT sharing

Fully shared storage requires identity of BED content/rows, marker order/orientation, AF, blocks, filter parameters, reconstructed values, diagonal, and projection. Shared boundaries alone permit independent numerical storage. Fully independent blocks/reference/filters are representable by one view per trait, provided one canonical marker domain remains.

## 21. Future MT transformed data

Every hard-truncated trait must consume its matching projected `wy`; ridge traits consume original `wy`. Every trait consumes its matching runtime diagonal and filtered matrix. Phase 17L will require BED descriptors, rows, marker columns, AF, block boundaries, filter parameters, and sharing policy beyond the current public MT CSR contract.

## 22. Complexity

Storage is `sum_b s_b(s_b+1)/2` floats plus two `m`-integer mappings and `m` doubles. Construction uses `n*max(s_b)` floats and quadratic double workspaces/eigensystems. A marker update costs `O(s_b)`; sweep/rebuild cost `O(sum_b s_b^2)`. Hard truncation does not reduce runtime storage rank.

## 23. Noncanonical MT research route

`mtblr_eigen()` uses legacy sparse-row inputs, a separate Gibbs loop, marker conditional, covariance/retention logic, and positional output. It neither builds nor consumes this contract, is unsupported, is not an oracle, and is not eligible as the Phase 17L basis.

## 24. Evolution policy

New scalar/MT block-filtered execution must extend this owner/view contract rather than introduce a parallel representation. Format changes require numerical references, explicit versioning if persisted, and documented ownership/alignment migration.
# Phase 17L MT consumption

An internal MT adapter now consumes one canonical `BlockEigenView` and its matching builder-transformed `wy` per trait. Shared, shared-boundary-only, and independent operators are representable; public provenance remains future work.

# Phase 17M public consumption

The public route supplies only validated same-BED descriptors. Diagnostics are
captured during the canonical build, native starts remain zero-based, and the
formatted fit adds one-based diagnostic starts. No operator is rebuilt for
inspection.
