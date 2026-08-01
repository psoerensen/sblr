# Retained low-rank scalar block-eigen operator

## Status and scope

This document is the Stage A mathematical contract for a retained low-rank
likelihood operator used by scalar `sbayesc`, `sbayesr`, and `sbayesrc`.
Effects, activity indicators, mixture membership, marker-specific prior
probabilities, annotations, and marker order remain in SNP space. The operator
changes only the representation of the likelihood. Multivariate sampling,
block discovery, cross-block LD, and hybrid operator selection are excluded.

The implementation baseline is sblr commit
`26bc36db474e24801bc12c088b9fc5dd788608eb`. The external design references
are qgg commit `bfac8b2c388afb7ae1c88019bcfef8588f81aedb` and GCTB commit
`cc7fa7d765c83a89c6375946cf77fe50ba1a317e`.

## Current sblr reconstructed-dense design

`build_block_eigen()` decodes one standardized BED block at a time and forms
the double-precision marker-by-marker cross-product `A = X'X`. It normalizes
`A` to the correlation matrix `C = D^{-1/2} A D^{-1/2}`, where `D = diag(A)`.
Hard truncation currently retains correlation eigenvalues at least
`max(eigen_tau, 0.01)`; ridge modes modify the full block. The implementation
then reconstructs a dense marker-by-marker block, stores its upper triangle as
`float`, and projects the SNP-space score into the reconstructed range for hard
truncation. MCMC residuals remain marker-length score residuals. This is a
reconstructed dense approximation, not a retained-factor runtime operator.

The historical route and its absolute-threshold, fixed-ridge, and Ledoit-Wolf
policies must remain unchanged as `dense_reconstructed`.

## qgg retained-factor design

At the pinned qgg revision, `genomic_bayes.R` forms a marker correlation matrix
`B`, obtains an R `eigen()` decomposition in descending order, uses
`cumsum(values) / sum(values) < eigen_threshold`, and stores
`Q = diag(sqrt(lambda)) V'`, an eigen-by-marker matrix. The corresponding
transformed marginal statistic is `w = diag(lambda^{-1/2}) V' b_scaled`.
`src/gbayes.cpp::sbayes_reg_eigen()` maintains `r = w - Q beta`, computes a
marker score from `Q[,j]'r`, and updates the reduced residual by a scaled Q
column. qgg uses correlation/marginal-effect scaling plus an explicit marker
sample-size multiplier; it is not directly interchangeable with sblr's
cross-product scale. Its retained-subspace products, rather than raw
eigenvectors, are the comparison target because eigenvector signs and bases of
repeated eigenspaces are non-identifiable.

## GCTB retained-factor design and rank boundary

At the pinned GCTB revision, Eigen returns eigenvalues in ascending order and
the operator is eigen-by-marker. `Qblocks` stores
`diag(sqrt(lambda)) U'`; `wcorrBlocks` is initialized as
`diag(lambda^{-1/2}) U' b_GWAS` and then maintained as `w - Q beta`.
`whatBlocks` stores the fitted reduced value. The sampler keeps marker effects
and priors in SNP space and uses an explicit `nGWASblock / vareBlk` multiplier.
The exact scale and variance crosswalk is in
`stblr_low_rank_gctb_crosswalk.md`.

The retained low-rank operator follows the GCTB/SBayesRC eigenspace likelihood
strategy, represented in sblr cross-product units with a global projected
residual-variance contract. This statement covers retained-rank,
transformed-score, and marker-conditional compatibility; it does not claim
GCTB's block-specific variance procedure.

The pinned construction routine treats values strictly greater than `1e-10`
as positive, initializes the descending accumulation with the largest
eigenvalue separately, and walks the remaining values from largest to
smallest. Its ordinary boundary behavior is equivalent to retaining the
smallest leading set whose cumulative positive mass is **strictly greater
than** `eigen_prop`. Consequently, if `eigen_prop` is exactly a cumulative
boundary, the boundary eigenvector alone is insufficient and the next smaller
positive eigenvector is retained. The largest eigenvalue is counted in mass
but separately from the loop counter.

The pinned code has unsafe/anomalous fallback indexing when no cumulative tail
is at or below the requested proportion, and it assumes at least two input
eigenvalues. sblr does not reproduce out-of-bounds behavior. The production
policy is the well-defined strict-mass rule above, with minimum retained rank
one, `0 < eigen_prop < 1`, and a clear error for positive rank zero. This agrees
with the pinned loop wherever its ordinary boundary is reached and preserves
its exact strict-boundary behavior. As `eigen_prop` approaches one, all
positive eigenvectors are retained once every proper leading cumulative mass
is less than or equal to the requested value.

## Canonical runtime scale

Production uses the **general cross-product scale**. For block `b`,

```text
A_b = X_b' X_b
D_b = diag(A_b)
C_b = D_b^{-1/2} A_b D_b^{-1/2}
```

BED dosages use `(dosage - 2p) / sqrt(2p(1-p))`; missing BED genotypes map to
zero and are therefore mean-imputed after standardization. The same selected
BED rows, columns, allele frequencies, and order must construct both the score
and operator. Cross-products and eigendecomposition are computed in double
precision.

Let positive eigenvalues be the finite values `lambda > 1e-10`, ordered
descending for selection. If `k` is the first rank for which

```text
sum(lambda[1:k]) / sum(lambda_positive) > eigen_prop,
```

retain those `k` eigenpairs. Define

```text
Q_b = diag(sqrt(lambda_bk)) V_bk' D_b^{1/2}
w_b = diag(lambda_bk^{-1/2}) V_bk' D_b^{-1/2} s_b
s_b = X_b' y.
```

`Q_b` has retained-rank rows and block-marker columns. It satisfies
`Q_b' Q_b = A_tilde_b` and `Q_b' w_b = s_tilde_b`, the score projected using
the same retained eigensystem. At full positive rank these equal `A_b` and
`s_b`, up to the numerical tolerance implied by the eigensolver and stored Q.
No standard errors are used to construct Q or w. GWAS sample size enters the
residual-variance degrees of freedom and reported variance scaling, not the
marker conditional, because Q and w already have cross-product units.
Reference sample size is implicit in `A_b`; the canonical route currently
requires same-sample/by-construction summary provenance for `updateE = TRUE`.
External summary statistics can define an algebraic projected score but do not
by themselves establish a coherent `yy - w'w` constant or same-sample
likelihood, so low-rank `updateE` must reject unresolved external provenance.

The exact mapping to GCTB's normalized correlation scale is

```text
Q_sblr = Q_gctb D^{1/2}
w_sblr = w_gctb                    when b_GWAS = D^{-1/2} s,
```

with sblr likelihood precision one in cross-product units. For the common
standardized case `D = n I`, `Q_sblr = sqrt(n) Q_gctb`; GCTB instead leaves Q
normalized and multiplies marker likelihood terms by `nGWASblock / vareBlk`.
Mixing one Q/w scale with the other's conditional is invalid.

## Reduced residual and marker conditionals

Each chain owns one reduced residual per block:

```text
r_b = w_b - Q_b beta_b.
```

For local marker `j`, with `q_j = Q_b[,j]`,

```text
d_j = q_j' q_j
u_j = q_j' r_b + d_j beta_j.
```

`u_j` is the unscaled score conditional on all other effects. After drawing a
new SNP-space effect,

```text
delta_j = beta_new - beta_old
r_b <- r_b - q_j delta_j.
```

For prior variance `v_j` and residual variance `v_e`, the non-null Gaussian
conditional has precision `d_j / v_e + 1 / v_j`, mean
`(u_j / v_e) / precision`, and variance `1 / precision`. Existing BayesC,
BayesR, BayesRC, MAF scaling, mixture, annotation, and RNG policies supply
`v_j` and prior odds unchanged. Each marker operation costs O(retained rank of
its block).

The operator must expose `diag`, `corrected_rhs`, `apply_difference`, `rebuild`,
`quadratic_form`, `projected_score_dot`, `residual_norm_squared`, and
`fitted_norm_squared`. A periodic rebuild sets each residual directly to
`w_b - Q_b beta_b` without dense reconstruction.

## Residual-variance contract

Production selects **`sblr_global_projected_variance`**, not GCTB's blockwise
variance sampler. With compatible same-sample Q, w, score, and `yy = y'y`,

```text
SSE = yy - sum_b w_b'w_b + sum_b r_b'r_b
    = yy - 2 beta' Q'w + beta' Q'Q beta.
```

The first expression is the reduced-space implementation; the second proves
equivalence to sblr's existing `yy - beta'(s_tilde + residual_score)` update.
At full positive rank it is exactly `||y - X beta||^2` when the block-diagonal
operator is the complete X'X model. Under truncation it defines a coherent
projected likelihood: the retained coordinates are modeled by Q and the
orthogonal response sum of squares `yy - sum w'w` remains constant. sblr keeps
its existing inverse-chi-square degrees of freedom `n + nue` and prior scale
`nue * sse_prior`; it does not import GCTB's per-block degrees of freedom or
fallback gates. `updateE` cannot be enabled until the construction validates
finite nonnegative projected constants within a scale-aware tolerance.

GCTB instead computes, per block, `sse_b = nGWASblock_b * ||wcorr_b||^2`, uses
`df + numEigenvalBlock_b`, adds `df * scale`, draws a block residual variance,
and accepts it only when both `ssqBlocks/vargBlocks` and `sample/vary` pass its
gates; otherwise it assigns the fixed phenotypic variance. The global value is
the mean of block values. That complete contract is intentionally not combined
with sblr's global update.

## Genetic and reported variance quantities

The low-rank genetic variance is

```text
v_g = sum_b ||Q_b beta_b||^2 / n
    = beta' A_tilde beta / n.
```

This is the variance represented by the truncated operator; at full rank it is
the ordinary `beta' X'X beta / n`. The projected score product is
`beta' Q'w`. The linkage-equilibrium quantity is
`vle = sum_j d_j beta_j^2 / n`, and `vld = vg - vle`; both therefore describe
the retained operator when truncated. `vbs`, mixture/component effect-prior
variances, SNP PIPs, component probabilities, and annotation quantities stay
in SNP/prior units and do not rotate. `ves` follows the projected likelihood
above. `log_cpo` and posterior-predictive quantities that currently require a
marker-length residual are not automatically valid on reduced coordinates;
they must either receive a proved reduced-space definition or be explicitly
unavailable. No dense formula may be silently reused.

## Storage, ownership, and construction

Immutable fit-level storage owns block start/size/rank, column-major float Q
constructed from double intermediates, double w for each input trait/score,
double runtime diagonals, marker-to-block and local-offset maps, eigenvalue
diagnostics, scale/sample metadata, timings, workspace estimates, and storage
bytes. Each chain owns only its reduced double residuals plus existing effects,
states, and traces. Q is neither copied per chain nor returned in ordinary
scientific fits. A development-only inspector may expose Q and w for small
fixtures.

Construction decodes each block once, uses the current missingness and
standardization rules, forms A/C and decomposes them in double precision,
selects rank, constructs Q/w/diagonals, records diagnostics and timings, then
releases decoded genotypes and dense temporaries. MCMC does not reconstruct A,
reread BED, or repeat eigendecomposition.

Per block metadata records start, size, positive and retained rank, both rank
and mass fractions, positive/retained mass, smallest retained and largest
omitted values, negative count/mass, tolerance, storage, and construction
measurements. Negative mass is the sum of absolute values of eigenvalues below
`-1e-10`; values in `[-1e-10, 1e-10]` are numerical zero. Positive-rank-zero
blocks and non-positive retained diagonals are rejected.

## Alignment and block contract

The first public block begins at marker one. Starts are strictly increasing,
sizes are positive, and derived half-open ranges cover every marker exactly
once with no gap or overlap. Marker IDs/order, BED columns/files, allele
frequencies, and selected rows must match the summary provenance. Within-block
and whole-block permutations are valid only when scores, metadata, annotations,
effects, and BED provenance are permuted together. There is no cross-block
operator action.

## Equivalence, limitations, and extension points

Same-sample full-positive-rank Q and w reproduce A and s, so deterministic
conditionals, SSE, and genetic variance match the dense/CSR target within
float-aware tolerance. Truncated comparisons must reconstruct a dense oracle
from the exact stored Q, not independently recompute eigenvectors.

Low rank is beneficial only when retained rank is materially below marker
count. CSR can be preferable for sparse blocks, while packed dense storage can
be preferable for small dense full-rank blocks. External-summary and
reference-panel mismatch remains a scientific limitation.

The immutable factor/view and separate chain state are intentionally scalar-
trait-neutral building blocks for future MT reuse, but an MT likelihood and
residual-covariance contract must be derived separately. A future hybrid
planner may select CSR, packed dense, or retained low rank per block using
measured rank/sparsity and storage/runtime costs; no automatic selection is
part of this contract.
