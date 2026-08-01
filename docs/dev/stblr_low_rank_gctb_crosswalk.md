# sblr retained low rank and pinned GCTB crosswalk

This crosswalk is based only on GCTB commit
`cc7fa7d765c83a89c6375946cf77fe50ba1a317e`, specifically `scr/eigen.cpp`,
`scr/data.cpp`, `scr/data.hpp`, `scr/model.cpp`, and `scr/model.hpp`. It does
not infer unresolved quantities from names.

The retained low-rank operator follows the GCTB/SBayesRC eigenspace likelihood
strategy, represented in sblr cross-product units with a global projected
residual-variance contract. Retained rank, transformed score, and marker
conditionals are compatible under the scale mapping below. GCTB's
block-specific residual-variance procedure remains a distinct contract and is
not reproduced by the production sampler.

## Object map

| sblr object | Pinned GCTB object | Shape/orientation | Initialization and update |
|---|---|---|---|
| retained factor | `Qblocks[b]` | retained eigenvalues x block SNPs | GCTB: `diag(sqrt(lambda)) U'`; sblr: same correlation factor right-multiplied by `D^{1/2}` |
| transformed score | input `w_b` | retained eigenvalues | GCTB stores the initial value in `data.wcorrBlocks`: `diag(lambda^-1/2) U' b_GWAS`; sblr uses `diag(lambda^-1/2) V' D^-1/2 s` |
| chain residual | `wcorrBlocks[b]` in a model | retained eigenvalues | copied from data, then `w - Q beta`; update adds `Q_i(old-new)` |
| fitted reduced value | `whatBlocks[b]` | retained eigenvalues | zeroed each sweep and accumulated as `sum Q_i beta_i` for sampled non-null effects |
| SNP effects | `SnpEffects::values` | block SNPs in global marker order | never rotated; sampled marker by marker |
| block sample size | `nGWASblock[b]` | scalar per block | median/assigned GWAS N; multiplies normalized likelihood precision |
| retained rank | `numEigenvalBlock[b]` | scalar per block | number of rows of Q |
| block residual variance | `vareBlk.values[b]` | scalar per block | blockwise inverse-chi-square draw with gates; global `vare` is their mean |

## Scale

GCTB stores a correlation-scale factor

```text
Q_G = diag(sqrt(lambda)) V'
w_G = diag(lambda^-1/2) V' b_GWAS.
```

For `b_GWAS = D^-1/2 s`, sblr's general cross-product representation is

```text
Q_S = Q_G D^1/2
w_S = w_G.
```

Thus `Q_S'Q_S = D^1/2 V lambda V'D^1/2` and `Q_S'w_S` is the correspondingly
projected cross-product score. When `D = nI`, GCTB leaves Q normalized and uses
`nGWASblock/vareBlk`; sblr absorbs `sqrt(n)` into Q and uses `1/ve`. The two
marker likelihood precisions then agree if block N, diagonal N, and residual
variance refer to the same sample/units.

## Rank selection

Eigen's self-adjoint solver supplies ascending eigenvalues. GCTB counts mass
from the largest downward, with positivity `lambda > 1e-10`. The largest value
is initialized separately. At an ordinary boundary it retains the smallest
leading set whose cumulative positive mass is strictly greater than the
requested cutoff; equality retains one more eigenvector. Its fallback branch
contains unsafe/anomalous indexing for very small cutoffs and rank-one inputs.
sblr preserves the verified strict-boundary rule but defines the fallback
safely as minimum rank one and rejects rank zero.

## Marker update

GCTB uses normalized columns and computes

```text
rhs_G = (Q_i' wcorr + beta_old) * nGWASblock / vareBlk.
```

The `+ beta_old` assumes the normalized marker diagonal is one. sblr uses the
general form

```text
u_i = q_i' r + (q_i'q_i) beta_old
rhs_S = u_i / ve.
```

After sampling, both perform `residual += q_i(beta_old-beta_new)`. BayesC,
BayesR, and BayesRC differ only in the SNP-space prior variance/probability
terms; annotations and component memberships remain marker aligned.

## Rebuild and genetic variance

`ApproxBayesC::Rounding::computeWcorr_eigen()` resets GCTB residuals to the
stored data vectors and subtracts Q columns times SNP effects. sblr's rebuild
is exactly `r_b = w_b - Q_b beta_b`.

GCTB's `BlockGenotypicVar::compute()` sums `whatBlocks[b].squaredNorm()`.
sblr uses `sum ||Q_b beta_b||^2 / n`; the division is required because sblr Q
is in X'X cross-product units. These are the same correlation-scale quantity
after applying the scale mapping.

## Residual variance

Pinned GCTB computes, for each block,

```text
sse_b       = nGWASblock_b * ||wcorr_b||^2
df_tilde_b  = df + numEigenvalBlock_b
scale_b     = sse_b + df * scale
sample_b    = InvChiSq(df_tilde_b, scale_b).
```

It adopts `sample_b` only when
`ssqBlocks[b] / vargBlocks[b] > threshold` and
`sample_b / vary > 0.9`; otherwise it assigns `vary`. The reported global
residual variance is the mean of block values. This includes block-specific N,
rank degrees of freedom, prior scale, fixed fallback, genetic/effect-sum gates,
and block-to-global aggregation.

sblr intentionally selects `sblr_global_projected_variance`:

```text
SSE = yy - sum ||w_b||^2 + sum ||r_b||^2,
df  = n + nue,
scale = SSE + nue * sse_prior.
```

It is a coherent same-sample projected likelihood and is not described as
exact GCTB compatibility. No GCTB block-variance component is silently mixed
into it.
