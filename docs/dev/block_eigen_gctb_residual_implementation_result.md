# Block-eigen GCTB residual implementation result

## Decision

The scalar retained block-eigen SBayesR/SBayesRC implementation satisfies
decision gate **GCTB-R1**. Its default residual policy is `gctb_block` with
`block_ve_mode = "allMixVe"`. `fixed_block` is an explicit sensitivity policy,
and `global_projected_legacy` preserves the former experimental transition and
RNG behavior. BED and all CSR routes are unchanged.

The implementation was developed from `sblr` commit
`5fbc435f4db2f39b2351ef3e32316c6a4c4e28ad`, pinned SBayesRC v0.2.6 commit
`b95d3fcbad8ff358290922a58fff893439296138`, and pinned GCTB commit
`cc7fa7d765c83a89c6375946cf77fe50ba1a317e`. The read-only `sblrbench`
evidence commit was `536bde2a91b8497252e497d9f5cfc7fe189ed8c0`.

## Implemented contract

Each logical chain owns a residual vector and residual variance for every
retained block. Marker updates use their block's variance. A selected block
draw uses

```text
scale = ||r_b,sblr||^2 + nue * ((nue - 2) / nue) * Vy
df    = retained_rank_b + nue
Ve_b* = scale / chi_square(df)
```

and accepts the draw only when `Ve_b* / Vy > 0.7`; otherwise it resets that
block to `Vy`. The strict eligibility and ratio boundaries, the zero-based
`mixVe` trial interval, and all four official modes match the pinned source.
No global projected SSE is evaluated under a block policy.

The factor/score derivation is recorded in
`block_eigen_gctb_residual_contract.md`. In matched phenotype units,
`Q_sblr = sqrt(n) Q_official` and `r_sblr = sqrt(n) r_official`, so multiplying
the `sblr` residual norm by `n` would be wrong. The pinned `nbsq` wrapper uses
standardized-phenotype units, while `sblr` preserves the scale encoded by
`yy`: official effects scale by `sqrt(Vy)` and official `Ve` and `Vg` by
`Vy`. Heritability is invariant.

`fit$ves` remains present and is the retained mean block residual variance for
block policies. `fit$heritability_summary` is `sum(Vg_b) / Vy`. Compact output
also includes posterior mean and final block `Ve`, resampling/reset counts,
block sample sizes and ranks, policy metadata, and optional explicitly
requested retained block histories.

## Deterministic and regression validation

Independent tests cover factor, score, residual, marker-conditional, block
genetic-variance, effect-sum-of-squares, inverse-chi-square, and phenotype-unit
crosswalks. Exact boundary tests cover `Vg` and effect-SS at `1e-8`, the strict
`resam_thresh` comparison, the strict `0.7` minimum ratio, and the zero-based
`mixVe` selection iteration. One- and multi-block fixtures include overlapping
genotype-space projections. Diagnostics are RNG-neutral. Explicit
`global_projected_legacy` fits reproduce historical outputs exactly.

## Frozen large B0 validation

The frozen state had `n = 5,000`, `m = 37,991`, 76 blocks, and 1,945 initial
active markers. The legacy global projected scale reproduced `-4221.347`.
The 76 mapped block-local draw scales were all finite and nonnegative:

| statistic | block-local scale |
|---|---:|
| minimum | 754.807 |
| median | 900.628 |
| maximum | 1015.264 |

At iteration zero no block met the pinned `allMixVe` eligibility rule, so no
draw or minimum-ratio reset occurred and mean block `Ve` remained `2.045187`.
The registered 12-iteration, four-chain non-inferential SBayesR, fixed-alpha
SBayesRC, and learned-alpha SBayesRC smokes all completed with finite state.
Both annotation fits retained 12 alpha and 12 `sigmaSqAlpha` states per logical
chain. Their mean summary heritabilities were 0.5082, 0.5065, and 0.5612,
respectively. These are software smokes, not scientific fits.

## Pinned single-trajectory comparison

The 1,500-marker comparison is descriptive because the public official seed
interface does not provide independent native chains. All 100 modes were
retained in each 100-marker block. After applying the proved phenotype-scale
conversion (`Vy_sblr = 1.987863`):

| quantity | SBayesR | SBayesRC |
|---|---:|---:|
| PIP correlation | 0.9954 | 0.9760 |
| posterior-effect correlation | 0.9984 | 0.9973 |
| validation genetic-value agreement | 0.9985 | 0.9972 |
| effect RMSE after scale conversion | 0.00782 | 0.00845 |
| sblr summary heritability | 0.4887 | 0.4902 |
| official summary heritability | 0.4650 | 0.4609 |
| sblr mean block Ve | 1.98770 | 1.98761 |
| scale-converted official mean block Ve | 1.98765 | 1.98754 |
| sblr total block Vg | 0.9715 | 0.9745 |
| scale-converted official total block Vg | 0.9243 | 0.9161 |

The block-`Ve` scale is essentially identical and residual/heritability
behavior moves to the pinned contract without losing the strong effect/PIP
agreement. Active-count differences remain (176.6 versus 486.4 for SBayesR;
128.5 versus 78.9 for SBayesRC) and reflect other already documented sampler,
prior, and annotation-update differences rather than the residual-scale
mapping. Exact draw equality is neither expected nor claimed.

## Limitations and next step

Block `Ve` is not a BED-equivalent individual-level residual variance, and
`fit$ves` must not be compared across routes without its metadata. The official
minimum-ratio reset is retained solely for compatibility. This implementation
does not resolve annotation/allocation mixing or the remaining BED/block
calibration differences.

After review and a new package SHA, `sblrbench` should update only its
block-eigen method contract and semantic checkpoints, then run a separate
short design-validation review before any Study 06 requalification. This task
did not modify or run `sblrbench`.
