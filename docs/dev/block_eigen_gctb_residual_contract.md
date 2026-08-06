# Block-eigen GCTB-compatible residual contract

## Scope and pinned sources

This note fixes the scalar retained block-eigen residual-variance contract for
SBayesR and SBayesRC. It does not apply to BED, complete CSR, threshold/window
CSR, reconstructed dense block eigen, or multitrait samplers.

The primary implementation reference is `zhilizheng/SBayesRC` release v0.2.6,
commit `b95d3fcbad8ff358290922a58fff893439296138`, specifically
`src/BlockLDeig.cpp`, `src/SBayesRC.cpp`, `src/SBayesRC.h`, and
`R/sbr_eigen.r`. The GCTB cross-check is `jianzeng/GCTB` commit
`cc7fa7d765c83a89c6375946cf77fe50ba1a317e`, specifically
`scr/model.cpp`, `scr/model.hpp`, and `scr/data.hpp`.

## Factor and score units

For block $b$, let $R_b=U_b\Lambda_bU_b^\top$ be the retained eigendecomposition
of the standardized-marker correlation matrix and let $n_b$ be the GWAS sample
size. In the equal-diagonal standardized case, $X_b^\top X_b=n_bR_b$.

The pinned official code constructs

$$
Q_{b,o}=\Lambda_b^{1/2}U_b^\top,
\qquad
w_{b,o}=\Lambda_b^{-1/2}U_b^\top\hat\beta_b,
$$

where $\hat\beta_b=X_b^\top y/n_b$ on the official standardized-phenotype
scale. The official residual is
$r_{b,o}=w_{b,o}-Q_{b,o}\beta_b$.

`sblr` constructs the same retained eigenspace from the cross-product matrix.
Writing $D_b=\operatorname{diag}(X_b^\top X_b)$, its stored quantities are

$$
Q_{b,s}=\Lambda_b^{1/2}U_b^\top D_b^{1/2},
\qquad
w_{b,s}=\Lambda_b^{-1/2}U_b^\top D_b^{-1/2}X_b^\top y.
$$

Let $c_y$ convert official phenotype units to the units used by the supplied
`sblr` sufficient statistics. The pinned `nbsq` R wrapper sets the official
phenotype variance to one, whereas `sblr` deliberately retains the phenotype
scale encoded by `yy`; hence $c_y^2=V_{y,s}/V_{y,o}=V_{y,s}$ for that
comparison. Effects and transformed scores obey
$\beta_{b,s}=c_y\beta_{b,o}$ and
$w_{b,s}=c_y\sqrt{n_b}w_{b,o}$. The factor itself is phenotype-scale
independent. Thus, when $D_b=n_bI$,

$$
Q_{b,s}=\sqrt{n_b}Q_{b,o},\qquad
w_{b,s}=c_y\sqrt{n_b}w_{b,o},\qquad
r_{b,s}=c_y\sqrt{n_b}r_{b,o}.
$$

Consequently,

$$
Q_{b,s}^\top Q_{b,s}=X_b^\top X_b,
\qquad
c_y^2n_b\lVert r_{b,o}\rVert^2=\lVert r_{b,s}\rVert^2.
$$

When both implementations use the same phenotype units, $c_y=1$ and this
reduces to the simpler $n_b\lVert r_{b,o}\rVert^2=\lVert
r_{b,s}\rVert^2$ identity. The general `sblr` construction retains observed
cross-product diagonals rather
than replacing them by $n_b$. The official-unit conversion is therefore an
exact identity for the standardized same-sample construction and a diagnostic
crosswalk, not permission to overwrite observed diagonals.

## Marker conditional

The official marker update uses

$$
\frac{n_b}{V_{e,b}}\{q_{j,o}^\top r_{b,o}+\beta_j\}
$$

as its likelihood RHS and $n_b/V_{e,b}$ as its likelihood precision. Under the
conversion above,

$$
q_{j,s}^\top r_{b,s}+d_j\beta_j
=n_b\{q_{j,o}^\top r_{b,o}+\beta_j\},
\qquad d_j=q_{j,s}^\top q_{j,s}=n_b.
$$

Thus the existing `sblr` cross-product marker conditional is already the mapped
official conditional when it uses the block's $V_{e,b}$:

$$
\mathrm{rhs}_j=\frac{q_{j,s}^\top r_{b,s}+d_j\beta_j}{V_{e,b}},
\qquad
\mathrm{precision}_j=\frac{d_j}{V_{e,b}}+\mathrm{prior\ precision}.
$$

Within an `sblr` fit no factor, score, diagonal, or marker-effect rescaling is
required: $V_{e,b}$ is initialized and sampled in the phenotype units already
encoded by `yy`. When comparing output with the pinned `nbsq` wrapper,
official effects are multiplied by $c_y$ and official variance quantities by
$c_y^2$; heritability is invariant.

## Block variance quantities

The mapped official block genetic variance and effect sum of squares are

$$
V_{g,b,s}=c_y^2\lVert Q_{b,o}\beta_{b,o}\rVert^2
=\frac{\lVert Q_{b,s}\beta_{b,s}\rVert^2}{n_b},
\qquad
SS_{\beta,b}=\sum_{j\in b}\beta_j^2.
$$

The official-compatible residual scale in `sblr` units is

$$
S_{e,b}=\lVert r_{b,s}\rVert^2+\nu_es_{e,b},
\qquad
s_{e,b}=\frac{\nu_e-2}{\nu_e}V_y,
$$

and the draw is

$$
V_{e,b}^*=S_{e,b}/\chi^2_{q_b+\nu_e},
$$

where $q_b$ is retained block rank. The prior scale is computed once from the
initial phenotype variance. Multiplying the `sblr` residual norm by $n_b$ would
apply sample-size scaling twice and is incorrect.

The v0.2.6 compatibility safeguard accepts the draw only when
$V_{e,b}^*/V_y>0.7$; equality resets the block to $V_y$. This is an official
compatibility rule, not a generic clamp.

## Residual modes

The pinned v0.2.6 behavior is:

- `fixVe`: set every block to $V_y$ without sampling;
- `samVe`: sample every block every iteration;
- `allMixVe`: each iteration, sample a block only when
  $V_{g,b}\ge 10^{-8}$, $SS_{\beta,b}\ge 10^{-8}$, and
  $SS_{\beta,b}/V_{g,b}>\texttt{resam_thresh}$;
- `mixVe`: do not sample during zero-based iterations 0--49, select the
  permanent eligible block set at iteration 50 using the same strict tests,
  and sample that set from iteration 50 onward.

The strict `>` comparisons are retained. `updateE = FALSE` resolves to
`fixVe`. The default official-compatible mode for `updateE = TRUE` is
`allMixVe`, with `resam_thresh = 1.1` and `minimum_ve_ratio = 0.7`.

## Public policies and output meanings

Scalar retained block eigen exposes three policies:

- `gctb_block`: the block-specific contract above;
- `fixed_block`: $V_{e,b}=V_y$ for every block and iteration;
- `global_projected_legacy`: the previous global projected-residual sampler,
  preserved with its existing failure and RNG behavior.

`adjE` belongs only to `global_projected_legacy`; it is not part of the mapped
official marker conditional. Reconstructed-dense block eigen supports only the
legacy policy.

Under block policies, `fit$ves` is the retained arithmetic mean of block
residual variances. It is not a BED-equivalent global residual variance. The
official-compatible summary heritability is retained separately as

$$
h^2_{\mathrm{summary}}=\frac{\sum_bV_{g,b}}{V_y},
$$

not $V_g/(V_g+\overline V_e)$.
