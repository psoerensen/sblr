# Posterior-preserving block SBayesRC probit sandwich transition

## Augmented conditional

For one populated stick, let `X` be the eligible annotation design, `d` the
continuation outcomes, `z` the Albert--Chib latent vector, and `alpha` the
stick coefficients. At the ordinary BayesRC endpoint,

\[
 z\mid\alpha,d \sim N(X\alpha,I)
\]

restricted to the orthant encoded by `d`. Conditional on the current
`sigmaSqAlpha`, write the proper Gaussian prior as

\[
 \alpha\sim N(m,B^{-1}).
\]

`B` contains the resolved proper intercept precision and the existing
non-intercept precision `1 / sigmaSqAlpha`. Define

\[
 P=X^T X+B,\quad u=X^Tz,\quad h=u+Bm.
\]

The ordinary exact blocked conditional is

\[
 \alpha\mid z,d,\sigmaSqAlpha\sim N(P^{-1}h,P^{-1}).
\]

Allocations enter by fixing the orthant and the eligible rows. Their strong
feedback with alpha is mediated partly by the persistent latent magnitudes:
ordinary Albert--Chib draws `z | alpha,d`, then alpha is regressed on that same
latent scale.

## Scale sandwich

Positive scaling preserves every latent sign. Integrating alpha from the
augmented density gives, up to constants independent of `z`,

\[
 \log f_Z(z)=
 -\tfrac12\{z^Tz+m^TBm-h^TP^{-1}h\}.
\]

Propose `ell' = ell + epsilon`, with `epsilon ~ N(0,s_px^2)`, on the
log-scale orbit and transform `z' = exp(ell' - ell) z`. The production step
starts each sandwich application at `ell = 0`. Symmetry on log scale and the
`n`-dimensional transformation Jacobian give

\[
 \log r = n\log g-
 \tfrac12\{a(g^2-1)-2b(g-1)\},
\]

where

\[
 a=z^Tz-u^TP^{-1}u,
 \qquad b=u^TP^{-1}Bm,
 \qquad g=\exp(\epsilon).
\]

Accepting with `min(1, exp(log r))` defines a reversible Markov kernel on the
latent orbit with invariant marginal `f_Z`. Drawing alpha exactly from its
blocked Gaussian conditional after this sandwich step yields the original
`p(alpha | d, sigmaSqAlpha)`. The existing
`p(sigmaSqAlpha | alpha)` update then completes an ordinary Gibbs composition
for the unchanged joint hierarchy posterior.

The argument is conditional on `sigmaSqAlpha`; no working-scale parameter is
retained and no model variance is introduced. Empty sticks retain their
existing posterior-valid prior-only update. The transition is initially
restricted to coupling one, so no tempered offset is involved.

## Correctness requirements

- the scale ratio must match direct evaluation of the integrated latent
  density plus its Jacobian;
- `P` must be positive definite and solved by Cholesky without inversion;
- the post-sandwich alpha draw must be the full multivariate Gaussian
  conditional, not one coefficient sweep;
- diagnostics may record scale acceptance and movement but must not consume
  RNG;
- the ordinary coefficient-wise path must remain byte-for-byte unchanged when
  the method is disabled.

This is a posterior-preserving sampler change, not the different
ApproxBayesRD annotation-selection model identified in the GCTB audit.
