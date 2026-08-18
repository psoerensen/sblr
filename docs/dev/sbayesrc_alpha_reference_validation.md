# Standard SBayesRC alpha reference validation

## Question

Is the standard SBayesRC alpha conditional itself mathematically correct and
estimable when the sequential-stick hierarchy is known? This note records a
small permanent package guard. It does not reopen Study 06, change the standard
model, or introduce a blocked production sampler.

## Reference provenance

The first available reference directory in the required search order was the
ignored repository-local `research/local_reference/` directory. The files read were:

| Role | File | SHA-256 |
|---|---|---|
| Compact mathematical reference | `test_sbayesrc_continuous_reference.R` | `bba54addbec93a2219f01bd0269714fa2a0032a1327092e5ee58732637b41c16` |
| Repeated isolated-alpha evidence | `test_sbayesrc_alpha_recovery_v2.R` | `9f1edd471d019aaf1e65a63be289abe47f5e162460c58ef5d4639730bcf082e3` |
| External/reference provenance | `sbayesrc_reference_jian_zeng_2024.R` | `50c393dd579459ce611fbc58f884d56a032531437b1d190411a5d7db842ab683` |

The files remain ignored and were not copied into the package.

## Exact mathematical target

For one stick, conditional on fixed Albert--Chib latent responses,

\[
z = X\alpha + \epsilon, \qquad \epsilon \sim N(0,I).
\]

The independent reference uses a flat intercept and
\(N(0,\tau^2)\) non-intercept priors. With

\[
B=\operatorname{diag}(0,\tau^{-2},\ldots,\tau^{-2}),
\quad P=X^\top X+B,
\]

the exact conditional is

\[
V=P^{-1},\qquad m=V X^\top z,\qquad
\alpha\mid z,X\sim N(m,V).
\]

The package default is the later, proper intercept generalization. Its first
prior-precision entry is finite and its prior-mean contribution is added to the
right-hand side. The test-only oracle supports this generalized term solely to
compare the current production formula; its primary fixed-\(z\) regression
retains the supplied flat-intercept reference target.

## Production source audit

The permanent statistical owner is
`src/st_bayesrc_annotation_prior.h`:

- `st_bayesrc_build_step_indicators()` defines continuation at zero-based
  stick `j` as `component > j`;
- stick 1 uses every marker, while later stick `j` is eligible exactly when it
  continued at `j - 1`;
- the same `idx` selects the outcome, annotation rows, latent responses, and
  coefficient cross-products;
- `st_bayesrc_scalar_conditional()` adds likelihood and prior precision once;
- the residual update is `residual += x * (old - new)`;
- the intercept uses the proper resolved prior by default and the exact zero
  precision only in explicit legacy-flat mode;
- only non-intercept coefficients enter the `sigmaSqAlpha` sum of squares;
- `st_bayesrc_compute_snp_pi()` applies the probit continuation map, the
  configured probability floor, and row normalization.

Packed-BED calls this owner from `src/blr_bed_bayesrc_core_impl.h`. CSR and
retained block-eigen routes call it from
`src/blr_csr_sbayesrc_core_impl.h`; block eigen changes the fitted likelihood
operator and residual policy, not the annotation hierarchy. The alpha
conditional is therefore shared rather than duplicated across these scalar
routes.

Annotation preprocessing in `R/stblr-bed-bayesrc-internal.R` recognizes a
single all-ones column, moves it to column 1 when necessary, rejects multiple
intercepts, and can add one at the R boundary. For an eligible set of size
\(n_k\), its direct cross-product is consequently
\(x_0^\top x_0=n_k\). The native regression fixture confirms eligible counts
8, 6, and 4 and continuation counts 6, 4, and 2 for a three-stick allocation.

For coordinate `j`, the audited implementation is

\[
V_j=(x_j^\top x_j+\tau^{-2})^{-1},\qquad
m_j=V_j\{x_j^\top r+(x_j^\top x_j)\alpha_j^{old}\},
\]

with the appropriate resolved intercept-prior terms for `j = 0`. Sign and
residual direction agree with the independent implementation.

**Production scalar alpha conditional matches independent oracle: YES.**

**Eligible-set intercept contract: PASS.**

## Package regression guard

`tests/testthat/helper-sbayesrc-alpha-reference.R` is an independent,
test-only implementation of:

1. the exact Gaussian posterior;
2. a coordinate Gibbs sweep with separate residual bookkeeping;
3. a precision-Cholesky blocked draw;
4. the sequential stick-to-component probability map.

`tests/testthat/test-sbayesrc-alpha-reference.R` checks that scalar and blocked
draws reproduce the same exact mean and covariance, compares the independent
coordinate conditional directly with the existing unexported native
fixed-latent hook, verifies later-stick eligibility and intercept counts, and
compares the independent probability map with `sbayesrc_marker_pi()`.

The existing `test-alpha-hierarchy-conditionals.R` and
`test-bayesrc-probit-intercept-prior.R` continue to guard the proper-intercept,
empty-stick, zero-information, separation, and `sigmaSqAlpha` contracts. No new
production/debug interface was needed. A second observed-stick recovery
simulation was deliberately not copied into package tests: the existing tests
cover executable Albert--Chib behavior, while repeated scientific recovery
belongs in `sblrbench`.

## Jian Zeng 2024 reference detail

The supplied reference restricts later-stick outcome and annotation records to
the eligible set but shows the total marker count `m` in the intercept update.
The sibling isolated-alpha comparison found only modest differences between
eligible-\(n\) and total-\(m\) updates on that fixture. This is provenance, not
a claim about production GCTB correctness. `sblr` retains the mathematically
matched eligible-set conditional.

## External scientific evidence and interpretation

The committed sibling Study 06 isolated-alpha addendum reports maximum
fixed-\(z\) mean and covariance errors of 0.00428 and 0.00084. Known-outcome
blocked and scalar posterior means differed by at most 0.00936. Across 20
Study-06-like simulated hierarchies, empirical 95% coverage was 0.90--1.00,
maximum R-hat was 1.00199, and uncertainty increased as eligible counts fell
down the hierarchy.

Repeated known-allocation simulations therefore show approximately unbiased
and calibrated alpha recovery across all sticks. Precision declines down the
hierarchy, but the alpha regression itself remains well behaved. The difficult
learned-alpha behavior identified in Study 06 concerns joint alpha--allocation
posterior exploration rather than an incorrect alpha conditional.

## Validation result

On the permanent 180-observation fixed-\(z\) fixture, the scalar Gibbs maximum
absolute mean and covariance errors were 0.00155 and 0.00045. The blocked
reference errors were 0.00147 and 0.00022. The direct native coordinate moments
matched the independent proper-prior formula to the test tolerance of
\(10^{-14}\), and its near-flat limit matched the supplied flat-intercept
reference. The eligible-set and probability guards passed.

The complete source-tree suite passed with no failures. The built-package suite
reported 4,809 passes, three expected skips, and the existing single MT
covariance warning. `R CMD check --no-manual --as-cran` completed with zero
errors, zero warnings, and the five pre-existing environment/package notes.

## Status

Standard SBayesRC remains unchanged and is the frozen validated baseline. The
new helper and tests consume RNG only inside the test process; production RNG,
defaults, output structure, priors, and likelihoods are unchanged.

SBayesRC-S is a separate future annotation-selection model and will have its
own model identity and implementation path.
