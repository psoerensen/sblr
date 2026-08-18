# Annotation log-variance BLR implementation contract

> **Status: CURRENT QUALIFICATION RECORD.** This page records the qualified
> version-1 LV implementation and its reduction/trajectory gates. General
> architecture ownership and future restructuring are defined by
> [`README.md`](README.md) and the annotation-prior architecture documents.

## Frozen version-1 model

This implementation follows
`research/local_reference/annotation_log_variance_blr_development_plan_v1.md` and the
two executable `validate_bayesc_bayesr_logvar_*.R` references. The version-1
prior is fixed by default at

\[
\theta \sim N(0, 0.7^2 I).
\]

Users may supply another positive finite `theta_prior_sd`; it is not estimated.
There is no annotation intercept. For the processed marker-by-annotation design
\(X\),

\[
\eta_j=x_j^T\theta,\qquad q_j=\exp(\eta_j).
\]

Binary columns are centered without scaling. Continuous columns are centered
and divided by their sample standard deviation. Processing rejects non-finite
values, constant columns (including an all-ones intercept), and exact duplicate
or otherwise rank-deficient designs. Thus every column mean is zero,
\(\operatorname{mean}(\eta)=0\), the geometric mean of \(q\) is one, and the
global \(v_b\) remains the overall variance scale. R performs this transform
once; native code validates dimensions and finite values but does not normalize
again.

For BayesC-LV, \(d_j\sim\operatorname{Bernoulli}(\pi)\), \(b_j=0\) when
\(d_j=0\), and

\[
b_j\mid d_j=1 \sim N(0,v_bq_j).
\]

The inclusion probability remains global. Conditional on the current state,
the theta log likelihood (excluding its Gaussian prior) is

\[
-\tfrac12\sum_{j:d_j>0}\left[\eta_j+
\frac{b_j^2}{v_b\exp(\eta_j)}\right].
\]

For BayesR-LV, component probabilities remain global, \(\gamma_0=0\), and

\[
b_j\mid c_j=k>0 \sim N(0,v_b\gamma_kq_j).
\]

One theta vector is shared by every non-null component. Its conditional log
likelihood is

\[
-\tfrac12\sum_{j:c_j>0}\left[\eta_j+
\frac{b_j^2}{v_b\gamma_{c_j}\exp(\eta_j)}\right].
\]

Theta is updated jointly by Gaussian-prior elliptical slice sampling. When the
active set is empty it is drawn directly from its Gaussian prior. The sampler
does not clamp \(q\): non-finite theta, eta, likelihood, or exponential states
fail with a diagnostic error. Minimum and maximum observed log-q are recorded.

## Production layout and schedule

The shared annotation mathematics belongs in
`src/st_logvar_annotation_prior.h` and contains no CSR or retained-block
residual code. BayesC-LV and BayesR-LV use distinct type/core headers and
translation units named `blr_csr_logvar_*` and
`st_cpg_omp_csr_logvar_*.cpp`. Existing BayesC and BayesR likelihood,
marker-scale, variance-update, residual/operator, chain-seeding, LD-swap, and
aggregation conventions are reused without changing their scientific behavior.

The frozen iteration order is marker/component sweep, global probability
update, scale-aware \(v_b\) update, theta ESS update, rebuild \(q\), residual
variance update, and collection. The scale terms are \(b_j^2/q_j\) for BayesC
and \(b_j^2/(\gamma_{c_j}q_j)\) for BayesR.

## Output contract

The formatted fit retains the standard stable fields, including present-but-
`NULL` `ld_swap` and `ld_swap_chains` where unavailable, and adds `theta`,
`theta_summary`, `annotation_variance_ratio`, `annotation_transform`, and a
memory-safe posterior marker-scale mean. Extended convergence output may retain
theta traces, but marker-by-iteration q histories are not retained by default.
Theta summaries contain mean, SD, median, interval limits, R-hat, bulk and tail
ESS, and MCSE when supported. ESS diagnostics contain update counts, likelihood
evaluation and bracket-contraction means/maxima, and min/max log-q.

For a centered binary column, `exp(theta)` is the annotated-to-unannotated prior
variance ratio conditional on other annotations. For a standardized continuous
column it is the ratio for a one-SD increase.

## Validation gates

1. Native eta/q, both conditional log likelihoods, empty-active prior draws,
   numerical guards, ESS moments, theta-zero reductions, and fixed-q reductions.
2. BayesC-LV CSR: theta-zero ordinary-BayesC concordance, fixed-q fixed-marker
   concordance, learned-theta R-oracle agreement, then the full regression suite.
3. BayesR-LV CSR: corresponding ordinary/fixed-scale/learned-theta/component
   concordance, then the full regression suite.
4. Public `stblr_csr_annot(annotation_model = "log_variance")` dispatch and R
   preprocessing validation for both methods.
5. CSR versus effectively exact retained block-eigen concordance for both models.
6. Theta convergence, q stability, numerical, output-schema, and memory checks.
7. BED only if the existing scalar architecture permits narrow kernel reuse.

## Explicit version-1 non-goals

Version 1 excludes sampled `maf_effect_s` jointly with theta, component-specific
theta vectors, hierarchical or sparse/shrinkage annotation priors, multitrait
theta sharing, sparse annotation matrices, GPU work, and broad unrelated core
refactors. BayesC/R-LV is an alternative to SBayesRC, not a claimed replacement
or superiority result.

## Implemented endpoint and retained block-eigen support

The binding-neutral BayesC and BayesR engines now expose only an optional
marker-prior-scale provider and one post-`vb` hook. Their ordinary no-op
policies consume no random numbers. BayesC additionally exposes its CSR
operator and canonical raw conversion adapters. BayesR exposes a narrow policy
interface and canonical CSR Rcpp adapter. Ordinary same-seed trajectories were
frozen before extraction and reproduced exactly afterward; the full package
suite also passed after each extraction.

Separate CSR translation units implement BayesC-LV and BayesR-LV. Both use the
shared ESS/log-variance kernel. The production order preserves each established
sblr sampler's ordinary order and inserts theta/q after the scale-aware `vb`
update and before `ve`. This differs from the illustrative placement of the
global probability update in the independent R oracle, but every update remains
a full conditional and learned-theta posterior concordance was verified against
that oracle.

The public CSR interface preprocesses annotations exactly once in R and
dispatches `annotation_model = "log_variance"` with method `"sbayesc"` or
`"sbayesr"`. It returns the version-1 theta, variance-ratio, transformation,
marker-scale, convergence, and ESS diagnostic contract. Gates 0, 1, R1, 2, R2,
3, 4, and 6 are qualified by focused tests and the full regression suite.

Retained block-eigen construction now resolves Rcpp inputs once and passes a
plain `BlockEigenExecutionInput` to the shared constructor in
`src/st_block_eigen_execution.h`. That constructor owns no model policy and
continues to use the existing packed-BED decoder, retained/dense builders, and
operator representation. Thin internal `with_policy` entry points attach either
the ordinary no-op policy or an LV policy to the existing BayesC/BayesR block
engines. The public ordinary native signatures and R interfaces are unchanged.

Before extraction, deterministic two-chain low-rank trajectories were frozen
for ordinary BayesC and BayesR, including chain traces, residual outputs,
component quantities, block-Ve histories, and reduced-residual diagnostics.
Their normalized hashes remained byte-identical after extraction. Existing
dense/low-rank, residual-policy, block-Ve, rebuild, initial-scale, SBayesRC, and
operator-reduction tests also passed.

`stblr_block_eigen()` accepts `annotation_model = "log_variance"` with method
`"sbayesc"` or `"sbayesr"`. It reuses the CSR annotation preprocessor exactly
once and returns the same theta, q, transformation, convergence, and ESS
diagnostic contract. Effectively exact dense-block fixtures reproduce CSR theta,
q, component probabilities, and chain summaries exactly; remaining genomic and
variance differences are at floating-point operator-construction scale. The
canonical reduced-rank paths are covered separately. Gate 5 is qualified.

BED LV support remains deliberately unattempted and is not part of this Gate 5
implementation.
