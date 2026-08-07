# SBayesRC-S Phase 1 mathematical specification

## 1. Model identity and naming

SBayesRC-S means **SBayesRC with annotation Selection**. The suffix does not
refer to the `S` parameter in SBayesS. A future internal model identifier may
be `sbayesrc_selection`, but Phase 1 adds no production identifier, wrapper,
backend, or public API.

This document specifies and validates an isolated fixed-latent-variable
reference posterior. It does not claim novelty, production readiness, or
superiority to standard SBayesRC.

## 2. Relationship to standard SBayesRC

Standard SBayesRC continuously includes every supplied non-intercept
annotation. SBayesRC-S is a different posterior model: relevance of each
non-intercept annotation is uncertain and summarized by an annotation
posterior inclusion probability (annotation PIP).

Standard SBayesRC remains the committed validated baseline. No standard
SBayesRC production source, prior, likelihood, RNG path, or result contract is
changed by this reference work.

## 3. Shared annotation inclusion model

For marker (i), probit stick (k=1,\ldots,K), and (J) selectable
non-intercept annotations,

\[
z_{ik}=\alpha_{0k}+\sum_{j=1}^{J} A_{ij}\alpha_{jk}+\epsilon_{ik},
\qquad \epsilon_{ik}\sim N(0,1),
\]

with continuation (d_{ik}=I(z_{ik}>0)). The intercept is always included.
Each annotation has one selection state shared across all sticks,

\[
\delta_j\in\{0,1\},\qquad
\alpha_{jk}\mid\delta_j,\tau_k^2 =
\begin{cases}
0,&\delta_j=0,\\
N(0,\tau_k^2),&\delta_j=1.
\end{cases}
\]

Sharing δ pools evidence for relevance across the hierarchy while allowing
the included coefficients to differ in sign and magnitude by stick. This is
attractive because stick 1 carries the null-versus-active information and
later sticks have smaller eligible sets. It is a scientific hypothesis for a
future full model, not a claim that full-model mixing will improve.

## 4. Phase 1 priors

Phase 1 fixes

\[
\delta_j\stackrel{ind}{\sim}\operatorname{Bernoulli}(\pi_A),
\qquad 0<\pi_A<1,
\]

and fixes every slab variance τ²_k. The intercept uses the same common
flat-prior convention as the supplied standard-SBayesRC mathematical
reference. There is no selection indicator for the intercept.

Phase 1 does not sample `pi_A`, `tau2`, or `sigmaSqAlpha`.

## 5. Fixed-z target

The hard reference target is

\[
p(\delta,\alpha_0,\alpha\mid z,A,\pi_A,\tau^2).
\]

Latent (z) and the nested eligible set for each stick are fixed. There are
no marker effects, LD matrices, BayesR allocation updates, or variance
components in this target. This isolation is intentional: it qualifies the
new selection mathematics, not a future genomic sampler.

## 6. Exact model enumeration

For small (J), all (2^J) binary selection states are enumerated. With
(J=3), the state labels are `000`, `100`, `010`, `110`, `001`, `101`,
`011`, and `111` in the R reference's deterministic ordering.

For a selected set (S), the prior is

\[
p(S)=\pi_A^{|S|}(1-\pi_A)^{J-|S|}.
\]

Because the same state is shared across sticks,

\[
p(S\mid z)\propto p(S)\prod_{k=1}^K p(z_k\mid S).
\]

Weights are accumulated in log space and normalized by log-sum-exp.

## 7. Flat-intercept marginal likelihood

For stick (k), let

\[
Z_{S,k}=[\mathbf 1,A_{S,k}],\qquad
B_{S,k}=\operatorname{diag}(0,D_{S,k}^{-1}),
\]

where (D_{S,k}=\tau_k^2I). Define

\[
P_{S,k}=Z_{S,k}^{\mathsf T}Z_{S,k}+B_{S,k},\qquad
h_{S,k}=Z_{S,k}^{\mathsf T}z_k.
\]

Integrating the common flat intercept and the proper Gaussian slopes gives,
up to a stick-specific constant common to every model,

\[
\log p(z_k\mid S)=C_k-½\log|D_{S,k}|-½\log|P_{S,k}|
+½ h_{S,k}^{\mathsf T}P_{S,k}^{-1}h_{S,k}.
\]

For the empty model, the design contains only the intercept and there is no
log-determinant term for (D). The common improper-intercept constant is
the same for all selection states, so posterior model ratios are defined.

The implementation uses a Cholesky factor of (P) for its solve and
log-determinant. It never forms a determinant directly.

### Correction to the provisional local reference

The supplied provisional SBayesRC-S files conditioned exact enumeration on
fixed intercept values. That is a different conditional target from the
specified Phase-1 posterior. The tracked reference instead integrates the
flat intercept jointly with selected slopes. The provisional optional redraw
also sampled slopes conditional on a fixed intercept; the corrected redraw
samples the intercept and selected slopes jointly.

## 8. Collapsed delta update

For annotation (j), condition on the intercept and all other annotation
coefficients. At stick (k), define

\[
r_{jk}=z_k-\mathbf1\alpha_{0k}
-\sum_{\ell\ne j}x_{\ell k}\alpha_{\ell k},\quad
s_{jk}=x_{jk}^{\mathsf T}x_{jk},\quad
t_{jk}=x_{jk}^{\mathsf T}r_{jk}.
\]

Integrating α_jk under its (N(0,\tau_k^2)) slab gives

\[
BF_{jk}=(1+\tau_k^2s_{jk})^{-1/2}
\exp\left\{\frac{\tau_k^2t_{jk}^2}
{2(1+\tau_k^2s_{jk})}\right\}.
\]

Since δ_j is shared,

\[
\operatorname{logit}P(\delta_j=1\mid-)
=\operatorname{logit}(\pi_A)+\sum_k\log BF_{jk}.
\]

The log Bayes factor was independently checked against the dense Gaussian
marginal ratio with covariance (I+\tau_k^2xx^{\mathsf T}). The maximum
absolute discrepancy on the qualification fixture was
(2.00\times10^{-15}).

## 9. Conditional alpha update

After drawing δ_j=1, each stick-specific coefficient is regenerated from

\[
V_{jk}=(s_{jk}+\tau_k^{-2})^{-1},\qquad
m_{jk}=V_{jk}t_{jk},\qquad
\alpha_{jk}\sim N(m_{jk},V_{jk}).
\]

When δ_j=0, every α_jk is set exactly to zero.

## 10. Blocked Gaussian redraw

After a complete selection sweep, each stick redraws the intercept and all
currently included slopes jointly:

\[
\theta_{S,k}\mid z_k,S\sim
N(P_{S,k}^{-1}h_{S,k},P_{S,k}^{-1}).
\]

This is an exact Gibbs step. In the fixed-(z) reference there are only
(2^J) models, so their Cholesky factors are cached. Caching is a computational
optimization only and does not change the transition.

## 11. Partially collapsed ordering

For a single annotation, the Bernoulli draw integrates only that annotation's
coefficient block. The coefficients are then immediately regenerated from
their conditional distribution. Together these operations sample the joint
full conditional of
((\delta_j,\alpha_{j1},\ldots,\alpha_{jK})), conditional on the current
intercepts and other annotation coefficients. No marginalized value is used
after this joint update.

After all annotation blocks are updated, the complete coefficient redraw is
an ordinary Gibbs step conditional on the final δ state. The composition of
these invariant kernels preserves the fixed-(z) posterior.

Global model marginal odds and the single-annotation Gibbs odds are distinct:
the former integrate the intercept and every selected coefficient, while the
latter condition on the current intercept and other coefficients. The tests
check both contracts without equating them.

## 12. Exact annotation PIPs and alpha summaries

The annotation PIP is

\[
P(\delta_j=1\mid z)=\sum_{S:j\in S}P(S\mid z).
\]

Within a model,

\[
\theta_{S,k}\mid z_k,S\sim N(m_{S,k},V_{S,k}).
\]

The oracle averages the within-model means over exact model probabilities,
using zero for excluded slopes. It reports unconditional alpha means,
intercept means, and slope means conditional on inclusion.

## 13. Exact continuation and mixture-probability summaries

For marker (i), stick (k), and model (S), the posterior linear predictor
is Gaussian with mean μ_ik,S and variance v_ik,S. The probit-normal identity

\[
E\{\Phi(\eta)\}=\Phi\left(\frac{\mu}{\sqrt{1+v}}\right)
\]

gives the exact posterior mean continuation probability. Conditional on (S),
sticks are independent, so exact means of the stick-breaking products are
products of their per-stick expectations. The four-component mapping is

\[
\pi_{i0}=1-q_{i1},\quad
\pi_{i1}=q_{i1}(1-q_{i2}),\quad
\pi_{i2}=q_{i1}q_{i2}(1-q_{i3}),\quad
\pi_{i3}=q_{i1}q_{i2}q_{i3}.
\]

Every exact and sampled row is finite, lies in [0,1], and sums to one within
floating-point tolerance.

## 14. Reference fixture

The deterministic primary fixture has 180 marker-like observations, three
selectable annotations (`enriched_binary`, `continuous_signal`, and
`null_annotation`), and three nested eligible sets of sizes 180, 90, and 45.
The generating selection state is `(1,1,0)`, with fixed π_A=0.35 and
τ²=(0.8,0.8,0.8). The generating slopes are modest so several of the
eight posterior models retain appreciable probability; truth recovery is not
used as a qualification criterion.

Four chains use seeds 20260901--20260904, 12,000 iterations, 2,000 burn-in,
and initial states `000`, `111`, `101`, and `010`.

## 15. Phase-1 qualification result

The exact model probabilities were:

| State | Probability |
|---|---:|
| 000 | 0.154089 |
| 100 | 0.370734 |
| 010 | 0.133743 |
| 110 | 0.320583 |
| 001 | 0.004619 |
| 101 | 0.009951 |
| 011 | 0.001978 |
| 111 | 0.004304 |

Exact annotation PIPs were 0.705571, 0.460608, and 0.020852. Pooled MCMC
PIPs were 0.712000, 0.463025, and 0.020625. The qualification metrics were:

| Metric | Result | Gate |
|---|---:|---:|
| Maximum annotation-PIP error | 0.006429 | 0.020 |
| Model total-variation distance | 0.006655 | 0.030 |
| Maximum model-state probability error | 0.005542 | 0.020 |
| Maximum unconditional alpha-mean error | 0.004232 | 0.035 |
| Maximum conditional-inclusion alpha-mean error | 0.003157 | 0.035 |
| Maximum intercept-mean error | 0.001204 | 0.035 |
| Maximum q-mean error | 0.001469 | 0.012 |
| Maximum mixture-probability mean error | 0.001469 | 0.012 |

Chain-specific PIPs overlapped closely. Per-chain switching counts were
2,689--2,769 for the enriched annotation, 4,977--5,040 for the continuous
annotation, and 380--442 for the null annotation. The four-chain primary run
took 14.98 seconds on the recorded development machine.

Additional guards passed:

- a zero annotation column had exact PIP 0.35, equal to π_A, and sampled
  PIP 0.3439;
- exact column-permutation error was (1.11\times10^{-16}), with maximum
  permuted-chain PIP error 0.00771;
- duplicated annotations had exact symmetry error
  (5.55\times10^{-17}), with sampled PIP difference 0.0171.

The standalone reference is `tools/sbayesrc_s_reference.R`. It writes only a
compact ignored summary under `results/local/sbayesrc_s_reference/`. Installed
package tests use the separate independent helper
`tests/testthat/helper-sbayesrc-s-reference.R` and never source `tools/`.

## 16. Qualification gates

The hard reference gates are exact probability normalization, independent
Bayes-factor agreement, maximum PIP error at most 0.02, model-distribution TV
at most 0.03, maximum model-state error at most 0.02, alpha/intercept mean
errors at most 0.035, and q/π mean errors at most 0.012. Zero-information,
permutation, and duplicate-column guards are also required.

The exact oracle, not simulated truth recovery or R-hat, is decisive.

## 17. Deferred pi_A hierarchy

A future phase may place

\[
\pi_A\sim\operatorname{Beta}(a_\pi,b_\pi)
\]

and use its beta conditional or a collapsed beta-binomial model prior. Phase 1
does not implement either extension.

## 18. Deferred slab-variance hierarchy

Future work may sample τ²_k (or a carefully specified selection-model
analogue of `sigmaSqAlpha`) conditional on included coefficients. Phase 1
keeps all slab variances fixed to avoid confounding inclusion correctness with
variance learning.

## 19. Deferred stick-specific selection and production boundary

Stick-specific δ_jk, grouped annotations, annotation clusters, and
hierarchical inclusion probabilities are outside Phase 1. No C++ file, Rcpp
export, public model option, or production wrapper is added. A future backend
must use an explicit selection-model identity and remain separate from
standard SBayesRC.

## 20. Later evaluation plan

Phase 2 should couple this shared-δ model to observed continuation outcomes
through Albert--Chib latent updates, still in standalone R and still without
C++. Only after that small posterior is trusted should production integration
be designed. Subsequent scientific behavior belongs in `sblrbench`; no Study
07 work is part of Phase 1.

## Reference provenance

All five ignored files under `local_reference/` were read completely. The
continuous-alpha and recovery files establish the standard SBayesRC baseline;
the Jian Zeng 2024 script supplies external implementation provenance and the
known total-marker intercept detail; the two SBayesRC-S files are provisional
research inputs. They remain ignored and unmodified.

The Jian Zeng reference restricts later-stick rows to eligible markers but
uses total marker count `m` in its displayed intercept coordinate update.
Phase 1 does not adjudicate that separate source contract. Its exact reference
uses the explicitly specified flat-intercept Gaussian integration on each
eligible set.
