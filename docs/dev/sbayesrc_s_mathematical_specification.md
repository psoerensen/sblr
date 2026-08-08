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

## Phase 2: observed continuation outcomes

### Status and scope

Phase 1 remains the validated fixed-z exact baseline. Phase 2 adds only the
Albert--Chib latent-variable layer for observed continuation outcomes. The
shared selection state, fixed `pi_A = 0.35`, fixed `tau2 = 0.8`, flat
always-included intercept, collapsed annotation update, and blocked Gaussian
coefficient redraw are unchanged.

Phase 3 hyperparameter learning, production APIs, C++, genomic marker effects,
LD, and Study 07 remain deferred.

### Observed-d target

For observation (i) eligible at stick (k),

\[
d_{ik}=I(z_{ik}>0),\qquad
z_{ik}\mid\eta_{ik}\sim N(\eta_{ik},1),
\]

where

\[
\eta_{ik}=\alpha_{0k}+\sum_{j=1}^J A_{ij}\alpha_{jk}.
\]

The normal distribution is truncated to z>0 when d=1 and to z<=0 when
d=0. The augmented target is

\[
p(z,\delta,\alpha,\alpha_0\mid d,A,\pi_A,\tau^2),
\]

and the required marginal target integrates z.

### Stick eligibility

Stick 1 uses all observations. Stick k+1 uses exactly the rows whose outcome
at stick k is one. A validation guard reconstructs each later eligible index
from the preceding outcome and requires exact identity. Outcomes, latent
variables, annotations, predictors, Bayes factors, and Gaussian redraws all
use that same index vector.

### Truncated-normal update

The base-R inverse-CDF implementation works on tail probabilities. For d=1,
it samples a standard-normal survival probability uniformly between zero and
the survival mass above -eta. For d=0, it samples a lower-tail probability
between zero and Phi(-eta). It then shifts by eta. Focused tests cover eta
from -6 to 6 and require finite draws with the exact observed sign.

### Gibbs schedule and invariance

Each iteration applies:

1. draw every z conditional on d and the current coefficients;
2. for each annotation j, sample the joint block
   (delta_j, alpha_j1:K) using the validated Phase-1 collapsed/regenerated
   conditional at the new z;
3. redraw each stick's intercept and all included slopes jointly from their
   exact Gaussian conditional.

Step 2 integrates only the coefficient block that it immediately regenerates.
No marginalized coefficient is conditioned upon before regeneration. Step 3
is an ordinary Gibbs update conditional on the completed selection sweep.
The composition therefore preserves the augmented Phase-2 target, and its
marginal transition preserves the observed-d posterior.

### Independent observed-d comparator

The primary sampler uses the Phase-1 rank-one log Bayes-factor helper and a
precision-Cholesky coefficient draw. The comparator independently evaluates
the inclusion marginal as

\[
-\frac12\log\tau_k^2
-\frac12\log(s_{jk}+\tau_k^{-2})
+\frac12\frac{t_{jk}^2}{s_{jk}+\tau_k^{-2}},
\]

and samples coefficients through a covariance-Cholesky construction. These
routes share only the unavoidable latent-normal generator and low-dimensional
Gaussian model construction. Their maximum primary-fixture differences were
0.00890 for annotation PIPs, 0.0312 posterior SD for alpha means, 0.0204 for
conditional alpha means, 0.00118 for intercept means, and 0.00285 for both q
and mixture-probability means.

### Nested-model bridges

With all delta values fixed to one, SBayesRC-S reduces to the corresponding
fixed-tau continuous-alpha probit hierarchy. Comparison with the committed
standard-SBayesRC R reference gave maximum differences of 0.00932 for alpha
means, 0.00380 for alpha SDs, 0.00712 for q means, and 0.00446 for mixture
probability means.

With all delta values fixed to zero, every non-intercept coefficient remained
exactly zero, intercept draws were finite, q depended only on the intercepts,
and mixture probabilities normalized to machine precision.

### Primary fixture and diagnostics

The deterministic fixture has 420 observations; three sticks with 420, 215,
and 128 eligible rows; continuation counts 215, 128, and 65; and the three
Study-06-analogue annotations. The generating selection state is (1,1,0).
Four chains start from `000`, `111`, `101`, and `010`, use independent seeds,
and retain 7,500 of 9,000 iterations.

Pooled annotation PIPs were 0.9991, 0.9579, and 0.0009. Chain-specific ranges
were 0.0012, 0.0145, and 0.0008. Maximum alpha and intercept R-hats were
1.00021 and 1.00014. Classical alpha ESS estimates ranged from 6,758 to
30,000. The intermediate continuous annotation made 92--108 transitions in
each direction per chain; near-certain annotations appropriately switched
less often.

The primary and comparator runtimes were 38.96 and 36.08 seconds. Every
retained q and mixture-probability summary was finite and normalized.

### Structural and signal guards

- An exactly zero column had pooled PIP 0.3525 versus prior 0.35; its four
  chain PIPs were 0.3620, 0.3598, 0.3470, and 0.3413.
- Permuting annotation columns changed corresponding PIPs by at most 0.000283
  and alpha means by at most 0.00265.
- Two duplicate columns had pooled PIPs 0.5089 and 0.5449. Their difference
  was 0.0360; chain-specific results exposed the expected redundant-model
  uncertainty rather than selecting an arbitrary truth.
- A strong-signal fixture initialized fully excluded reached PIPs 1.000 and
  1.000 for the informative annotations.
- A weak-signal fixture produced an intermediate PIP of 0.178 with repeated
  state switching.
- Across ten intercept-only datasets with exchangeable standardized null
  columns, mean PIPs were 0.0028, 0.0049, and 0.0265; the between-column range
  was 0.0237. Finite-data null evidence moved probabilities below the prior,
  without a persistent index-specific preference.

### Repeated-simulation calibration

Twenty independently regenerated datasets used 300 observations and four
chains each. Informative annotation PIPs had mean 0.9389, median 0.9986, 5th
and 95th percentiles 0.6641 and 1.0000; 95% exceeded 0.5 and 85% exceeded 0.9.
Null PIPs had mean 0.00631, median 0.00275, and 5th and 95th percentiles
0.00050 and 0.03376; none exceeded 0.5.

Conditional informative-alpha bias averaged 0.0153 and empirical 95% coverage
was 0.9667. Conditional null-alpha bias averaged 0.00673. One of 20 short
replicate fits exceeded the preregistered between-chain PIP-range diagnostic;
this was within the allowed development gate of two and is retained as a
finite-chain limitation. The repeated experiment took 105.49 seconds.

### Phase-2 decision

All Phase-1 regression tests, truncated-normal guards, eligibility checks,
primary/comparator agreement gates, nested-model bridges, zero/permutation/
duplicate guards, convergence checks, and repeated-simulation gates passed.

**SBS2-R1: the fixed-hyperparameter shared-delta SBayesRC-S posterior is
validated through the observed continuation/probit hierarchy in standalone
R.**

The next permissible task is Phase 3 reference work for inclusion-probability
and slab-variance learning. It must remain R-only until separately validated.

## Phase 3: hierarchical hyperparameter learning

### Status and production boundary

Phase 3 extends only the standalone observed-outcome R reference. Phase 1 and
Phase 2 remain validated and unchanged. Standard SBayesRC production code,
the public API, C++, Rcpp exports, and genomic marker-effect samplers are not
part of this phase.

Phase 3 was qualified sequentially:

1. **3A:** learn the global annotation inclusion probability while holding
   stick slab variances fixed;
2. **3B:** learn one slab variance per stick while holding the inclusion
   probability fixed;
3. **3C:** learn both hierarchies jointly only after 3A and 3B passed.

### Phase 3A: inclusion-probability hierarchy

The hierarchy is

\[
\pi_A\sim\operatorname{Beta}(a_\pi,b_\pi),\qquad
\pi_A\mid\delta\sim\operatorname{Beta}
\{a_\pi+M,b_\pi+J-M\},
\]

where (M=\sum_j\delta_j). Integrating (\pi_A) gives

\[
p(\delta)=\frac{B(a_\pi+M,b_\pi+J-M)}{B(a_\pi,b_\pi)}.
\]

Consequently, the collapsed prior log odds for annotation (j), before its
likelihood Bayes factors, are

\[
\log(a_\pi+M_{-j})-
\log\{b_\pi+J-1-M_{-j}\}.
\]

The reference implements two invariant routes. The explicit route conditions
each selection update on the current sampled (\pi_A), then draws (\pi_A) from
its beta full conditional. The collapsed route uses the beta-binomial odds,
then draws (\pi_A) conditionally only for posterior output. A chain never uses
both sampled and integrated (\pi_A) in one selection probability.

For fixed z and small J, the exact finite-state oracle replaces the independent
Bernoulli model prior with the beta-binomial expression above. Explicit and
collapsed chains reproduce the exact model probabilities and annotation PIPs.
The identity

\[
E(M\mid y)=\sum_jP(\delta_j=1\mid y)
\]

is a permanent numerical invariant.

The primary Phase-3A exact PIPs were 0.999914, 0.999637, and 0.002728.
Maximum exact-versus-MCMC errors were 0.000161 for the explicit route and
0.000606 for the collapsed route. Their maximum mutual difference was
0.000444. On observed outcomes, their maximum PIP difference was 0.00967 and
maximum q-mean difference was 0.00522. The beta(1,9) sparse prior reduced the
posterior expected number selected relative to beta(1,1), as expected; this is
prior sensitivity, not a correctness criterion. **SBS3A-R1 passed.**

### Phase 3B: stick-specific slab-variance hierarchy

The R reference uses the explicit inverse-gamma convention

\[
p(\tau_k^2)=\frac{b_\tau^{a_\tau}}{\Gamma(a_\tau)}
(\tau_k^2)^{-a_\tau-1}\exp(-b_\tau/\tau_k^2),
\]

written `IG(shape, scale)`. Conditional on the selected coefficients,

\[
\tau_k^2\mid\alpha,\delta\sim\operatorname{IG}
\left(a_\tau+\frac M2,
b_\tau+\frac12\sum_{j:\delta_j=1}\alpha_{jk}^2\right).
\]

Excluded zero coefficients do not enter this conditional. When M=0, the
conditional is exactly the prior. The qualification prior is IG(3,1.6), whose
mean is 0.8 and whose variance is finite. A more concentrated IG(10,7.2)
prior, also with mean 0.8, provides the declared sensitivity comparison.

The schedule is z update, joint selection/regeneration conditional on the
current tau values, blocked coefficient redraw, then tau updates conditional
on the selected coefficients. Thus no coefficient is integrated under one
tau value and subsequently used as if it had been integrated under another.

For J=1 and K=1, deterministic quadrature in log tau integrates the exact
flat-intercept Gaussian model likelihood against the inverse-gamma prior.
The quadrature PIP was 0.999512 versus 0.999526 from MCMC; the marginal
posterior tau mean was 0.65754 versus 0.66316. An empty model produced a mean
0.79849 under the prior mean 0.8. On the observed-d fixture, posterior tau
means were 0.661, 0.577, and 0.625, with R-hats below 1.00001 and ESS above
20,000. **SBS3B-R1 passed.**

Across ten additional observed-d datasets with fixed pi_A and learned tau,
the mean posterior tau was 0.657. Because these fixtures use fixed generating
slopes rather than draws from the slab, this is empirical regularization, not
tau truth recovery. Mean informative and null PIPs were 0.802 and 0.0766, and
the largest tau R-hat was 1.00284.

### Phase 3C: joint hierarchy and Gibbs schedule

The full reference target is

\[
p(z,\delta,\alpha,\alpha_0,\pi_A,\tau^2\mid d,A).
\]

One explicit iteration is:

1. draw z given observed continuation outcomes and current coefficients;
2. draw each joint (\delta_j,\alpha_{j1:K}) block using current (\pi_A)
   and (\tau^2);
3. jointly redraw each stick intercept and all selected slopes;
4. draw (\pi_A) given the completed selection state;
5. draw each (\tau_k^2) given the completed selected coefficient state.

The independent collapsed route replaces step 2's sampled-pi prior odds with
the beta-binomial odds and regenerates (\pi_A) only after the selection sweep.
Both are compositions of exact Gibbs or partially collapsed/regenerated Gibbs
blocks and target the same joint posterior.

### Sparsity/scale dependence and diagnostics

The hierarchy can trade fewer, larger annotation effects against more,
smaller effects. Therefore the reference retains pi_A, tau by stick, M,
annotation states, and coefficients, and reports chain-specific PIPs, switch
counts, R-hat/ESS, and correlations among pi_A, tau, and M.

On the three-annotation primary fixture, explicit and collapsed routes differed
by at most 0.00340 in PIP, 0.00139 in tau means, and 0.00330 in q means.
Primary PIPs were 0.99696, 0.97028, and 0.01236. The posterior mean pi_A was
0.59696 with 95% interval [0.185,0.937]; posterior tau means by stick were
0.659, 0.574, and 0.623. All pi_A, tau, and M R-hats were within 0.00008 of
one. Correlation(pi_A,M) was 0.203; absolute pi_A/tau and M/tau correlations
were at most 0.058 on this fixture. These small-fixture values are diagnostic,
not universal claims about the model.

### Moderate-J and correlated-annotation reference

A second fixture uses J=12, N=260, three informative annotations, eight null
annotations, and a strong proxy for one informative continuous annotation.
With beta(1,9), posterior mean pi_A was 0.174, posterior mean M was 2.82, and
tau means were 0.547, 0.583, and 0.611. Mean PIP for the three generating
signals exceeded the mean across the null annotations.

For the signal/proxy pair, posterior probabilities were 0.6015 for signal
only, 0.3709 for proxy only, 0.0276 for both, and zero (at retained precision)
for neither. This deliberately records substitutable correlated evidence; the
reference does not force an arbitrary unique annotation.

Pairwise inclusion probabilities and posterior delta correlations are
reference summaries only. Grouped selection remains deferred.

### Bayesian FDR reference summary

For threshold t and selected set (S_t={j:PIP_j\ge t}), the development-only
summary is

\[
\operatorname{BFDR}(t)=|S_t|^{-1}\sum_{j\in S_t}(1-PIP_j),
\]

when the selected set is nonempty. On the moderate fixture, t=0.8 selected one
annotation with BFDR 0.0214. No public package API is added.

### Prior predictive and structural guards

For pi_A~Beta(a_pi,b_pi), the implementation checks

\[
E(\pi_A)=\frac{a_\pi}{a_\pi+b_\pi},\qquad
E(M)=J\frac{a_\pi}{a_\pi+b_\pi}.
\]

Duplicate annotation columns remain exchange-symmetric and annotation-column
permutations preserve the posterior up to the same permutation. With all
annotations excluded, slopes remain exactly zero, pi_A follows its beta
conditional based on M=0, tau receives no likelihood contribution and reverts
to its prior, and the component probabilities normalize.

When pi_A is learned, an exactly zero annotation still has likelihood Bayes
factor one, but its marginal PIP need not equal a fixed pi_A value. Its odds
are determined entirely by the posterior sparsity hierarchy induced by the
other annotations. This replaces, rather than contradicts, the fixed-pi
Phase-1/2 zero-column identity.

Conditioning all deltas to one yields the corresponding continuous-alpha
hierarchical probit model with learned tau. This is a model bridge; standard
SBayesRC production code and its variance-prior parameterization remain
unchanged.

### Repeated hierarchical simulation

Twenty observed-d datasets used J=10, N=140, four deliberately initialized
short chains, pi_true=0.25, tau_true=0.8, and coefficients generated from the
stated hierarchy. Mean informative and null PIPs were 0.8836 and 0.0462.
Across replicate summaries, mean posterior pi_A was 0.258 and mean posterior
tau was 0.853. One of 20 deliberately short fits exceeded the preregistered
PIP-range/tau-R-hat diagnostic; the allowed gate was four. Small J implies
broad hyperparameter uncertainty, so these results establish coherent
calibration and exploration rather than high-precision recovery.

Across genuinely included annotation-stick effects, mean conditional-alpha
bias was -0.0100 and empirical 95% interval coverage was 0.972. The broad
95% posterior intervals covered pi_true in all 20 datasets and covered each
stick's tau_true in aggregate for all 20 datasets. These perfect small-sample
hyperparameter coverage proportions reflect deliberately broad posteriors and
must not be read as high-precision estimation.

### Phase-3 decision and remaining boundary

All Phase-1 and Phase-2 tests remain required. Phase 3A and 3B exact gates,
the joint explicit/collapsed comparison, primary and moderate-J diagnostics,
prior predictive checks, correlated-annotation summaries, and repeated
hierarchical simulation passed.

**SBS3-R1: the full shared-delta, global-pi_A, stick-specific-tau SBayesRC-S
hierarchy is validated in standalone R.**

The next permissible task is Phase 4 design of a separate production C++
SBayesRC-S backend using this R hierarchy as its oracle. This document does
not create such a backend, a public model identifier, Rcpp exports, or Study
07. Standard SBayesRC remains a separate validated posterior model.

## Phase-4 implementation finding: empty genomic sticks

The isolated C++ hierarchy passed Phase-4A R parity (`SBS4A-R1`). During the
first selection-enabled CSR genomic screen, however, the allocation sampler
entered a state with no SNP eligible for a later stick. Under the validated
flat intercept convention, that stick has zero likelihood precision and zero
prior precision for its intercept, hence an improper conditional. The native
implementation correctly refused this state rather than inventing a
numerical fallback.

Consequently Phase 4B is blocked at `SBS4B-R3` and the overall implementation
decision is `SBS4-R4`. The fixed-z and observed-d reference results remain
valid on their nonempty eligible sets, but the Phase-3 posterior is not yet a
proper genomic model over the complete allocation state space. A future
mathematical revision must specify and independently validate a proper
intercept prior or another genuinely proper empty-stick model before genomic
implementation resumes.
