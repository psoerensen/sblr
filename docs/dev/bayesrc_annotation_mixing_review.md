# BayesRC annotation mixing review

## Scope and decision

This review starts at `e7eba992ce3b3ce22e684d6fd96d0967c069cb30`, after the
proper probit-stick intercept prior was committed. It uses Study 06 replicate 1
only as read-only package-development evidence. No qualification capsule,
threshold, prior, annotation, seed, operator policy, or benchmark result was
changed.

The decision is **Outcome D: component allocation/occupancy movement is the
demonstrated bottleneck**. A blocked alpha update is not justified by the
available evidence, and a non-centred parameterization is not justified by the
observed alpha--variance dependence. A separate posterior-preserving allocation
move would require its own derivation and validation task.

## Conditional audit and correction

For a stick with Albert--Chib latent response `z`, the committed conditionals
remain:

```text
intercept precision = x0'x0 + 1 / tau0^2
intercept rhs       = x0' residualized-z + mu0 / tau0^2
annotation precision = xk'xk + 1 / sigmaSqAlpha
annotation rhs       = xk' residualized-z
sigmaSqAlpha = (sum_k alpha_k^2 + b) / chi-square(q + a)
```

The intercept is excluded from the `sigmaSqAlpha` hierarchy and `q` is the
number of non-intercept coefficients. Empty sticks draw the intercept and every
non-intercept coefficient from their conditional priors before the unchanged
variance update. Proper-prior separated non-empty sticks use the ordinary
latent-response likelihood.

One defect was found. For a non-empty eligible set, a non-intercept annotation
column can have zero norm on that subset. The coefficient-wise loop previously
skipped it and retained its old value. With no likelihood contribution its
exact full conditional is `Normal(0, sigmaSqAlpha)`. The shared kernel now makes
that draw. This removes an artificial persistent direction and prevents the
stale coefficient entering the following variance update. No prior, likelihood,
RNG stream construction, component definition, or operator contract changed.

## Unchanged-kernel Stage A evidence

Both BED runs used 1,000 unthinned iterations in each of the four original
Study 06 chains. The informative run took 1,172 seconds. All annotation and
variance states were finite; alpha ranged from -8.56 to 8.45 and
`sigmaSqAlpha` from 0.103 to 192.1. Chain locations still differed, including
the component-0 intercept means (-3.28 to -2.79) and sparsely occupied later
sticks.

The uninformative run used a 200-iteration diagnostic burn-in, retained 800
draws per chain, and took 860 seconds. Its worst registered results were:

| Quantity family | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS | Maximum relative MCSE |
|---|---:|---:|---:|---:|
| Component-0 alpha | 2.122 | 5.22 | 11.33 | 0.451 |
| Implied prior summaries | 1.870 | 5.72 | 19.10 | 0.429 |
| `sigmaSqAlpha` | 1.026 | 115.85 | 316.05 | 0.064 |
| Variance components | 1.033 | 118.64 | 636.68 | 0.091 |

The slowest alpha effective rate was about 0.36 bulk effective draws per
minute. Component-1 `sigmaSqAlpha` achieved about 8.1 bulk effective draws per
minute, while expected active count achieved about 0.43. Component-0 alpha
lag-50 autocorrelation reached 0.81.

The scientifically relevant implied probabilities did **not** converge while
raw alpha failed: expected active count and enriched/unannotated prior summaries
shared the same slow movement. Thus raw-parameter diagnostics cannot be removed
from the qualification contract on an invariance argument.

## Geometry and occupancy

The largest observed absolute alpha--`sigmaSqAlpha` correlation was about
0.43. Late-stick variance chains were broad but crossed scale regimes, with
R-hat 1.020--1.026. This is not strong evidence that a centred hierarchical
funnel is the primary bottleneck.

Final uninformative BED component occupancy differed substantially across
chains. Component 1 contained 51, 14, 11, and 9 markers; component 2 contained
3, 6, 1, and 2 markers. Component-0 alpha and implied active probabilities had
the strongest long-lag persistence. These facts support slow marker-allocation
movement, coupled to early-stick annotation probabilities, rather than a
conditional Gaussian coefficient-update bottleneck. Stages B and C were not
run because Stage A already rejects the premise needed to justify an alpha-only
sampler change.

## Retained block eigen

An informative replicate-1, one-chain retained-0.995 development control ran
for 1,000 iterations with the proper prior. Residual variance remained
0.769--1.147, genetic variance remained positive, and 1,001 projected-residual
rebuilds had maximum absolute drift `2.66e-15`. This confirms numerical operator
stability, but it is not a four-chain mixing result and supplies no
qualification evidence.

## Diagnostics and tooling

`tools/pilot_study06_proper_intercept_prior.R` now accepts configurable burn-in
and output locations and writes compact convergence, per-chain autocorrelation,
cross-parameter correlation, and final-occupancy CSV summaries outside the
benchmark repository. `tools/diagnose_study06_sbayesrc_operators.R` provides the
retained-eigen stability control and now reports its actual diagnostic history
length.

Production fits expose true alpha and `sigmaSqAlpha` convergence traces, but do
not yet aggregate changing eligible/continuation counts and occupancy ranges by
chain and stick. Adding those compact summaries should accompany a future
allocation-move task; storing iteration-level classification matrices by
default is not recommended.

## Limitations and recommendation

The informative Stage A run predates the compact diagnostic extension, and the
retained control is one chain. No 3,000- or 9,000-iteration current-posterior
stage was run because the 1,000-iteration evidence already identifies a
different bottleneck than coefficient-wise alpha Gibbs. Study 06 must not be
declared qualified.

The next `sblrbench` action is **none yet**. First derive and test a
posterior-preserving marker-component allocation move in `sblr`, or decide in a
separate method review that qualification-length chains should assess the
existing allocation kernel. Only after that package decision should
`sblrbench` pin a new SHA and perform its separate design-update review.
