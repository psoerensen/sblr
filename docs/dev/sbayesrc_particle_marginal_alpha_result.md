# Particle-marginal global alpha feasibility result

## Decision

The result is **PMA-R3: exact but computationally impractical**. The likelihood
estimator and selected-path PMMH construction passed the exact gates, and a
short 1,500-marker reference chain made genuinely coordinated alpha and
occupancy moves. It did not establish joint convergence, and measured cost on
representative 38k-marker blocks is too high to justify production integration
or a full scientific run.

No package default, production sampler, likelihood, prior, residual policy, or
public interface changed. The recommended next scientific task is Bayesian
annotation selection/annotation PIPs rather than further unrestricted
continuous-alpha sampler engineering.

## Provenance and immutable inputs

- `sblr`: `343ffc0cb085918e6276b8be5ba5ae0f368bf6ca`, package 0.2.0.
- `sblrbench` (read only): `fbe80603ff6fa09e0a611a56d09130cb4b2cbc8c`.
- Study 06 specification hash:
  `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`.
- Study 06 truth hash:
  `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.
- The 1,500-marker cache contains the same 15 full-rank retained blocks used by
  the committed BLOCK-MIX-R4 screen. The large-design audit used measured
  blocks 1, 20, 39, 58, and 76 of the 76-block frozen design.

## Exact target and estimator

For fixed non-alpha global state `theta`, the target is

```text
p(alpha | y, theta) proportional to
  p(alpha | sigmaSqAlpha) product_b L_b(alpha)
```

where `L_b` sums allocations and integrates effects in block `b`. The
fully-adapted marker-prefix proposal samples a component and its compatible
Gaussian effect. The product of normalized incremental factors is the exact
path target/proposal ratio. Averaging these weights over independent particles
gives `Zhat_b`, unbiased on the likelihood scale. Independent block random
variables make `product_b Zhat_b` unbiased.

The auxiliary state comprises standard-normal component, effect, and selected
path variables. A terminal path is selected using normalized path weights. On
acceptance, alpha, every block's selected allocations/effects, residuals, the
likelihood estimate, and auxiliary variables are replaced together; on
rejection all remain unchanged. The symmetric random-walk PMMH ratio is

```text
log r = log p(alpha') - log p(alpha)
      + log Zhat(alpha'; u') - log Zhat(alpha; u).
```

The alpha prior includes the proper stick-specific intercept normals and the
non-intercept normals conditional on `sigmaSqAlpha`. For correlated PMMH,
`u' = rho u + sqrt(1-rho^2) epsilon`; this Gaussian AR(1) kernel is reversible
with respect to the standard-normal auxiliary density.

## Tiny exact oracle

The five-marker, two-component oracle enumerated all 32 allocations and
integrated active Gaussian effects. Across 400 independent estimates:

| Particles | mean Zhat / exact Z | var(log Zhat) | final particle ESS |
|---:|---:|---:|---:|
| 4 | 0.99954 | 9.83e-5 | 4.00 |
| 8 | 1.00071 | 4.90e-5 | 8.00 |
| 16 | 0.99976 | 2.32e-5 | 15.99 |
| 32 | 1.00015 | 1.55e-5 | 31.99 |
| 64 | 1.00009 | 6.30e-6 | 63.97 |

Selected paths from the extended target agreed with exact PIPs within the
preregistered Monte Carlo tolerance. A discrete two-alpha PMMH oracle recovered
the exact alpha marginal. Fixed auxiliaries reproduced identical estimates and
the complete alpha-prior log ratio matched an independent calculation.

## Independent estimator scaling

For one retained large block, estimator noise was already small. At 500 markers
`sd(log Zhat_b)` ranged from 0.007--0.009 with eight particles to about 0.002
with 128 particles. Particle ESS remained essentially equal to particle count.

For the actual 15-block, 1,500-marker model:

| Particles | sd(log Zhat) | seconds / global estimate |
|---:|---:|---:|
| 8 | 0.536 | 0.65 |
| 16 | 0.235 | 0.85 |
| 32 | 0.245 | 1.09 |
| 64 | 0.188 | 0.99 |

The small non-monotonicity is compatible with 12 repetitions and timing noise;
the sum of measured block variances decreased from 0.239 to 0.024. Independent
PMMH noise was therefore not prohibitive at 1,500 markers.

## Correlated log-ratio screen

At 1,500 markers and 16 particles, `rho=0.99--0.999` reduced the measured SD of
the log-likelihood difference to approximately 0.03--0.18 for most registered
stick-wise and all-stick proposals. Alpha jumps of 0.05--0.25 changed expected
activity by up to about 16 markers. Both stick-wise and all-stick proposals had
non-negligible predicted acceptance for some preregistered directions. This
passed the narrow feasibility gate and authorized only the short reference
screen.

## Short selected-path PMMH screen

Four chains used 300 iterations, 100 burn-in iterations, eight particles per
block, `rho=0.99`, and a fixed stick-wise proposal SD of 0.10. The non-alpha
global state was held fixed, so this is a conditional PMMH reference, not a
replacement for the full SBayesRC sampler.

| Chain | Acceptance | Median alpha jump | Max alpha jump | Median absolute active jump | Max absolute active jump | Mean active |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.777 | 0.168 | 0.388 | 4 | 32 | 90.6 |
| 2 | 0.717 | 0.175 | 0.380 | 4 | 43 | 100.6 |
| 3 | 0.747 | 0.172 | 0.375 | 4 | 47 | 87.1 |
| 4 | 0.757 | 0.170 | 0.453 | 4 | 28 | 83.6 |

The move therefore solves the stale-allocation compatibility problem in the
narrow sense: accepted alpha moves install compatible full-genome particle
paths and can move realized occupancy by tens of markers.

It did **not** establish convergence. With only 200 retained draws, alpha R-hat
ranged from 1.27 to 3.34 and alpha ESS from 8.7 to 25.7. Realized active-count
R-hat was 1.14 (ESS 40.7); expected-active R-hat was 1.22 (ESS 32.2). Active and
expected-active lag-50 autocorrelations were near zero in most chains, but
between-chain locations still differed. `sigmaSqAlpha` was fixed in this
conditional reference, so no claim about its joint convergence is possible.
SNP-level safeguards were not promoted because the short latent chain did not
converge.

The four chains required 69--193 seconds each (547 seconds total). Timing
variation reflects the pure-R reference and host scheduling.

## Frozen 76-block feasibility

Five evenly spaced blocks were evaluated directly and their mean block
variance/cost projected to 76 blocks:

| Particles | projected sd(log Zhat), 76 blocks | projected seconds / likelihood |
|---:|---:|---:|
| 8 | 0.939 | 18.4 |
| 16 | 0.474 | 24.9 |
| 32 | 0.355 | 37.8 |
| 64 | 0.322 | 66.1 |

Thus likelihood noise is manageable, but reference cost is not. Registered
correlated proposal pairs projected to roughly 18--57 seconds each at 16
particles. Correlation `rho=0.99--0.999` reduced projected log-difference noise
to roughly 0.05--0.27 for most proposals, and some directions changed expected
activity by 100--370 markers. Direction-specific projected acceptance varied
from effectively zero to one, so those five-block extrapolations are overlap
diagnostics, not valid 76-block PMMH draws.

No worker-scaling experiment was run: the reference has no block-parallel
implementation, and production integration was barred once the cost and
convergence gates failed. Even ideal eight-way scaling would not establish the
required full joint convergence or SNP safeguards.

## Interpretation and limitations

The pseudo-marginal construction is not blocked by estimator variance or
selected-path exactness. It is blocked at the production gate by the combined
cost of repeated genome-wide particle systems and the absence of demonstrated
joint convergence in the short reference. A C++ implementation might be much
faster, but building it without evidence that a qualification-length chain is
scientifically and computationally defensible would violate the stopping rule.

The result does not diagnose a likelihood defect. It supports the narrower
conclusion that useful SNP prioritization can coexist with difficult full
posterior learning of an unrestricted continuous annotation hierarchy under
strong genome-wide alpha--allocation coupling.

## Recommendation

Do not add a production `particle_marginal` alpha kernel and do not launch a
full Study 06 run. Close the same-posterior sampler-development line at this
checkpoint. The next method-design task should specify Bayesian annotation
selection and annotation PIPs as a distinct model, with its own simulation,
estimands, priors, and qualification plan.

