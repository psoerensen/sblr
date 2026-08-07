# Retained block SBayesRC mixing development result

## Decision

**BLOCK-MIX-R4 — block transition works but global alpha remains stuck.**

The probit sandwich transition is posterior-correct but insufficient
(`PX-B2`). The fixed-alpha conditional particle-Gibbs block reference is also
posterior-correct and does not suffer weight or ancestry collapse at 500
markers. It moves many component labels locally, but it produces only small
net active-count changes. This does not yet address the demonstrated global
alpha/occupancy regime separation. The particle kernel therefore remains a
development reference and is not integrated into the production sampler.

## Provenance and scope

- source HEAD: `740b837ed5f794dc39d57f0706f5a47f718f76e4`;
- package: `sblr` 0.2.0;
- R 4.4.1, GCC/G++ 13.2, C++17, OpenMP enabled;
- pinned standalone SBayesRC: v0.2.6,
  `b95d3fcbad8ff358290922a58fff893439296138`;
- pinned GCTB: `cc7fa7d765c83a89c6375946cf77fe50ba1a317e`;
- read-only `sblrbench` HEAD:
  `fbe80603ff6fa09e0a611a56d09130cb4b2cbc8c`.

The scientific screen uses only retained block eigen with
`residual_policy = "gctb_block"` and `block_ve_mode = "allMixVe"`. BED and CSR
were not run or changed. No formal Study 06 qualification or 38k-marker
scientific fit was launched.

## GCTB/SBayesRC audit

Standalone SBayesRC and GCTB ApproxBayesRC both use stick-wise
Albert--Chib updates followed by ordinary SNP updates. Neither contains a
hidden global or LD-block alpha/allocation move. GCTB's named joint annotation
update belongs to ApproxBayesRD and adds a BayesC point-mass annotation
inclusion state shared across sticks. That is a different model and was not
ported. See [the method review](sbayesrc_block_mixing_method_review.md).

## PX exact validation

The implemented development option scales the sign-constrained latent vector
on a positive orbit and uses the exact integrated-latent MH ratio, followed by
a blocked Gaussian alpha draw and the unchanged `sigmaSqAlpha` conditional.
The direct and reduced log ratios agreed in deterministic tests. The existing
8-marker AA oracle used four chains, 30,000 iterations, 5,000 burn-in and
thinning five. Prespecified gates all passed:

| quantity | ordinary | PX | tolerance |
|---|---:|---:|---:|
| maximum absolute PIP error | 0.00452 | 0.00289 | 0.010 |
| maximum absolute beta-mean error | 0.00159 | 0.00149 | 0.015 |
| active-count distribution L1 error | 0.01899 | 0.01384 | 0.030 |
| maximum absolute alpha-mean error | 0.00515 | 0.01203 | 0.040 |
| PX scale acceptance | — | 0.53065 | descriptive |
| mean alpha jump | 0.69255 | 0.89748 | descriptive |

Ordinary mode remains byte-identical under fixed seeds after excluding wall
clock construction/timing diagnostics. Enabling PX diagnostics does not consume
sampler RNG.

## Matched 1,500-marker short screen

The safety-gated screen used the exact informative Study 06 v2 identity,
four registered chain seeds, 1,000 recorded iterations, a 300-iteration
burn-in and 700 descriptive retained draws per chain. It is not a
qualification-length result.

| diagnostic | ordinary | PX |
|---|---:|---:|
| runtime, seconds | 133.08 | 127.38 |
| passing quantities / 35 | 1 | 1 |
| maximum R-hat | 2.106 | 1.543 |
| minimum bulk ESS | 5.38 | 7.44 |
| minimum tail ESS | 14.03 | 9.66 |
| maximum relative MCSE | 0.476 | 0.380 |
| realized-active R-hat | 2.091 | 1.360 |
| realized-active bulk ESS | 5.43 | 9.10 |
| expected-active R-hat | 2.103 | 1.370 |
| expected-active bulk ESS | 5.42 | 8.99 |

PX brought three chains into similar active-count ranges, but one chain
remained higher. First-stick intercept R-hat changed from 1.452 to 1.543 and
first-stick continuous-alpha R-hat from 1.334 to 1.506. Thus the occupancy
improvement did not constitute joint hierarchy convergence.

Compact stick traces tell the same story. First-stick continuation R-hat
improved from 2.091 to 1.360 (bulk ESS 5.43 to 9.10), while second-stick stop
R-hat was 1.231 and third-stick stop R-hat was 1.401 under PX. The constant
first-stick eligible count is structurally 1,500 and has no defined R-hat.
Stick-specific `sigmaSqAlpha` R-hats changed from 1.289/1.038/1.054 to
1.044/1.117/1.091: PX helped the first variance but worsened later-stick
exploration at this horizon.

Scientific safeguards were stable at this short horizon: PIP AUPRC was 0.5582
ordinary versus 0.5577 PX; effect/truth correlation 0.8857 versus 0.8858;
validation genetic-value correlation 0.8961 versus 0.8950; and phenotype
prediction correlation 0.6323 versus 0.6314.

## Conditional particle-Gibbs validation

For fixed alpha and block variance, the development target is

```text
exp{-||w - Q beta||^2 / (2 Ve)}
prod_j pi[j,c[j]] p(beta[j] | c[j], vb, gamma).
```

The prefix target treats future effects as zero. Its fully adapted
component/effect proposal has the ordinary integrated marker conditional, so
the incremental importance weight is its normalizer. Conditional SMC fixes
the current path as particle one and selects a terminal path by its exact
weight. A 3-marker two-component test with 16 particles and 10,000 retained
updates matched exact allocation enumeration within a prespecified 0.025 PIP
tolerance. Diagnostics were RNG-neutral.

The state-based scale screen used the preserved large B0 block factor, not a
new fit. The factor was 500 by 500; the 100-marker test conditionally removed
the other 400 effects. Results average five refreshes at 100 markers and two at
500 markers:

| markers | particles | min weight ESS | median weight ESS | changed allocations | active-count jump | seconds |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 8 | 8.00 | 8.00 | 1.2 | 0.0 | 0.20 |
| 100 | 16 | 16.00 | 16.00 | 5.4 | -0.2 | 0.17 |
| 100 | 32 | 32.00 | 32.00 | 2.4 | 0.0 | 0.34 |
| 100 | 64 | 64.00 | 64.00 | 4.8 | 0.0 | 0.69 |
| 500 | 8 | 7.99 | 8.00 | 14.0 | -3.5 | 0.45 |
| 500 | 16 | 15.99 | 16.00 | 26.0 | -2.5 | 0.83 |
| 500 | 32 | 31.98 | 32.00 | 31.5 | -1.5 | 1.66 |
| 500 | 64 | 63.96 | 63.99 | 30.5 | -0.5 | 3.09 |

No resampling was triggered and terminal component/effect paths remained
diverse. This rules out ordinary particle-weight collapse for the tested
fixed-alpha state. It also shows the limitation: component changes largely
balance entries and exits, leaving net block sparsity nearly unchanged. A
local block refresh alone has no demonstrated mechanism for moving the global
alpha/expected-active regime.

## Implementation boundary and next task

The exact PX kernel is available only through dot-prefixed block development
controls and is disabled by default. The particle kernel remains an internal R
reference. No public signature, default, output schema, likelihood, prior,
residual policy, or ordinary RNG path changed.

The next method task should derive a global alpha transition with an exact
block-factorized marginal likelihood estimator, and first measure estimator
variance using common random numbers/correlated particle seeds. Production
particle integration is not justified until that global transition can make
meaningful expected-active moves at acceptable estimator variance. The GCTB
annotation-selection hierarchy remains a separate alternative-model question.
