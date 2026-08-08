# MCEM-SBayesRC Phase 5A

## Status and inference target

Phase 5A added an internal, single-trait CSR qualification path for what is now
named SBayesRC-EM. SBayesRC-EM estimates continuous
annotation coefficients at an observed-data MAP / empirical-Bayes solution
using a Monte Carlo E-step. It is not a sampler from the full joint alpha
posterior.

The existing joint SBayesRC, SBayesRC-S, and BayesR samplers are unchanged.
SBayesRC-S annotation selection is outside the scope of Phase 5A.

## Architecture

The internal `.stblr_mcem_sbayesrc_csr()` driver keeps the outer MCEM loop in
R and reuses the existing CSR BayesR/SBayesRC marker kernel for each genomic
block. The separation is:

```text
R outer MCEM driver
  -> construct SNP-specific continuation-stick priors from alpha
  -> existing C++ CSR genomic block with those priors fixed
  -> stream and average the existing conditional component probabilities
  -> R penalized soft-probit M-step
  -> damp alpha and test convergence
```

The C++ information-flow accumulator already normalizes the same `logp`
vector used by the unchanged categorical component draw. Normalization occurs
after the hard component draw, consumes no random numbers, and stores only an
`M x K` running sum. Phase 5A reuses its post-burn `rb_comp_prob` average; it
does not retain an iteration-by-marker-by-component cube.

The R/C++ boundary is intentional. Genomic work remains in the optimized
native kernel, while the inexpensive M-step and outer history stay directly
comparable with the gold-standard R implementation.

## M-step and convergence

Responsibilities are converted once per outer iteration to soft continuation
eligibility and success weights. Every SNP contributes one expected latent
observation; retained inner sweeps are averaged and never treated as repeated
annotation observations.

Each stick is fitted independently by a BFGS penalized probit optimization.
The intercept and non-intercept Gaussian priors use the established proper
SBayesRC prior inputs. Fixed stick-specific `sigmaSqAlpha` values provide the
non-intercept prior variance. The update is damped as

```text
alpha_next = (1 - damping) * alpha + damping * alpha_mstep
```

with Phase-5A defaults `damping = 0.5`, `tol_alpha = 1e-3`,
`tol_prior = 1e-3`, `min_outer = 3`, and `max_outer = 50`. Convergence requires
both the maximum alpha change and maximum SNP component-prior change to pass
their tolerances. The complete alpha and aggregate diagnostic history is
retained even if `max_outer` is reached.

## Internal interface and result

No public model identifier, export, or `NAMESPACE` entry is added. The
development entry point is `.stblr_mcem_sbayesrc_csr()`. It returns:

- `genomic`: the ordinary raw final CSR genomic block, sampled with the learned
  SNP priors frozen;
- `mcem$method`: `"SBayesRC-EM"` and `mcem$algorithm`: `"MCEM"`;
- `mcem$target`: the observed-data alpha MAP / empirical-Bayes target;
- `mcem$alpha_map` and `mcem$component_prior`;
- convergence state, iteration count, tolerances, damping, and block lengths;
- aggregate outer history and alpha history;
- optional final/E-step RB responsibilities for development validation.

The MAP is deliberately not written into joint-SBayesRC posterior-alpha
fields. The final genomic block is separate from the E-step blocks and retains
the existing raw genomic output semantics.

## Phase-5A restrictions

The qualified path is deliberately narrow:

- continuous annotations only;
- one trait and the CSR summary-statistics route;
- fixed genomic effect and residual variance parameters (`adjE = 0`);
- fixed `sigmaSqAlpha` during MCEM;
- no global-mixture update overriding the SNP-specific annotation priors;
- no annotation selection, `delta`, `pi_A`, or selection `tau2` hierarchy.

Releasing genomic hyperparameters is deferred and must be tested one at a
time. The internal interface rejects configurations outside the qualified
fixed-hyperparameter regime instead of silently enabling them.

## Qualification evidence

Permanent tests cover the annotation-to-component probability map, the soft
stick transformation, independent M-step parity from multiple starts, and RNG
neutrality of RB accumulation. They also run compact orthogonal and
correlated-LD end-to-end fixtures.

The development qualification tool
`tools/sbayesrc_mcem_qualification.R` reruns the full local gold-standard cases.
Against their independently optimized exact targets it produced:

| case | outer iterations | converged | alpha RMSE | max alpha error | SNP-prior RMSE | RB responsibility RMSE |
|---|---:|:---:|---:|---:|---:|---:|
| orthogonal | 28 | yes | 0.002277 | 0.004428 | 0.000253 | 1.05e-14 |
| correlated LD | 17 | yes | 0.002920 | 0.004698 | 0.000789 | 0.000816 |

These results reproduce the validated MCEM-R1 regime using the production CSR
genomic machinery. They qualify implementation parity, not broad genomic
hyperparameter learning or production readiness.

## Files and remaining work

The implementation lives in `R/sbayesrc-mcem.R`; compact independent guards
live in `tests/testthat/helper-sbayesrc-mcem-reference.R` and
`tests/testthat/test-sbayesrc-mcem.R`. No native sampler source was changed.

The next phase should release currently fixed genomic hyperparameters one at a
time and re-establish reference behavior before any MCEM extension to
SBayesRC-S annotation selection.

Phase 5B completed that gated block-eigen and genomic-hyperparameter work; see
[`sbayesrc_em_phase5b.md`](sbayesrc_em_phase5b.md). This document remains the
CSR reference qualification record.
