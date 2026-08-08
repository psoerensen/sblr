# SBayesRC-S-EM Phase 5C development record

## Method identity

SBayesRC-S-EM is the annotation-selection counterpart of SBayesRC-EM. It is
separate from joint SBayesRC-S:

- joint SBayesRC-S samples the full posterior of `delta`, `alpha`, `pi_A`,
  `tau2`, and the genomic state;
- SBayesRC-S-EM estimates a point `(delta, alpha)` annotation model by MCEM and
  performs final genomic inference conditional on that fitted model.

The Phase-5C algorithm is named **MCEM-Laplace**. MCEM describes the genomic
E-step and discrete/continuous MAP M-step. Laplace model-space probabilities
provide the explicitly conditional annotation inclusion summary described
below.

## Gate 1: statistical target

For Phase 5C, `pi_A` and stick-specific `tau2` are fixed hyperparameters. For a
shared nonintercept selection vector `delta`, the observed-data MAP target is

```text
log p(y | delta, alpha) + log p(delta | pi_A) +
sum_k log p(alpha_k | delta, tau2[k]).
```

The proper intercept prior is present for every stick and every allocation
state. The genomic allocations, marker effects, and other genomic latent state
are integrated by the Monte Carlo E-step. The M-step maximizes

```text
Q(delta, alpha | current) =
sum_i,k r[i,k] log component_prior[i,k](delta, alpha) +
log p(delta | pi_A) + log p(alpha | delta, tau2),
```

where `r` is one averaged Rao-Blackwell responsibility vector per SNP. A
single-bit model-space search plus exact continuous optimization within each
visited model avoids exponential enumeration in the fitted implementation.
Exact `2^J` enumeration is qualification-only.

`delta` and `alpha` are point-estimated. `pi_A` and `tau2` are neither sampled
nor optimized in Phase 5C. Reusing their joint-SBayesRC-S Gibbs updates would
mix a full-posterior hierarchy into a MAP inference line and is therefore not
coherent here.

## Meaning of `annotation_pip_eb`

At the converged E-step responsibilities, each annotation model is assigned a
proper Gaussian-prior marginal weight for the expected complete-data probit
target. The continuous coefficients are integrated with a Laplace
approximation and the shared-`delta` model probabilities are explored with a
single-bit model-space Gibbs chain. The reported quantity is

```text
annotation_pip_eb[j] =
  P_Q(delta[j] = 1 | final RB responsibilities, fixed pi_A, fixed tau2),
```

under this responsibility-conditioned Laplace target. It is a genuine
conditional inclusion probability for the defined outer annotation target,
but it is an approximation to an observed-data inclusion probability: genomic
uncertainty enters through the final averaged responsibilities, rather than by
jointly integrating the complete genomic and selection hierarchy.

`annotation_pip_eb` is not identical to the joint SBayesRC-S posterior
inclusion probability. The latter integrates posterior uncertainty in
`pi_A`, `tau2`, annotation coefficients, latent probit values, allocations,
and genomic parameters. Phase 5C fixes `pi_A` and `tau2`, uses a MAP annotation
model for the next genomic E-step, and uses the Laplace model distribution only
for conditional inclusion summaries and model-averaged coefficient summaries.

No posterior intervals are inferred from outer iterations. They are
optimization diagnostics, not posterior samples.

## Alternatives considered

1. **Fixed `pi_A` and `tau2` (selected).** This gives a clean MAP objective,
   preserves the validated slab/intercept definitions, and permits an
   independently checkable conditional inclusion target.
2. **Empirical-Bayes `pi_A` or `tau2`.** These require additional marginal or
   ECM derivations. They are not necessary for the first qualified method and
   are deferred rather than approximated silently.
3. **Joint Gibbs updates (rejected).** Sampling `delta`, `pi_A`, and `tau2`
   inside the outer loop would reproduce pieces of joint SBayesRC-S rather than
   define an MCEM/MAP method.

## Architecture and qualification gates

The genomic block is unchanged and reused from SBayesRC-EM:

```text
MAP delta/alpha
  -> SNP-specific continuation-stick priors
  -> fixed-prior CSR or block-eigen genomic block
  -> averaged Rao-Blackwell responsibilities
  -> shared-delta MCEM-Laplace annotation update
  -> repeat
  -> separate final genomic block under frozen fitted priors
```

CSR remains the transparent reference backend. The retained block-eigen route
remains the scalable backend and preserves the established `gctb_block`, B,
and E contracts. The implementation does not change BayesR, joint SBayesRC,
SBayesRC-EM, or joint SBayesRC-S transitions.

The permanent qualification covers independent quadrature for tiny model
spaces, null/signal/correlated annotations, proper empty-stick behavior,
shared selection across sticks, dispersed starts, CSR/block parity, learned B
and E, final-block separation, and regression of every existing inference
line. Broad scientific calibration across architectures remains Study 07
work.

## Phase-5C gate results

All gates passed in the consolidated development qualification:

| Gate | Result |
|---|---|
| G1 | fixed-hyperparameter MAP target and conditional EB-PIP defined |
| G2 | independent quadrature versus Laplace/model-space maximum PIP error 0.0202 in the permanent fixture and 0.0169 in the consolidated run |
| G3 | CSR reference stable; dispersed-start PIP difference 0.010 |
| G4 | CSR/block PIP difference 0, alpha difference 0.00164, component-prior difference 0.000656, RB difference 0.00363 |
| G5 | learned E/B/B+E PIP start differences 0.0050, 0.00667, and 0.0117 |
| G6 | `pi_A_mode = fixed`; `tau2_mode = fixed` |
| G7 | controlled EB-PIP behavior qualified |
| G8 | final genomic block is separate and freezes the fitted annotation model |
| G9 | internal family naming finalized without a public API change |

The controlled behavior screen gave mean PIPs near 0.011 for all-null
annotations, 0.032 for a weak signal, 0.807 for a moderate signal, and 1.000
for a strong signal. In a correlated signal/proxy case the marginal PIPs were
0.953 and 0.083; with two independent informative annotations they were 0.991
and 0.792. Three independent model-space streams per case produced small
ranges. These are implementation qualification results, not a universal
calibration or threshold study.

The final internal result uses:

```text
method = "SBayesRC-S-EM"
algorithm = "MCEM-Laplace"
inference = "mcem"
```

It returns `delta_map`, `alpha_map`, `alpha_model_average`,
`annotation_pip_eb`, the final SNP component prior, complete outer histories,
fixed selection-hyperparameter metadata, and explicitly separate last-E-step
and final-genomic responsibilities. Joint SBayesRC-S retains its existing
`annotation_pip` semantics and output path.

## Remaining boundary

Phase 5C fixes `pi_A` and `tau2`; learning either requires a separate marginal
or ECM derivation. The Laplace approximation is qualified on small cases but
its accuracy with many highly collinear annotations remains a scientific
benchmark question. The path is single-trait and internal. Broad calibration,
architecture sensitivity, runtime, and memory studies belong in `sblrbench`
Study 07 rather than further method expansion in this package.

**Phase-5C decision: SBS-EM5C-R2 — SBayesRC-S-EM complete with fixed
selection hyperparameters.** The implementation/model-development set is
complete enough to proceed to broad `sblrbench` Study 07 evaluation; this is
not itself broad scientific production validation.
