# SBayesRC-EM Phase 5B implementation

## Identity and status

SBayesRC and SBayesRC-EM share the continuous continuation-stick annotation
model but use different inference targets:

```text
SBayesRC
  model = sbayesrc
  inference = joint
  full joint Bayesian alpha inference

SBayesRC-EM
  model = sbayesrc
  inference = mcem
  observed-data alpha MAP / empirical-Bayes inference
```

The method name is **SBayesRC-EM** and its fitting algorithm is **MCEM**.
SBayesRC-EM is not a modified sampler for the joint alpha posterior.
SBayesRC-S and annotation selection remain outside Phase 5B.

## Architecture

The backend-neutral R outer engine in `R/sbayesrc-mcem.R` owns:

```text
alpha -> SNP priors -> genomic E-step block -> RB responsibilities
      -> soft-probit M-step -> damping -> convergence
```

Two genomic backends implement the fixed-prior inner block:

- CSR is retained as the transparent exact-reference and debugging backend;
- retained low-rank block-eigen is the intended scalable production backend.

Both reuse the existing SBayesRC marker/allocation kernel. The block binding
enables the existing information-flow accumulator through its internal block
configuration object, avoiding any new native-call argument or public API.
The same normalized `logp` used by the unchanged hard component draw is
accumulated after the draw, so RB capture consumes no RNG and stores only one
running marker-by-component matrix.

## Block residual contract

Block SBayesRC-EM preserves the established retained-block contract:

```text
residual_policy = gctb_block
block_ve_mode = fixVe       when E is fixed
block_ve_mode = allMixVe    when E is learned
```

No residual mathematics or default of the ordinary block sampler changed.
For fixed-parameter CSR/block comparisons, CSR residual variance is aligned to
the fixed block phenotype-variance contract; it is not compared under a
different residual model.

## Hyperparameter semantics

### Genomic residual and effect variance

The production block path can use the existing genomic updates for E, B, or
both while alpha and SNP-specific priors remain fixed within each E-step.
Their final states are carried between outer blocks, and their values are
recorded in every outer-history row. This changes the genomic conditional
integrated by the E-step, not the annotation MAP objective.

### Alpha prior variance

Phase 5B finalizes `sigmaSqAlpha` as a **fixed prior variance**. The M-step
maximizes the validated objective with a fixed Gaussian alpha prior. Sampling
`sigmaSqAlpha` using the joint-SBayesRC Gibbs update would mix joint-posterior
and MAP semantics. A future empirical-Bayes outer update would require a new,
explicitly derived objective and is not silently enabled here.

### Global mixture probability

There is no independent global `Pi` update. SNP-specific probabilities are
defined completely by `A %*% alpha`; continuation-stick intercepts provide the
baseline mixture propensity. The validated R MCEM reference fixes this same
parameterization. Adding a separate global mixture update would override or
double-count the annotation-derived prior, so Phase 5B records:

```text
mixture_prior_mode = annotation_stick_intercepts_no_global_Pi_update
```

## Final genomic block and output

After convergence or `max_outer`, the engine runs a distinct genomic block
with the learned component priors frozen and alpha updates disabled. Its raw
genomic summaries retain ordinary backend semantics and are conditional on the
estimated annotation priors. The last E-step and final-block responsibilities
are distinguished as:

```text
mcem$last_estep_responsibilities
mcem$final_genomic_responsibilities
```

The internal result records `method = "SBayesRC-EM"`, `algorithm = "MCEM"`,
`inference = "mcem"`, backend, convergence controls, alpha MAP, component
priors, B/E modes and final values, and compact outer histories. Phase-5A
responsibility aliases remain temporarily available for internal compatibility.
No public model switch or export is introduced in Phase 5B.

## Gate results

The consolidated qualification tool is
`tools/sbayesrc_mcem_qualification.R`.

### Gate 1: fixed-parameter block parity

PASS. With residual contracts aligned, the full qualification observed maximum
CSR/block differences of 0.00127 for alpha, 0.00070 for SNP priors, and 0.00113
for RB responsibilities. Dispersed block starts differed by at most 0.0391 in
alpha. RB diagnostics on/off produced identical scientific state.

### Gate 2: learned E

PASS. The existing `allMixVe` update remained finite, RNG-neutral under RB
capture, and start-robust. Relative to the fixed-E fit, maximum differences
were 0.00473 for alpha, 0.00258 for priors, and 0.00643 for responsibilities.

### Gate 3: learned B and learned B+E

PASS. Both configurations remained finite and start-robust. Maximum alpha
start differences were 0.0450. Relative to the fixed-B/E reference, maximum
alpha differences were 0.0642 for learned B and 0.0605 for learned B+E.

The inner-length screen found alpha/prior solutions stable from 250 through
1,000 retained sweeps on the qualification fixture. B itself was noisy at 250
sweeps, while 500 and 1,000 were close; shorter genomic blocks are therefore a
future performance study rather than a changed default.

### Gates 4–6

PASS. Fixed alpha-prior variance and no-global-Pi semantics are explicit and
tested. The final conditional genomic block is separate, freezes alpha/prior,
and returns ordinary genomic summaries plus explicitly named MCEM metadata.

## Current boundary and next work

The implementation remains internal and single-trait. CSR remains the
qualification backend; retained low-rank block-eigen is the intended scalable
backend. Broad architecture, runtime, and predictive validation belongs in
`sblrbench`, not package tests.

The next `sblr` model-development task may implement SBayesRC-S-EM using its
separately validated selection hierarchy. Study 07 should begin only after the
desired method set is complete.
