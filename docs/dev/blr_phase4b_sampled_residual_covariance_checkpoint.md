# Phase 4b sampled residual-covariance checkpoint

Status: implemented, qualification-only, pending independent verification.

## Scope

Phase 4b extends the internal two-trait common-sample packed-BED Cheng
MT-BayesC$\Pi$ route from the
[Phase 4a checkpoint](blr_phase4a_cheng_mt_bayesc_checkpoint.md). It adds one
optional sampled full residual covariance without creating a second sampler.
The route remains internal: no public MT dispatch, legacy MT kernel, public
default, or `NAMESPACE` export is changed.

The supported residual policies are:

- `fixed_full`: the Phase 4a supplied finite symmetric positive-definite
  covariance, with no residual-covariance RNG draw;
- `sampled_full`: an explicitly initialized full covariance updated from an
  inverse-Wishart full conditional.

Both policies retain exactly two ordered traits, complete common-sample
phenotypes, one common packed-BED resource, the four Cheng activity patterns,
null collapse, conditional completion, sampled $\boldsymbol{\Pi}$, and the
authoritative $V_b$ transition.

## Prior and conditional

The sampled policy uses the degrees-of-freedom/scale convention

$$
V_e\sim\operatorname{IW}_2(\nu_{e0},\Psi_{e0}).
$$

The scale $\Psi_{e0}$ and explicit initial covariance $V_e^{(0)}$ must be
finite, symmetric, positive-definite $2\times2$ matrices in declared trait
order. The distribution is proper when $\nu_{e0}>1$. A finite prior mean
requires $\nu_{e0}>3$; Phase 4b records this property but does not impose it.

At the completed marker-sweep boundary the stored residual is

$$
E=Y-XA.
$$

No intercept or covariate is sampled internally, so the conditional uses all
$N$ likelihood rows:

$$
S_e=E^\top E,
$$

$$
V_e\mid-
\sim
\operatorname{IW}_2\left(
\nu_{e0}+N,
\Psi_{e0}+S_e
\right).
$$

It does not use $N-1$, $N-p$, marker counts, active-marker counts, or
traitwise scalar updates. Phenotypes are expected to be centred or otherwise
residualized by the caller when that is required by the intended analysis.

## Update order and authoritative state

The sampled-policy completed-iteration order is fixed as:

1. marker sweep using the current $V_b$, $\boldsymbol{\Pi}$, and $V_e$;
2. update $\boldsymbol{\Pi}$;
3. update $V_b$;
4. update $V_e$ from the complete current residual;
5. capture unthinned completed-iteration convergence state;
6. capture retained state at the Phase 3 retention indices.

The new draw is the sole residual-covariance state. It supplies the precision
for the next marker sweep and is copied unchanged to convergence, retained,
posterior-summary, and final-state outputs. There is no residual moment
replacement, diagonal approximation, or post hoc overwrite.

In `fixed_full` mode step 4 is absent. The residual statistic is not formed for
an update, the inverse-Wishart helper is not called, and no additional RNG is
consumed. Sampled mode intentionally consumes one private inverse-Wishart
transition per chain per completed iteration after the $V_b$ update.

## Resolved specification

The resolved model uses `model$residual_policy = "fixed_full"` or
`"sampled_full"`. Sampled mode records under `prior$residual_covariance`:

- `degrees_of_freedom`;
- `scale`;
- required `initial_value`;
- `sampled = TRUE` and `fixed_value = NULL`;
- `parameterization = "degrees_of_freedom_scale"`;
- `proper` and `finite_mean` status.

It also records `mcmc$update_flags$residual_covariance = TRUE` and update-order
version 2. Fixed mode retains the Phase 4a four-field fixed prior descriptor,
update flag `FALSE`, and update-order version 1. Partial or policy-inconsistent
argument sets fail before native dispatch.

## Raw-v2 contract

Sampled mode adds:

- `draws$residual_covariance` with axes
  `draw × chain × trait_row × trait_col`;
- `draws$convergence$residual_covariance` with axes
  `iteration × chain × trait_row × trait_col` when traces are retained;
- `posterior$residual_covariance_mean` with axes
  `trait_row × trait_col`;
- `final$residual_covariance` with axes
  `chain × trait_row × trait_col`;
- exact initial/prior/update metadata in the resolved input and qualification
  diagnostics.

All covariance states are finite, symmetric, and positive definite. The
posterior mean is computed only from retained sampled draws. When convergence
is retained, its last completed state equals the final state chain by chain.
Size-one draw and chain dimensions are not dropped.

Fixed mode continues to return a final fixed covariance while
`draws$residual_covariance` and
`posterior$residual_covariance_mean` are required-present `NULL` fields. Raw
validation enforces these policy-specific presence rules independently of the
compatibility identifier.

## Execution and validation boundary

Each joint chain owns its own `std::mt19937`, effects, patterns, residual,
$V_b$, $\boldsymbol{\Pi}$, $V_e$, captures, and diagnostics. Static scheduling
remains across chains only. Worker diagnostics and convergence capture consume
no RNG.

Qualification uses:

- an independent R conditional based on $\Psi_{e0}+E^\top E$ and
  $\nu_{e0}+N$;
- an independent completed-active R sampler with sampled
  $\boldsymbol{\Pi}$, $V_b$, and $V_e$;
- direct reconstruction of $E=Y-XA$ from decoded standardized BED genotypes;
- fixed-policy identity checks;
- repeated-run and one-worker/two-worker comparisons;
- raw-v2 presence, axes, SPD, posterior-mean, and final-state checks.

Phase 4b remains qualification-only. Public promotion, sampled fixed effects,
missing phenotypes, summary-statistic MT likelihoods, heterogeneous providers,
overlap-aware likelihoods, MT-BayesR/RC, regional covariance, templates,
mixtures, factor models, and large-$T$ pattern systems remain deferred.
