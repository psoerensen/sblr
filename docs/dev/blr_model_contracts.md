# BLR model contracts

## Model semantics

The public model name combines prior kernel and data level:

| Model | Prior kernel | Data level | Baseline effect scale |
|---|---|---|---|
| `bayesc` | BayesC | individual | unit |
| `sbayesc` | BayesC | summary statistics | unit |
| `bayesr` | BayesR | individual | component |
| `sbayesr` | BayesR | summary statistics | component |
| `bayesrc` | annotation-informed BayesR | individual | component |
| `sbayesrc` | annotation-informed BayesR | summary statistics | component |

The `s` prefix means summary statistics and never activates MAF scaling.
`model_semantics_version = 2` records this convention.

## Selection-S policy

`maf_effect_s = NULL` means no MAF scale. A finite scalar activates

\[
q_j(S)=\{2p_j(1-p_j)\}^{S+1}.
\]

`maf_effect_s = -1` records explicit user intent but gives a unit numerical
scale. BayesC changes from `unit` to `maf_s`; BayesR/BayesRC changes from
`component` to `component_maf_s`. Valid scalar routes may sample S with the
documented bounded prior/proposal controls. MT sampled S is unsupported.

Frequency resolution is explicit: aligned `effect_maf`, GWAS-summary MAF,
analysis-genotype MAF by construction, then reference MAF only under
`allow_reference_maf_for_maf_effect_s = TRUE`. Source, population, alignment,
and fallback metadata are mandatory when the scale is active.

## Probability policies

Global BayesC/BayesR models use a global binary or mixture/simplex policy.
Annotation policies are:

- `fixed_marker`: fixed marker inclusion probabilities/scales;
- `group`: sampled group-level BayesC parameters;
- `learned_logistic`: learned inclusion/variance annotation coefficients. Its
  inclusion provider is the unclipped offset-logistic model; obsolete
  `pi_min`/`pi_max` controls and ambiguous abbreviations are rejected before
  dispatch, while exact `pi_marker` remains the supported initial global
  probability control. Arguments forwarded to this provider through `...`
  must have unique, nonempty names; unnamed, empty-name, `NA`-name, mixed, and
  duplicate forwarding is rejected before positional or ambiguous matching.
  The global BayesC inclusion
  probability has a nonconjugate logit-scale update when offsets are nonzero
  and the ordinary Beta reduction when they are exactly zero;
- `annotation_probit_stick`: annotation-dependent BayesR component
  probabilities.
- `log_variance`: learned marker-specific relative non-null variance
  \(q_j=\exp(X_j\theta)\) with global BayesC inclusion or global BayesR
  component probabilities.

MT BayesRC factorizes annotation-dependent component probabilities and global
conditional non-null trait-pattern probabilities. Annotations do not alter
covariance matrices, component multipliers, operators, or residual likelihoods.
Ordinary MT BayesR instead uses one global `joint_pi` simplex over complete
joint states; its component and pattern probabilities are marginals, not a
current factorized \(P\times H\) parameterization.

Current `log_variance` support is ST CSR and retained block eigen with
the public `sbayesc` and `sbayesr` methods. BED-LV, MT-LV, and simultaneous learned
probability-plus-LV variance are not implemented. See
[the annotation-prior capability matrix](annotation_prior_architecture_matrix.md)
for current and proposed scope.

## State and effect meanings

`beta` is latent effect state where defined; `b` is the effective activity-
and-component-scaled effect. `d` is binary trait activity. BayesR/BayesRC
component code is a separate ordered state with zero as null. `dm` is posterior
trait activity/non-null probability, never a component mean.

## Approved Phase 0 target policies

Current behavior above remains in force until migrated. The target model-
policy contract is defined by Sections 25--28 of
[the unified framework design](blr_unified_framework_design.md) and is not yet
implemented.

`single_trait` and `independent_traits` use scalar policies, with the latter
factorizing across traits. `joint_multitrait` uses an explicitly declared joint
state space and covariance policy; diagonal covariance alone does not establish
an independent-trait reduction. Probability inputs are named by their role:
inclusion mass, component mass, activity-pattern mass, or annotation policy.
They are not collapsed into generic `pi`.

For marker covariance,

$$
V_b\sim\operatorname{IW}_T(\nu_b,\Psi_b)
$$

means a density proportional to

$$
|V_b|^{-(\nu_b+T+1)/2}
\exp\left\{-\frac12\operatorname{tr}(\Psi_bV_b^{-1})\right\}.
$$

The canonical fields are `degrees_of_freedom = nu_b` and `scale = Psi_b`.
The distribution is proper for $\nu_b>T-1$; its mean is finite only for
$\nu_b>T+1$. A mean-calibration helper uses and records
$\Psi_b=(\nu_b-T-1)M_b$. No ambiguous prior-covariance shortcut crosses the
native boundary.

The first corrected Cheng MT-BayesC$\Pi$ slice is staged. Phase 4a uses a
supplied fixed symmetric positive-definite full $V_e$ for common-sample BED.
Phase 4b samples full inverse-Wishart $V_e$ and must finish before promotion as
the maintained MT replacement. Diagonal $V_e$ remains an explicit reduction;
independent no-overlap summaries do not identify off-diagonal residual
covariance. The current MT covariance hybrid is not the target and cannot be
migrated as authoritative covariance draws.
