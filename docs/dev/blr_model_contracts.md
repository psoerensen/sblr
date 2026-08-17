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

$$
q_j(S)=\{2p_j(1-p_j)\}^{S+1}.
$$

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
  $q_j=\exp(X_j\theta)$ with global BayesC inclusion or global BayesR
  component probabilities.

MT BayesRC factorizes annotation-dependent component probabilities and global
conditional non-null trait-pattern probabilities. Annotations do not alter
covariance matrices, component multipliers, operators, or residual likelihoods.
Ordinary MT BayesR instead uses one global `joint_pi` simplex over complete
joint states; its component and pattern probabilities are marginals, not a
current factorized $P\times H$ parameterization.

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

## Approved target policies and Phase 1 boundary

Current native behavior above remains in force. The target model-policy
contract is defined by Sections 25--28 of
[the unified framework design](blr_unified_framework_design.md). Phase 1
resolves maintained ST wrapper semantics into that contract and converts only
source-understood ST raw fields. It does not change a statistical policy or
implement a joint-MT transition.

`single_trait` and `independent_traits` use scalar policies, with the latter
factorizing across traits. `joint_multitrait` uses an explicitly declared joint
state space and covariance policy; diagonal covariance alone does not establish
an independent-trait reduction. Probability inputs are named by their role:
inclusion mass, component mass, activity-pattern mass, or annotation policy.
They are not collapsed into generic `pi`.

## Phase 2 provider/operator boundary

Phase 2 adds shared global-marker, immutable operator-resource, and likelihood-
provider contracts around maintained ST routes. The resolved specification
owns each resource once; providers reference it by ID and carry their own
ordered trait set, local-to-global marker map, sufficient statistics, sample
size, likelihood regime, scales, and provenance. Missing provider markers
contribute no likelihood information and are not interpreted as null states or
zero effects.

This representation does not change a model policy or native transition.
Multiple independent providers and a coupled common-sample multi-trait BED
provider are structurally representable, but Phase 2 does not combine them in
a new posterior. Existing native CSR, BED, and block-eigen views remain the
authoritative operator implementations. Phase 3 layers versioned task seeds,
retention, and worker observations around qualified ST routes without changing
likelihood or model policies. The current MT covariance hybrid remains legacy
and is not routed through the Phase 3 scheduler.

Phase 3 is active for ordinary CSR, packed-BED, and block-eigen BayesC/BayesR,
plus fixed-marker and learned-logistic BayesC. Scheduled CSR, log-variance,
group, BayesRC/SBayesRC, and current MT routes retain version 0. A version-1
seed changes a deterministic trajectory but not the posterior target or the
within-chain transition.

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

The first corrected Cheng MT-BayesC$\Pi$ slice is implemented as an internal
qualification-only route. Phase 4a uses a supplied fixed symmetric
positive-definite full $V_e$ for common-sample BED, joint categorical updates
over $(0,0),(1,0),(0,1),(1,1)$, null collapse, conditional completion, a
Dirichlet activity-pattern update, and one authoritative inverse-Wishart
$V_b$. The exact contract and evidence are in
[the Phase 4a checkpoint](blr_phase4a_cheng_mt_bayesc_checkpoint.md).
For this joint activity policy, `posterior$pleiotropic_probabilities` is a
marker-axis posterior summary in $[0,1]$ and is required to equal the
markerwise probability of the declared all-traits-active pattern $(1,1)$,
identified by `1_1`. Models for which pleiotropy is not scientifically defined
do not manufacture or require this field.
Phase 4b extends the qualification route with sampled full inverse-Wishart
$V_e$. Its update is
$\operatorname{IW}_2(\nu_{e0}+N,\Psi_{e0}+E^\top E)$ after the complete marker
sweep, and the resulting matrix is the sole next-sweep, convergence, retained,
posterior-summary, and final residual-covariance state. The fixed policy skips
the update with zero residual-covariance RNG. See
[the Phase 4b checkpoint](blr_phase4b_sampled_residual_covariance_checkpoint.md).
Phase 5A generalizes the same qualification target to dynamic modest $T$. All
$2^T$ binary patterns are declared explicitly with the first trait changing
fastest, and arbitrary inactive coordinates are completed from the Gaussian
Schur conditional. The unique $(1,\ldots,1)$ row defines pleiotropic
probability from binary metadata rather than a two-trait string. General
inverse-Wishart propriety and finite-mean conditions are respectively
$\nu_0>T-1$ and $\nu_0>T+1$. See
[the Phase 5A checkpoint](blr_phase5a_general_t_cheng_mt_checkpoint.md).
Both computation and mandatory markerwise pattern probabilities scale as
$M2^T$. A checked preflight covers chain state, retained/convergence output,
native-to-R coexistence, raw construction, material copies, the shared aligned
packed-BED owner, and any conditional source-row decode buffer before provider
construction. Packed ownership contributes
$M\operatorname{round\_up}(\lceil N/4\rceil,64)$ bytes and is not multiplied by
chain count.
Dense $2^T\times2^T$ transition diagnostics are not scientifically required
and are not allocated or represented as posterior quantities.
Public promotion remains gated. Diagonal $V_e$ remains an explicit reduction;
independent no-overlap summaries do not identify off-diagonal residual
covariance. The current MT covariance hybrid is not the target and cannot be
migrated as authoritative covariance draws.

Phase 5B promotes the corrected complete-pattern Cheng policy through
`mtblr_bed()` for common-sample packed BED only. Public $V_b$ and sampled
$V_e$ priors use explicit inverse-Wishart degrees-of-freedom/scale arguments;
fixed $V_e$ remains a zero-update, zero-residual-covariance-RNG policy. The
legacy covariance hybrid is unreachable from public dispatch. MT BayesRC,
annotation, overlap, and missing-phenotype policies remain
unavailable. See
[the Phase 5B checkpoint](blr_phase5b_public_cheng_mt_checkpoint.md).

Phase 6A qualified, and Phase 6B publicly promotes, the independent-summary
residual policy through `mtblr_csr()` and `mtblr_block_eigen()`.
Each singleton-trait provider supplies a fixed positive scale $\phi_p$ and
contributes its own $s_p=X_p^\top y_p$ and $C_p=X_p^\top X_p$. Scores and
diagonal precisions add only over providers containing the current marker.
There is no full $V_e$ state in this likelihood, and declared provider overlap
is rejected rather than treated as independent. The Cheng pattern,
conditional-completion, Dirichlet, and authoritative inverse-Wishart $V_b$
contracts are unchanged. See
[the Phase 6A checkpoint](blr_phase6a_summary_mt_checkpoint.md).
The public boundary is recorded in
[the Phase 6B checkpoint](blr_phase6b_public_summary_mt_checkpoint.md).

Phase 7 publicly supports a factorized pattern-by-scale `bayesr` policy for
BED, CSR, and block-eigen MT routes. It retains the activity-pattern simplex
including one null state and adds a positive-scale simplex conditional on
non-null activity. Component assignments use `-1` only for collapsed null
markers. The base covariance statistic removes $\gamma_kq_j$ exactly once.
See [the Phase 7 checkpoint](blr_phase7_pattern_scale_mt_checkpoint.md).

Resolved Phase 1 objects accept only registered exact model-policy values and
declared descriptor structures. Compatibility identifiers do not relax model,
prior, state, covariance, probability, or axis validation. Scalar-variance and
matrix-covariance priors are checked under their declared parameterizations;
fixed full residual covariance must be finite, symmetric, positive definite,
and ordered by the resolved trait IDs.

For Phase 1 ST conversion, `posterior$pips` is the posterior non-null marker
summary, while inclusion- or component-probability parameter draws and final
states occupy separate explicit fields. BayesR/BayesRC component assignment
probabilities remain distinct from PIPs. If a native v1 route does not expose
a contracted chain state, the route remains v1 rather than receiving an
inferred value.
