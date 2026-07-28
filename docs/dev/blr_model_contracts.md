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
- `learned_logistic`: learned inclusion/variance annotation coefficients;
- `annotation_probit_stick`: annotation-dependent BayesR component
  probabilities.

MT BayesRC factorizes annotation-dependent component probabilities and global
conditional non-null trait-pattern probabilities. Annotations do not alter
covariance matrices, component multipliers, operators, or residual likelihoods.

## State and effect meanings

`beta` is latent effect state where defined; `b` is the effective activity-
and-component-scaled effect. `d` is binary trait activity. BayesR/BayesRC
component code is a separate ordered state with zero as null. `dm` is posterior
trait activity/non-null probability, never a component mean.
