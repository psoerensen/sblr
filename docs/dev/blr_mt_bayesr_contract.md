# MTBLR BayesR and SBayesR contract

## Purpose and reused contracts

Phase 19 extends the Phase 18 joint small-trait BayesC model without changing
its operator likelihood, masked-latent/effective-effect semantics, covariance
update ordering, chain topology, convergence mathematics, or raw schema
version. CSR and block eigen continue to share the summary-operator core;
packed BED retains its sample-space residual core. The three adapters share the
state specification and marker-state probability implementation.

The public component-multiplier name is `mixture_var`, matching STBLR. Its
default is `c(0, 0.01, 0.1, 1)`. Values are finite, unique, and strictly
ascending; the first and only zero is the null multiplier and at least one
positive multiplier is required.

## Joint state space and ordering

Let the supplied unique binary trait patterns contain one unique null pattern.
Let `P` be the number of non-null patterns and `K` the number of positive
components. There are `J = 1 + P*K` joint states, ordered as:

1. `null`;
2. supplied non-null pattern order;
3. within each pattern, ascending positive-component order.

Names are `null` and
`<pattern_name>__component_<k>`. There is no null state per pattern. The native
descriptor stores the binary pattern, zero-based component, multiplier, and
name index. Public names never expose native indices. More than 4096 joint
states fails before native execution; no truncation occurs.

## Prior and marker update

Phase 18 uses latent `beta_j ~ N(0, B)` and effective
`b_j = D(a_j) beta_j`. Phase 19 preserves that masked-latent interpretation.
For positive component multiplier `gamma_k` and marker scale `q_j`, the latent
prior is

`beta_j | a_j,k,B,S ~ N(0, gamma_k q_j B)`.

Both BayesR and SBayesR use `q_j = 1` when `selection_s = NULL`; the S prefix
denotes summary-statistics data. When a fixed `selection_s` is supplied, either
public model uses `h_j = 2 p_j (1-p_j)` and `q_j(S) = h_j^(S+1)`. One fixed
scalar S is shared by the joint MT model. `selection_s = -1` therefore gives
`q_j = 1` exactly while remaining semantically distinct from a null request.
Sampled S is intentionally unsupported in Phase 19: deriving and validating a
joint shared-S MH transition is deferred, and `estimate_selection_s = TRUE`
fails explicitly.

Public mapping is `bayesr` for packed BED and `sbayesr` for CSR/block eigen.
All three dispatch the same `bayesr` prior kernel. New raw/fit metadata records
model-semantics version 2, data level, prior kernel, effect-scale policy, and
selection-MAF provenance; reference-panel MAF fallback is opt-in.

For every joint state the shared kernel uses prior precision
`B^{-1}/(gamma_k q_j)`, the operator-specific score and residual precision,
the Cholesky log determinant, and the posterior quadratic term. It normalizes
in the log domain, draws one categorical state, and samples one latent effect.
The null state has multiplier zero and returns exact zero latent and effective
effects. Pattern, component, and effective state are updated atomically.

## Base covariance and probabilities

`B` remains one shared base covariance. Known component/S scale is removed from
active marker cross-products before B updates; no component- or pattern-specific
B matrices exist. Null states contribute no marker covariance information.

One global probability vector is defined over the J joint states. With
`updatePi = TRUE`, a Dirichlet-multinomial update uses joint-state counts and a
positive `joint_pi_prior`; with `updatePi = FALSE` the initialized vector is
fixed. The default prior is one for every state. `joint_pi` initializes the
joint vector. BayesC retains `pimodels`; fields are never overloaded across
methods.

## Initialization and output

BayesR/SBayesR accept `beta`, `b`, binary `state`, and integer `component`.
Component zero is equivalent to the null pattern; positive components require
a non-null supplied pattern. Inactive effective effects are zero. The default
is the all-null, zero-effect state.

The additive BayesR raw and fit fields are `component_final` (zero is null) and
`component_probabilities` (marker by component, including zero, rows summing to
one). `dm[j,t]` is the posterior probability that marker j is active for trait
t, marginalized over components and patterns. `pi_final`, `pi_mean`, and
`pi_trace` are named joint-state probabilities. `model_parameters$mixture`
contains multipliers, names, supplied patterns, joint-state names, derived
pattern/component probability marginals, and fixed-S metadata. A marker by
joint-state posterior matrix is not retained.

## Operators, traces, and chains

CSR uses the aligned LD diagonal, block eigen uses its represented diagonal,
and packed BED uses selected-marker heterozygosity and its actual standardized
marker map. All retain every-iteration `vbs`, `vgs`, `ves`, `vle`, and `vld`,
with `vld = vgs - vle`. Formal diagnostics reuse the Phase 18 scalar engine
only. Joint pi, components, S, markers, and covariance off-diagonals are not
core convergence quantities.

Every route prepares immutable operator data once, dispatches one complete
joint model per logical chain with worker-independent seeds, aggregates in
chain order, pools retained posterior draws, takes final mutable state from
chain 1, and captures convergence traces independently of compact-chain
retention. BayesC takes the unchanged Phase 18 path.

## Validation, memory, and failure policy

The three public signatures use the same ordered controls: `mixture_var`,
`joint_pi`, `joint_pi_prior`, `component`, and the established selection-S
controls. BayesR rejects every S-only control. SBayesR requires aligned allele
frequencies strictly inside (0,1). Sampled S, BayesRC/SBayesRC, annotations,
and joint-probability convergence diagnostics fail or remain unavailable
rather than being silently mapped to another model.

The analytical estimate separately charges shared state descriptors, private
worker state, per-chain component/pi results, and formatted component output.
The `J <= 4096` state guard is evaluated before native allocation. No pattern,
component, chain, or output is truncated to satisfy memory thresholds. Native
failures retain logical chain, operator, and state context and never return a
partial aggregate.
