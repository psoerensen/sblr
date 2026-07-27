# MTBLR BayesRC/SBayesRC contract

## Scope

Phase 20 extends the Phase 19 MT trait-pattern-by-component sampler with
annotation-dependent component probabilities.  `bayesrc` is the individual
`packed_bed` model and `sbayesrc` is the summary-statistics `csr` and
`block_eigen` model.  Both use `prior_kernel = "bayesrc"` and
`annotation_policy = "annotation_probit_stick"`.  The `s` prefix denotes the
data level; MAF-dependent scaling is controlled only by `selection_s`.

## State space and prior

The Phase 19 descriptor is unchanged: state zero is the unique null state,
followed by supplied non-null trait patterns in supplied order and positive
components in ascending multiplier order.  For marker annotation row `A_j`,
the established probit stick-breaking convention produces
`theta_j = (theta_j0, ..., theta_jK)`.  A separate global vector `omega`
contains conditional non-null trait-pattern probabilities.  Thus

```
P(null | A_j) = theta_j0
P(pattern p, component k | A_j) = theta_jk * omega_p, k > 0.
```

Both vectors are normalized and the factorized joint prior sums to one.  There
is no repeated null state and no annotation coefficient indexed by trait
pattern.  Annotations affect state probabilities only; the Phase 19 base `B`,
component multipliers, optional `q_j(selection_s)`, residual covariance, and
operator are unchanged.

## Annotation input and preprocessing

`annotations` is a finite numeric matrix or data frame with unique nonempty
column names.  External rows require unique marker identifiers in row names;
rows are reordered only by an exact match to final public marker IDs.  Missing,
extra, duplicated, or ambiguous identifiers fail before native execution.
By-construction row order is accepted only when explicitly established by the
adapter.  Preprocessing occurs once per fit, never mutates the input, preserves
an explicit intercept, optionally adds one intercept, identifies binary
columns deterministically, standardizes continuous columns when requested,
and records centers, scales, names, binary status, and alignment.

## Coefficients and updates

There is one chain-private annotation coefficient matrix with one column per
stick, shared across traits.  The established ST probit augmentation is used:
markers reaching a stick receive correctly truncated latent-normal draws,
coefficients receive Gaussian conditional updates, non-flat coefficient
variances receive the established scaled inverse-chi-square update, and the
first processed column is flat only when `intercept_flat = TRUE`.
`alpha_update_every` is a positive iteration schedule independent of thinning,
chain retention, and workers.  `updateAlpha = FALSE` holds coefficients and
variances fixed.  `updatePi` updates only `omega` from active-marker pattern
counts under a positive Dirichlet prior.

## Canonical iteration checkpoint

All operators use: marker pattern/component/effect sweep; optional `omega`
update; scheduled optional annotation update; optional base-`B` update;
genetic covariance calculation; optional `E` update; completed-iteration
`vbs`, `vgs`, `ves`, `vle`, and `vld` capture; posterior accumulation.
Trace capture is post-burn and unthinned for convergence, while posterior means
retain the ordinary `nthin` contract.

## Output

The common marker, component, covariance, variance-trace, convergence, chain,
memory, and semantics fields remain unchanged.  BayesRC-specific values live
in `fit$model_parameters$annotations`: preprocessing metadata, coefficient and
variance final/mean values, explicit `pattern_pi_*`, component/stick/pattern
names, posterior-mean `prior_component_probabilities`, and
`prior_active_probabilities`.  `component_probabilities` continues to mean
posterior component-state probabilities.  Global joint-state `pi_*` is not
redefined; it is present but `NULL` for BayesRC/SBayesRC.

## Execution, memory, and limitations

Operator and annotation preparation occur once.  Logical chains own RNG,
marker state, annotation state, latent workspace, and accumulators.  Immutable
operator data, processed annotations, and state descriptors are shared; no R
or Rcpp work occurs in OpenMP workers.  Analytical pre-execution memory counts
shared annotations, private annotation state/workspace, per-chain results,
posterior accumulators, and formatted prior-probability output separately.

Phase 20 does not implement annotation-dependent trait-pattern probabilities,
annotation-dependent covariance matrices or multipliers, sampled MT
`selection_s`, formal annotation diagnostics, Tier 2 diagnostics, selected
marker traces, LD-swap, prediction, fine-mapping moves, or factor models.
