# MT BED extended convergence contract

## 1. Purpose

This document fixes the implementation contract for MT packed-BED convergence
scope beyond the current Tier 1 B/G/E diagonal diagnostics. It specifies future
Tier 2 covariance/probability and later Tier 3 selected-marker trace capture;
Phase 17W adds no traces or diagnostics.

## 2. Current implementation boundary

Only post-burn per-chain `B_diag`, `G_diag`, and `E_diag` traces exist. There is
no Tier 2/Tier 3 recorder, native route, public `extended` mode, output group,
or warning. Existing convergence mathematics and raw schema version 1 remain
unchanged.

## 3. Current sampler inventory

| Quantity | Iteration trace | Final value | Posterior mean | Owner/location | New trace | Storage | Tier |
|---|---:|---:|---:|---|---:|---|---|
| B diagonal | yes | B | covb | `vbs`, core | no | double | 1 |
| B off-diagonal | no | B | covb | `result.B` | yes | double | 2A |
| G diagonal | yes | G | covg | `vgs`, derived | no | double | 1 |
| G off-diagonal | no | G | covg | `genetic_covariance` | yes | double | 2A |
| E diagonal | yes | E | cove | `ves`, core | no | double | 1 |
| E off-diagonal | no | E | cove | `result.E` | yes | double | 2A |
| pi null/active mass | no | pi | pim | `result.pi` | yes | double | 2B |
| selected/all pi pattern | no | pi | pim | `result.pi` | yes | double | 2B |
| selected marker effective b | no | b | bm | `result.b` | yes | double | 3 |
| selected marker inclusion d | no | d | dm | `result.d` | yes | int32 | 3 |
| selected marker latent beta | no | beta | none | core-local `beta` | deferred | double | later |
| marker model-pattern index | no | none | none | marker draw | deferred | int32 | later |

Marker residual scores, all-marker traces, and full matrix/simplex diagnostics
are also deferred.

## 4. Canonical checkpoint

Future extended capture occurs once at the end of every completed Gibbs
iteration, after the marker sweep, post-sweep `samplePi`, final `sampleB`,
derived G calculation, and enabled E update, and immediately before the
post-burn covariance/pi accumulator block or the next iteration. In
`run_mt_bed_bayesc_core()` this is the boundary after the `ves` write and before
`if (iteration >= execution.nburn)` in `src/blr_mt_bed_core_impl.h`.

The current sequence is: per-set B/latent preparation; marker sweep updating
latent beta, effective b, d, residual, and model counts; thinned marker-summary
accumulation; pi update; final B update and `vbs`; derived G and `vgs`; E update
and `ves`; unthinned post-burn covariance/pi accumulators; next iteration.
Extended traces use every post-burn checkpoint and ignore `nthin`.

## 5. Tier hierarchy

- Tier 1 remains B/G/E diagonals.
- Tier 2A adds raw B/G/E off-diagonal covariance entries.
- Tier 2B adds active/null mass and explicitly selected or complete pattern
  probabilities.
- Tier 3 later adds explicitly selected marker effective b and binary d.

All-marker, latent-beta, residual-score, categorical-pattern, matrix-valued,
and simplex-valued joint diagnostics are excluded.

## 6. Covariance scalarization

The target is each raw covariance entry, matching public `vb`, `vg`, and `ve`.
Correlations, Fisher transforms, Cholesky factors, and partial correlations are
deferred. Tier 2A adds only off-diagonals; Tier 1 diagonals are never copied.

## 7. Lower-triangle ordering

For T traits, `Qoff=T*(T-1)/2`. Enumerate the strict lower triangle in
column-major order: column 1 with rows 2..T, column 2 with rows 3..T, and so on.
Names are `B[row_trait,column_trait]`, `G[...]`, and `E[...]`, preserving fitted
trait order rather than alphabetical order.

## 8. B semantics

With `updateB=TRUE`, B off-diagonals are stochastic and captured after the
post-sweep B update. With `updateB=FALSE`, stable rows have `not_updated`,
`captured=FALSE`, unavailable metrics, no trace allocation, and no flag solely
from being fixed.

## 9. G semantics

G is derived from phenotype and the current residual every iteration. Its
off-diagonals are eligible whenever T>1, marked `derived=TRUE`, and captured
even when B or E updates are disabled.

## 10. E semantics

For full residual covariance, updated off-diagonals are stochastic. Fixed E
uses `not_updated` without storage. Diagonal residual covariance includes
stable E off-diagonal summary rows with `structural_zero`, `structural=TRUE`,
`captured=FALSE`, all metrics unavailable, zero trace bytes, and no flags or
warning solely for structural status.

## 11. Probability timing

`samplePi(cmodel, result.pi, rng)` runs every iteration only when `updatePi`.
Probability capture is post-update at the canonical checkpoint. Model
probabilities are normalized by their sampler and validated/normalized at the
binding boundary. With `updatePi=FALSE`, probability rows are `not_updated`,
uncaptured, and unflagged.

## 12. Null/active mass

The model matrix must contain exactly one all-zero row. Store that row's
probability once under physical trace key `pi_mass:null_active`; expose the
primary diagnostic row as `pi_active=1-pi_null` and record `pi_null` as its
named complement metadata. The complement is not a second summary row, trace,
warning, or overview count.

## 13. Pattern probabilities

Pattern diagnostics are marginal scalar diagnostics. Dependence and the
sum-to-one constraint are expected; passing marginal rows does not establish
joint simplex convergence.

## 14. Pattern selection

Selection uses either public model names or one-based model-pattern indices in
one homogeneous vector. Names come from the existing authoritative model
names. Unknown/out-of-range/missing/duplicate selections fail and requested
order is preserved. A selected null pattern aliases the mass trace and scalar
result when mass is also requested, using one `diagnostic_key`.

## 15. All-pattern policy

All-pattern scope is explicit-only, never `auto` or the extended default, and
is never truncated. Capture is estimated before execution. A resolved hard
capture limit defaults to `min(2, max(0.25, 0.25*memory_warning_gb))` GiB
(`2` GiB when the ordinary threshold is infinite); a positive explicit
`trace_limit_gb` replaces it. Exceeding it is an error unless
`allow_large_traces=TRUE`. The ordinary memory warning remains independently
advisory. This protects against the capture/bundle/retention copies identified
by the Phase 17W feasibility grid.

## 16. Marker b

Tier 3 diagnoses effective b, not latent beta: b is zero when inactive and the
sampled effect when active. Names are `b[marker_id,trait]`, group `marker_b`.
Its zero-inflated, tied marginal posterior uses the existing rank-safe scalar
engine without jitter.

## 17. Marker d

d is the trait-specific binary inclusion state, named
`d[marker_id,trait]`, group `marker_d`. Constant all-zero/all-one traces are
`constant`, not perfect convergence; a constant chain against varying chains
is `constant_chain_mismatch`. Tail ESS may be unavailable when a quantile
indicator is constant. Mean ESS/MCSE concern posterior inclusion probability.

## 18. Marker selection

Exactly one of `marker_ids` and `marker_indices` is accepted. IDs resolve after
all BED-file, chromosome, column, and row alignment against final public marker
metadata; unknown or ambiguous duplicate source IDs fail. Indices are one-based
positions in that final selected marker order, not BED columns, global Glist
indices, or native indices. Missing, empty, duplicate, or out-of-range requests
fail; requested order is fixed before sampling. There is no posterior top-k.

## 19. Marker-trait scope

The initial Tier 3 contract diagnoses every requested marker for all fitted
traits in trait order, scaling as K*T. Pair-specific selection is deferred.

## 20. Binary diagnostic limitations

Average ranks handle ties deterministically. Constant indicators can make one
or both tail ESS values unavailable. No transition-count substitution or
jitter is introduced, and no constant trace is called converged.

## 21. Trace ownership

Every live chain exclusively owns its optional covariance, probability,
selected-b, and selected-d buffers and resolved integer selection maps. Workers
write no shared trace buffer and perform no R/Rcpp activity. After dispatch,
chain-private traces are copied in ascending logical-chain order into a
read-only deterministic bundle for the existing R engine.

## 22. Recorder architecture

The Phase 21 internal Tier 2 checkpoint should add Rcpp-free `MtBedExtendedTraceSpec` and
`MtBedExtendedChainTrace` types. The spec contains scope flags and resolved
model/marker indices. The chain trace contains post-burn draw count and separate
double covariance/probability/b buffers plus int32 d. It is optional, allocates
only post-burn storage, records selected K directly at the canonical checkpoint,
and is independent of compact chains and public trace retention.

## 23. Storage types

Covariance, probability, and marker b use IEEE double (8 bytes). Marker d uses
`std::int32_t` (4 bytes). R converts one d scalar at a time for diagnostics;
retained d arrays remain integer-valued.

## 24. Bundle schema

Use a new additive class `mtblr_extended_convergence_trace_bundle`, version 1,
rather than changing the existing Tier 1 bundle version. It contains `core`
plus separate typed `covariance`, `probability`, and `markers` sections, each
with descriptors and iteration-fastest, then chain, then quantity values.

## 25. Quantity descriptors

Descriptors contain `quantity_index`, `tier`, `group`, `trait_index`,
`trait2_index`, `marker_index`, `model_index`, `updated`, `captured`, `derived`,
`structural`, `storage`, and `diagnostic_key`. Native non-applicable indices use
`-1`. Public trait names, marker IDs, and model names are resolved after workers;
strings are not repeated per draw.

## 26. Summary schema

The convergence result and raw schema remain version 1. Existing identity
columns already represent `trait2`, `marker_id`, and `model_name`. Extended
groups are `B_cov`, `G_cov`, `E_cov`, `pi_mass`, `pi_pattern`, `marker_b`, and
`marker_d`. Descriptor metadata carries tier/capture/derivation/key/structure,
so the rectangular summary columns need not change. Any future need to change
columns requires an explicit convergence-result version decision first.

## 27. Public trace compatibility

Existing `fit$convergence_traces$scope`, draw count, core quantities, and core
values remain present. Add extended sections only under
`fit$convergence_traces$extended$covariance`, `$probability`, and `$markers`.
d stays integer; burn-in stays absent; no extra thinning is introduced.

## 28. Future public controls

Add `convergence="extended"`; `auto` remains Tier 1 only. The nested
`convergence_control$extended` accepts covariance `off_diagonal|none`,
probability `mass|none|selected|all`, `probability_models`, exactly one marker
selection representation, marker quantities from unique `b|d`, logical
`allow_large_traces`, and optional positive finite `trace_limit_gb`. Defaults
are off-diagonal covariance, mass probability, no patterns/markers, both marker
quantities once markers are selected, no override, and derived trace limit.

## 29. Conflict rules

Extended controls with another convergence mode fail. Selected probability
requires models; other probability modes forbid them. Marker IDs and indices
conflict. Explicit marker quantities without selection fail; selection with an
empty/unknown/duplicate quantity set fails. Unknown/duplicate control names,
malformed override, and invalid limits fail. Contradictions are never ignored.

## 30. Memory formulas

Let C chains, N post-burn draws, T traits, Qoff=T(T-1)/2, P patterns,
Psel selected unique non-reused patterns, and K markers. Capture bytes are:

```
B = updateB * 8*C*N*Qoff
G =           8*C*N*Qoff
E = full*updateE * 8*C*N*Qoff
mass = updatePi * 8*C*N
selected pi = updatePi * 8*C*N*Psel
all pi = updatePi * 8*C*N*P
b = 8*C*N*K*T
d = 4*C*N*K*T
```

Overlapping null traces are counted once. Structural/fixed rows need metadata
only. Workspace processes one scalar and is O(C*N), never multiplied by row
count. Summary rows are `3*T + 3*Qoff + probability rows + 2*K*T` when both
marker quantities are selected. Retention conservatively adds another captured
trace copy until ownership transfer is proven. Genotype and phenotype bytes are
not diagnostic memory. All figures are analytical, not measured RSS.

## 31. Large-request safety

All-pattern, large selected-pattern/marker, and retained-trace requests use the
same pre-execution analytical capture total and resolved hard limit in section
15. Above-limit requests fail before allocation unless explicitly overridden;
the override does not suppress the ordinary memory warning. There is no silent
truncation, selection reordering, or scope reduction.

## 32. Complexity

For S=C*N usable draws per scalar, sorting, FFT ESS, and workspace remain
O(S log S), O(S log S), and O(S); post-rank R-hat is O(S). Quantity counts are
3T for Tier 1, 3Qoff covariance descriptor rows, one mass row, selected/all
pattern rows after deduplication, and one K*T block per requested marker
quantity. Total time is approximately linear in scalar count; 1,000 markers
across many traits is not negligible.

## 33. Warning policy

At most one advisory warning reports overall status, total and per-tier/group
flag counts, extrema and quantities, mismatch count, constant/unavailable
selected-marker count, fewer-than-four-chain advisory, and
`fit$convergence`. It never enumerates all markers. Unrequested, structural,
and fixed rows do not affect warnings. `diagnostic_key` deduplicates complement
and overlap aliases in overview and warning counts.

## 34. Matrix and simplex limitations

Scalar covariance-entry diagnostics do not establish joint covariance-matrix
or positive-definite-geometry convergence. Marginal pattern diagnostics do not
establish joint simplex convergence. Selected-marker results say nothing about
unselected markers. No output uses definitive `converged` terminology.

## 35. Phase 21 internal Tier 2 checkpoint

This preserves the older Phase 17X specification, now superseded as an
immediate phase by the unified roadmap. Implement the binding-neutral
chain-private recorder, canonical checkpoint,
strict-lower B/G/E and active/null plus selected-pattern capture, extended
bundle construction/validation, and internal reuse of the existing scalar
engine. Do not change public formals or add marker traces.

## 36. Phase 21 public Tier 2 checkpoint

This preserves the older Phase 17Y specification as a later Phase 21
checkpoint. Activate `convergence="extended"`, nested Tier 2 controls, memory guard,
public Tier 2 summaries/optional traces, and deduplicated aggregated warnings,
preserving default Tier 1 behavior.

## 37. Phase 21 later Tier 3 checkpoint

In a separate contract/implementation phase add final-order marker selection,
chain-private selected b/d recording, high-dimensional safety, public summary
and trace formatting, and applicability tests. Do not combine all-marker or
latent-beta capture with this scope.

The immediate broader roadmap is Phase 18 unified STBLR/MTBLR alignment and
cleanup across CSR, block-eigen, and packed-BED; Phase 19 MTBLR BayesR across
supported operators; Phase 20 MTBLR BayesRC and annotation-informed models;
and Phase 21 extended diagnostics for both STBLR and MTBLR using this contract.
