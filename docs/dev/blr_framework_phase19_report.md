# Phase 19 — MTBLR BayesR and SBayesR

## 1. Executive summary

Phase 19 implements one MT pattern-by-component mixture contract for CSR,
block-eigen, and packed-BED operators. Its first completed validation used the
temporary meaning `sbayesr = BayesR + maf_s`; this continuation corrected the
public vocabulary before commit. Packed BED is `bayesr`; CSR and block eigen
are `sbayesr`; all share prior kernel `bayesr`, and fixed `selection_s` is an
independent option. The numerical mixture implementation is unchanged.

## 2. Repository baseline

The phase started clean on `9dc01acaf46edf8b5ad06723e8f5934cb44dec1e`
using R 4.4.1, GCC/G++ 13.2, OpenMP-capable Rtools44, and 12 reported logical
processors. The source suite and package check passed before production edits;
the two established optional skips and three classified package-check notes
were unchanged.

## 3. Existing contracts reused

Phase 18 naming, operator ownership, logical chains, worker-independent seeds,
static dispatch, deterministic aggregation, five core traces, convergence,
fit layout, memory vocabulary, and warnings are reused.

## 4. BayesR state space

For P non-null trait patterns and K positive components the model has exactly
`1 + P*K` states: one null plus every non-null pattern/component pair.

## 5. Deterministic state ordering

`null` is first, followed by supplied non-null pattern order and ascending
positive component order within each pattern.

## 6. Component multiplier contract

`mixture_var` defaults to `c(0, .01, .1, 1)`, contains one leading zero, and
then finite unique strictly increasing positive multipliers.

## 7. BayesR effect prior

The retained masked-latent MT BayesC covariance is multiplied by the known
component multiplier. Effective effects remain `pattern * beta`.

## 8. SBayesR scale policy

The S prefix denotes summary-statistics data, not MAF scaling. Both `bayesr`
and `sbayesr` use component scaling alone when `selection_s = NULL`; either
uses `[2p_j(1-p_j)]^(selection_s+1)` when the independent fixed-S option is
requested with aligned marker frequencies.

## 9. Sampled-S policy

One fixed scalar S is supported. `estimate_selection_s=TRUE` fails explicitly;
no unvalidated joint-MT MH transition is present.

## 10. Joint-state probability model

One named global probability vector uses joint-state counts and a positive
Dirichlet prior. Fixed-pi execution retains the supplied normalized vector.

## 11. B update under known scaling

Known component and MAF-S scale is removed from marker effects before sampling
the single shared base B. No component-specific covariance matrices exist.

## 12. Marker update

One shared binding-neutral kernel computes log weights, determinant and
quadratic terms, log-sum-exp probabilities, one state draw, and its conditional
latent effect. Null effects are exactly zero.

## 13. Initialization

`beta`, `b`, binary `state`, and zero-based `component` are validated jointly.
Component zero is equivalent to the null pattern; the default is all null.

## 14. Raw schema

`mtblr_raw` remains version 1. BayesR/SBayesR add model mixture metadata,
component final state/probabilities, and joint pi trace; BayesC rejects those
additions.

## 15. Public output

All Phase 18 fields remain. Additions are `component_final`,
`component_probabilities`, `pi_trace`, and `model_parameters$mixture`.

## 16. Component probabilities

The marker × component matrix includes null component zero and every row sums
to one. Pattern and component global pi marginals are cheap derived metadata.

## 17. Core variance traces

Every operator retains genuine `vbs`, `vgs`, `ves`, `vle`, and `vld` with
`vld = vgs - vle` at the completed iteration checkpoint.

## 18. Convergence

The unchanged shared Phase 18 engine diagnoses only those five post-burn,
unthinned per-chain quantities. Pi, components, S, markers, and off-diagonal
covariances remain outside core scope.

## 19. MT CSR implementation

One prepared CSR owner/view set is shared across deterministic joint chains;
the summary core invokes the shared mixture kernel and aligned LD diagonal.

## 20. MT block-eigen implementation

One reconstructed owner per required provenance is prepared before dispatch.
The same state order, priors, and summary core are used as CSR.

## 21. MT packed-BED implementation

One packed matrix and marker map are prepared per fit. The sample-space core
uses the same joint states with full or diagonal residual covariance.

## 22. Multichain execution

Each chain owns RNG, effects, states, residuals, accumulators, and traces.
Seeds and output order are independent of worker assignment; chain 1 supplies
final mutable state.

## 23. Memory

Estimates separately identify shared descriptors, private worker state,
per-chain results, component output, convergence capture/workspace, compact
chains, and retained convergence traces. They are analytical, not RSS.

## 24. BayesC reduction

The generalized conditional state contract reduces to BayesC with one unit
positive component and matching probabilities; existing BayesC references
remain protected. Trajectory identity is not claimed when categorical RNG
consumption differs.

## 25. One-trait ST/MT BayesR reduction

The shared component ordering, null semantics, effective-effect meaning,
operator diagonal, and LE/LD definitions reduce. Scalar versus joint
covariance-prior updates remain explicitly distinct where applicable.

## 26. One-trait ST/MT SBayesR reduction

The same qualification applies to the summary-statistics public model
`sbayesr`. Matched fixed S is tested independently of that public name.

## 27. S=-1 reduction

Exact tests across all three MT operators confirm `q_j(-1)=1` and equality of
all stochastic output fields between `selection_s = NULL` and explicit
`selection_s = -1`, while metadata preserves the semantic distinction.

## 28. CSR/block-eigen reduction

The unfiltered exact-reconstruction owner compares marker summaries, all five
variance traces, pi fields, and component probabilities at tight tolerance.

## 29. BED/summary reduction

The shared state and scale contracts are exercised on exactly derived data.
Posterior equality is not claimed when sample-space and marginal-summary
residual likelihoods are not mathematically identical.

## 30. Update-control reductions

Fixed B, E, and pi and updated configurations are covered; fixed quantities
retain `not_updated` convergence semantics.

## 31. Reproducibility

Phase 18 deterministic chain mapping is reused for default/explicit seeds,
serial/OpenMP dispatch, compact-chain retention, and convergence capture.

## 32. Architecture audit

The Phase 19 audit proves a single state descriptor/kernel, unique null,
single base-B scaling contract, three-operator support, version-one raw schema,
shared convergence, and mixture-memory categories. Every required guard passed,
as did the retained Phase 18 architecture audit.

## 33. Mutation sensitivity

Forty-nine scientific, state, operator, seed, output, convergence, schema,
scope, CI, and data-level/scale-separation mutations are checked by the
permanent audit. All 49 mutations were detected and restored; the retained
Phase 18 mutation audit also passed all 60 guards.

## 34. Benchmark

The benchmark records preparation-inclusive elapsed time, analytical mixture
memory, and fit size for both models on all three operators. Results are
regression signals and do not claim linear speedup. In the local small fixture,
elapsed times were 0.27/0.37/0.04 seconds for BayesR CSR/block-eigen/BED and
0.09/0.17/0.31 seconds for SBayesR. Estimated mixture bytes were 3,620, 3,636,
and 3,000 by operator; the benchmark reported 76 private-worker bytes, 908
per-chain-result bytes, and 108 component-output bytes in each case.

## 35. Public API and documentation

MT packed BED accepts `bayesc`/`bayesr`; MT CSR and block eigen accept
`sbayesc`/`sbayesr`. Each fit records prior kernel, individual versus summary
data level, effect-scale policy, selection-MAF provenance, and
`model_semantics_version = 2`. Reference-panel MAF fallback is opt-in.
BayesRC and annotations are not accepted.

## 36. Focused tests

Focused scientific, operator, convergence, raw-schema, and Phase 18 protection
owners passed with zero failures, errors, or test warnings. These included the
new state-kernel/model and three-operator owners plus Phase 17R/S/U/V, Phase 18,
unified convergence, raw-schema, reduction, and reproducibility protection.

## 37. Fast tier

Fast CI owns the state oracle, small operator/S reductions, public/raw and
shared convergence contracts, plus package checking. The exact workflow filter
passed in 109.1 seconds with zero failures or errors and the established opt-in
peak-RSS skip.

## 38. Full source suite

The final source suite passed in 158.8 seconds with zero failures, zero errors,
zero test warnings, and two established opt-in skips (fresh-process extended
coverage and peak-RSS measurement).

## 39. Installed tests

Installed-package tests passed from the built tarball.

## 40. Package check

The final replacement built-package check passed in 467.7 seconds with zero errors and
zero warnings. Its three classified pre-existing notes were the long Phase 17C
fixture path, installed DLL size, and legacy `std::cout` symbols. An initial
completion check exposed one missing packed-BED BayesR roxygen argument block;
that documentation defect was corrected, Rd validation passed, and the
replacement check was clean. During the final semantics continuation, the first
check attempt also exposed the missing block-eigen selection-MAF Rd arguments
and 20 installed historical expectations that still used `bayesc`/`summary`
for annotation-summary routes. Those contract migrations were corrected and
focused-tested before the clean replacement check; the intermediate
`1 ERROR, 1 WARNING, 3 NOTEs` result is therefore retained here rather than
hidden.

## 41. Diff and artifact hygiene

`git diff --check` and EOL-insensitive diff inspection pass. Generated compiler,
DLL, and check-directory artifacts are absent; the tracked baseline package
tarball is unchanged and no task-generated tarball change remains.

## 42. Deviations and blockers

Sampled MT S is explicitly unsupported. Full ST/MT posterior equality is not
asserted where scalar and joint covariance priors or likelihoods differ.

## 43. Recommended next phase

> implement MTBLR BayesRC and SBayesRC as annotation-informed extensions of the Phase 19 pattern-by-component model, initially using factorized active, trait-pattern, and component probability policies across the aligned CSR, block-eigen, and packed-BED operators.

## 44. Readiness marker

PHASE 19 COMPLETE — MTBLR BAYESR/SBAYESR DATA-LEVEL SEMANTICS ALIGNED
