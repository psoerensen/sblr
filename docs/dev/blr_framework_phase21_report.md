# Phase 21 report: extended STBLR and MTBLR diagnostics

## 1. Executive summary

Phase 21 adds the shared extended-mode contract and implements ST and MT
low-dimensional and explicitly selected-marker diagnostics. Every formal ST
quantity now comes from a genuine trait-by-chain native trace: BayesR mixture
probabilities, learned/group/probit-stick annotation parameters, and selected
effective effects, binary activities, and component states. No diagnostic is
reconstructed from pooled output, posterior means, final states, or reruns.

## 2. Repository baseline

Work began from `9baaaaa4f8e395636910478d17add67bcf22a83d` on `master` with a clean
tree. R 4.4.1/UCRT, GCC 13.2, OpenMP, 12 reported processors, and the default R
BLAS/LAPACK were used. Baseline source tests passed 6025 expectations with two
established opt-in skips. Baseline installed tests passed; package check had
zero errors, zero warnings, and three known notes (long Phase 17C path,
installed size, and legacy `std::cout` symbols).

## 3. Existing convergence mathematics

The Phase 17U/18 scalar engine is unchanged: Blom normalization, tied average
ranks, rank/folded split R-hat, FFT autocovariance, Geyer positive/monotone
sequences, stability-bounded bulk/tail/mean ESS, and MCSE of the mean.

## 4. Diagnostic inventory

The authoritative sampled/derived/fixed/structural inventory is in
`blr_extended_diagnostics_contract.md`. Native owners and update controls were
audited before trace capture was added.

## 5. Public convergence modes

All seven interfaces accept `auto`, `none`, `core`, and `extended`. Defaults
remain `auto`; automatic behavior remains core-only.

## 6. Extended convergence controls

The additive controls are `extended_groups`, `selected_markers`,
`selected_marker_quantities`, `full_probability_states`, `max_trace_gb`, and
`allow_large_traces`. Conflicts fail before native execution.

## 7. Tier structure

Tier 1 is the five core trait traces; Tier 2 is low-dimensional covariance,
probability, S, and annotation state; Tier 3 is explicit marker state.

## 8. Trace checkpoint and ownership

ST and MT capture occurs in logical-task-private native state at each completed
iteration, after enabled parameter updates and before the next iteration.
Burn-in is excluded, capture is unthinned and observational, and parameter
states updated on a schedule are repeated between update events.

## 9. Trace-bundle schema

The generic schema remains version 1. Descriptor identities were extended
additively and values remain iteration × chain × quantity.

## 10. Status vocabulary

The shared vocabulary and distinctions among fixed, inapplicable, structural,
constant, mismatch, and computed states are unchanged.

## 11. MT covariance diagnostics

Raw strict-lower B/G/E entries are captured without diagonal duplication. B/E
respect update controls; G is derived every iteration.

## 12. Covariance ordering

Ordering is `(2,1), (3,1), (3,2), ...` and public labels preserve trait order.

## 13. BayesC probability diagnostics

The ST adapter captures one active probability where its genuine per-chain
trace exists. MT captures the established pattern simplex.

## 14. BayesR probability diagnostics

ST BayesR/SBayesR captures the current global component simplex for CSR,
block-eigen, and packed-BED tasks in canonical component order. A binary
null/active simplex stores one physical scalar. MT captures component and
conditional-pattern marginals from each current joint probability state; full
joint states remain explicit opt-in.

## 15. BayesRC/SBayesRC probability diagnostics

Only conditional pattern probabilities are captured by the probability group.
Marker-average annotation priors are deliberately excluded.

## 16. Sampled selection-S diagnostics

Genuine ST task-private S traces are supported. Fixed S is `not_updated`,
absent S is `not_applicable`, and sampled MT S remains rejected.

## 17. Annotation coefficient diagnostics

ST and MT BayesRC/SBayesRC alpha and sigmaSqAlpha states are captured in
annotation-fast, stick-second order with explicit intercept identity.

## 18. Group/learned annotation diagnostics

Fixed-marker policies allocate no marker-expanded parameter traces. ST group
models capture current group inclusion probabilities and variance multipliers;
learned-logistic models capture distinct inclusion and variance coefficients.
Ordering and update ownership are native and deterministic.

## 19. Selected-marker contract

Only explicit homogeneous IDs or indices are accepted; shortcuts, duplicates,
missing values, unknown IDs, and out-of-range indices fail. Request order is
preserved.

## 20. Selected effective-effect diagnostics

ST and MT `b` are read directly from the current effective-effect state for
the pre-resolved selected indices only.

## 21. Selected binary-state diagnostics

ST and MT `d` are current binary activities per selected marker and trait,
never component codes. Ties are not jittered.

## 22. Selected component-state diagnostics

ST and MT BayesR/BayesRC use the established ordered component code. BayesC
component-only requests fail before native execution.

## 23. Descriptor identities and keys

Stable keys include family/model/operator/group/parameter and every applicable
identity dimension. Binary complements are not duplicated.

## 24. Retained traces

Retained bundles include schema, descriptors, quantity labels, and dimensioned
values at `fit$convergence_traces`. Compact chains remain independent.

## 25. Summary table

The table adds tier, parameter, trait pair, marker, component, pattern,
annotation, stick, intercept, update, structure, capture, and key columns.

## 26. Group and fit overviews

Group summaries cover core, covariance, probability, selection S, annotations,
and selected markers; fit statuses remain advisory `ok`, `warning`, `partial`,
`unavailable`, or `not_requested`.

## 27. Warning policy

At most one convergence advisory is emitted on the main R thread. Memory,
OpenMP fallback, and MAF/annotation overlap remain separate categories.

## 28. Memory accounting

Numeric/int32 capture, numeric conversion workspace, retained numeric traces,
descriptors, and summary output are separate analytical categories.

## 29. Hard memory guard

Capture plus retained-trace bytes are compared before sampling. Requests above
the threshold fail unless the explicit override is true, in which case one
warning is emitted without truncation.

## 30. ST execution integration

Trait × chain ownership now covers core traces, BayesC pi, BayesR component pi,
sampled S, low-dimensional group/learned/probit-stick parameters, and direct
selected-marker state. Allocation is proportional to requested quantities and
selected markers; no pooled trace is used for R-hat.

## 31. MT execution integration

CSR, block eigen, and BED receive one immutable resolved capture spec. Operator
and annotation preparation remain once per fit; workers write chain-private
plain C++ buffers.

## 32. Raw schema

Raw schema version 1 and model semantics version 2 remain unchanged. Empty
capture sections are additive and no Tier 2/3 values are populated under
`none`.

## 33. RNG neutrality

Focused ST and MT comparisons show exact posterior/final-state, probability,
annotation, and variance-output equality between `none` and extended capture.

## 34. Retention/thinning independence

Capture indexes every post-burn iteration independently of summary thinning,
compact-chain retention, and retained-trace output. Retention changes only
output and memory ownership.

## 35. Covariance reduction tests

Tiny diagonal-E tests verify strict-lower identity, fixed B status, derived G,
and structural E behavior.

## 36. Probability oracle tests

Tests verify ST/MT deterministic component/pattern/joint naming and ordering,
simplex normalization, binary-complement deduplication, update ownership,
opt-in joint state capture, and BayesRC exclusion of marker-average priors.

## 37. Selection-S tests

Synthetic task-private traces verify exact burn-in removal and no extra
thinning.

## 38. Annotation tests

ST and MT tests verify alpha/sigma ordering, learned coefficient identity,
group ordering, names, intercept identity, and every-iteration recording across
scheduled annotation updates.

## 39. Selected-marker oracle tests

Tiny CSR, block-eigen, and BED tests compare selected traces directly with
native task results and verify requested order, effective-b semantics, binary-d
coding, component applicability, and retained trace shape.

## 40. Reproducibility

Existing Phase 18–20 seed and OpenMP protections remain active; capture adds no
random draws.

## 41. Architecture audit

`blr_phase21_extended_diagnostics_architecture_audit.R` passes all 32 guards,
including native ST component, annotation/group, and selected-marker ownership.

## 42. Mutation sensitivity

`blr_phase21_mutation_sensitivity.R` detects all 62 mutations, including the 18
ST ownership, ordering, memory, thinning, retention, and RNG-neutrality repairs,
and does not modify the working tree.

## 43. Benchmark

`blr_phase21_extended_diagnostics.R` measures trace-plan resolution, analytical
capture/workspace/retention bytes, and scalar diagnostic compute time without
claiming linear scaling. The final 720-row grid required 0.07 seconds for plan
resolution and 0.11 seconds for scalar diagnostics; maxima were 54,816,000
capture bytes, 192,000 workspace bytes, 72,416,000 retained bytes, and 868,992
summary bytes. No linear-speedup claim is made.

## 44. Public API and documentation

All seven roxygen interfaces and maintained developer/Quarto notes describe
extended mode, explicit selection, and advisory interpretation.

## 45. Focused tests

The four central Phase 21 owners and the four repaired ST native owners passed
with zero failures, errors, warnings, or skips. Direct oracles cover component
pi, group/learned/probit-stick parameters, fixed-marker selected state, CSR and
block-eigen equality, packed-BED state, RNG neutrality, and thinning/retention
independence. Phase 18--20 focused protection owners and audits also passed.

## 46. Fast tier

The exact maintained fast filter, including the repaired BayesR and annotation
owners, passed with zero failures, errors, or warnings. Its only skip was the
established opt-in peak-RSS child-process case.

## 47. Full source suite

The required post-production-change source rerun passed all expectations with
zero failures, errors, or test warnings. The two established skips were
extended fresh-process reproducibility and peak RSS, both opt-in.

## 48. Installed tests

Installed tests from the built tarball passed with zero failures or warnings.

## 49. Package check

The final built-tarball check completed with zero errors and zero warnings.
Its three established notes were the long Phase 17C fixture path, installed
library size, and legacy `std::cout` use in scalar backends. An intermediate
check exposed undocumented internal `.convergence_spec` formals and a temporary
compile log; the owning roxygen blocks were corrected, the verified task log
was removed, and both findings are absent from the final check.

## 50. Diff and artifact hygiene

`git diff --check` passes after generated Rd whitespace normalization. Native
objects and the development DLL were removed after validation. The tracked
baseline source tarball was restored unchanged, no untracked tarball remains,
and no `.Rcheck` directory was created in the checkout.

## 51. Deviations and blockers

No scientific blocker remains. Fixed marker-specific priors intentionally have
no marker-expanded Tier 2 trace; BayesC component states remain inapplicable;
sampled MT selection-S remains unsupported. These are explicit model contracts,
not missing trace owners.

## 52. Recommended next phase

> consolidate permanent scientific test ownership, retire obsolete phase-specific scaffolding, review experimental and research-only routes, simplify maintained workflows and documentation, and stabilize the long-term BLR framework without changing posterior targets.

## 53. Readiness marker

PHASE 21 COMPLETE — EXTENDED STBLR AND MTBLR DIAGNOSTICS ACTIVE
