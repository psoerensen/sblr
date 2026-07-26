# Unified BLR Framework Phase 17W report

## 1. Executive summary

The extended MT BED convergence contract is implementation-ready. Phase 17W
defines Tier 2 covariance/probability and separately staged Tier 3 selected-
marker tracing without adding any production trace or output.

## 2. Repository baseline

The audit started from clean `master` at
`227c785e1faea0c624337765307506bad472da43`. R 4.4.1 UCRT used GCC/G++/Fortran
13.2 with `-fopenmp`; runtime reports 12 processors. R's bundled BLAS/LAPACK
were linked and relevant BLAS/OpenMP environment variables were unset. Hosted
CI status was not visible locally. Baseline source tests passed 5,523
expectations (5,521 passes, two established skips), and the built-package check
had zero errors/warnings and three established notes.

## 3. Phase 17V verification

Public `mtblr_bed()` retains `auto|none|core`, exactly-one route selection,
Tier 1 B/G/E diagonals, validated summaries, optional post-burn core traces,
and one aggregated advisory warning. Production hashes are unchanged.

## 4. Extended quantity inventory

B/G/E off-diagonals, probabilities, selected b/d, latent beta, and categorical
pattern state exist only as iteration-local/final/posterior-summary values; none
has an iteration trace. The complete availability/owner/storage table is in the
extended contract section 3.

## 5. Sampler sequence

`run_mt_bed_bayesc_core()` performs per-set B/latent preparation, marker sweep
(beta/b/d/residual/model counts), thinned marker accumulation, pi update, final
B update/trace, derived G/trace, E update/trace, then unthinned post-burn
covariance/pi accumulation.

## 6. Canonical checkpoint

One future checkpoint is fixed after all enabled pi/B/G/E updates and their
current diagonal writes, immediately before the post-burn accumulator block or
next iteration. Selected b/d already represent the completed marker sweep.

## 7. Tier hierarchy

Tier 1 remains diagonals; Tier 2A is raw covariance off-diagonals; Tier 2B is
active/null mass and explicit pattern scope; Tier 3 is explicit selected-marker
effective b and binary d. Latent/all-marker/categorical/matrix/simplex scope is
deferred.

## 8. Covariance parameterization

Raw covariance entries were selected because they map directly to `vb`, `vg`,
and `ve`. Correlation, Cholesky, partial-correlation, and Fisher-z diagnostics
are deferred.

## 9. Lower-triangle ordering

For `Qoff=T(T-1)/2`, descriptors traverse strict-lower columns then rows. Names
use fitted trait order: `B[row,column]`, `G[row,column]`, `E[row,column]` with
groups `B_cov`, `G_cov`, `E_cov`.

## 10. B update semantics

Updated B is stochastic and captured after final `sampleB`; fixed B rows are
`not_updated`, uncaptured, unflagged, and allocate no trace.

## 11. G derivation semantics

G is derived each iteration from phenotype and current residual and is eligible
for T>1 regardless of B/E update controls.

## 12. E mode semantics

Full updated E is stochastic; fixed full E is `not_updated`. Diagonal-mode
off-diagonals are stable `structural_zero` rows: uncaptured, unavailable,
unflagged, and zero-byte.

## 13. Probability update semantics

Pi is captured after `samplePi` at the canonical checkpoint. Fixed pi uses
`not_updated` without storage. The binding/model contract guarantees one null
pattern and normalized probabilities.

## 14. Null/active deduplication

One physical null trace under `pi_mass:null_active` supports primary row
`pi_active`; `pi_null` is complement metadata. It contributes one summary,
overview, trace, and warning key.

## 15. Pattern-probability scope

Selected/all rows are scalar marginal simplex components only. They do not
establish joint model-probability-vector convergence.

## 16. Pattern selection

Exactly one homogeneous vector of authoritative public names or one-based
pattern indices is accepted. Missing, unknown, out-of-range, and duplicate
requests fail; order is preserved; null overlap aliases existing storage.

## 17. All-pattern safety

All-pattern scope is explicit-only. Pre-execution capture above the resolved
hard limit requires `allow_large_traces=TRUE`; there is no truncation or scope
reduction. The ordinary memory warning remains separate.

## 18. Marker-effect scope

Tier 3 targets effective b (zero inactive, sampled effect active), never latent
beta, marker score, or automatic top-k.

## 19. Inclusion-state scope

Tier 3 d is the trait-specific 0/1 inclusion state; mean ESS/MCSE diagnose the
posterior inclusion-probability estimate.

## 20. Marker selection

IDs resolve against final selected marker metadata and ambiguous IDs fail.
Indices are one-based final selected-order positions. Representations conflict;
unknown/missing/duplicate/out-of-range inputs fail; request order is fixed.

## 21. Marker-trait scope

Every requested marker is initially diagnosed for all traits in fitted trait
order. Pair-specific selection is deferred.

## 22. Binary applicability

Deterministic binary experiments retained tie-safe behavior. Constant binary
draws returned `constant`; a constant chain against varying chains returned
`constant_chain_mismatch`; rare nonconstant inclusion supported mean ESS/MCSE,
while tail ESS can be unavailable when quantile indicators are constant.

## 23. Trace ownership

Optional buffers and resolved selection maps are chain-private. Workers write
no shared traces and use no R/Rcpp. Deterministic post-dispatch bundling follows
logical-chain order.

## 24. Recorder and bundle

Future Rcpp-free spec/chain-trace types capture post-burn only at one checkpoint.
A new version-1 extended bundle contains the unchanged core bundle and separate
covariance, probability, and marker sections in iteration/chain/quantity order.

## 25. Storage types

Covariance, probability, and b use 8-byte double. d uses 4-byte int32 and
remains integer in retained R output.

## 26. Output compatibility

Raw and convergence result remain version 1; current summary identity columns
support new groups. Descriptor metadata owns tier/captured/derived/structural/
diagnostic-key fields. Existing core trace fields remain stable, with extended
sections additive below `$extended`.

## 27. Future public API

The selected design adds `convergence="extended"` and nested controls:
off-diagonal/none covariance; mass/none/selected/all probability; selected
models; one marker selection representation; b/d quantities; large-trace
override; and optional limit. `auto` never becomes high-dimensional.

## 28. Conflict validation

Extended controls require extended mode. Selected probability requires models;
other modes forbid them. Marker representations conflict; explicit quantities
require selection. Unknown/duplicate controls, invalid logicals/limits, empty
or duplicate quantity sets all fail closed.

## 29. Memory formulas

Exact capture, reuse, retention-copy, O(CN) workspace, and summary-row formulas
are fixed in contract section 30. Diagnostics exclude genotype and phenotype
memory and do not scale with `ncores` beyond chain ownership.

## 30. Large-request policy

With default `memory_warning_gb=8`, the hard capture limit resolves to 2 GiB;
otherwise it is `min(2,max(.25,.25*memory_warning_gb))` GiB, or an explicit
positive `trace_limit_gb`. The 27,000-case grid found 3,764 cases over 2 GiB,
including approximately 2.06--2.08 GB selected-marker cases, supporting a
separate pre-allocation override gate. Maximum capture was 224.911 GiB.

## 31. Complexity

Per scalar, ranks and FFT ESS are O(CN log(CN)), R-hat is O(CN), and workspace
is O(CN). Counts scale as 3T, 3Qoff, deduplicated probability rows, and K*T per
marker quantity. Time therefore scales approximately with scalar count.

## 32. Warning policy

One advisory warning summarizes total and group/tier flags, extrema, mismatches,
selected-marker constant/unavailable counts, chain advisory, and result path.
Diagnostic keys deduplicate complements/overlaps. Structural, fixed, and
unrequested rows never warn solely from those states.

## 33. Matrix/simplex limitations

Scalar covariance rows do not establish joint SPD-matrix convergence; marginal
pi rows do not establish joint simplex convergence; selected markers do not
characterize unselected markers.

## 34. Deterministic contract fixtures

The 79-expectation owner passed nt=1/2/3/5 covariance ordering and update modes,
null/active and selected/all pattern cases, ID/index marker validation,
continuous/bounded/zero-inflated/binary applicability, overlap reuse, exact
memory formulas, safety decisions, conflicts, and grouped warnings.

## 35. Optional reference checks

With optional `posterior` 1.6.1, covariance, bounded probability,
zero-inflated, and binary fixtures matched the existing engine: maximum
absolute difference `7.11e-15`, maximum relative difference `2.08e-16`, and
binary difference zero. No runtime dependency was added.

## 36. Feasibility benchmark

The 27,000 analytical cases covered C=2/4/8, N=100/1000/5000, T through 50,
patterns through 16,384, selected patterns through 1,000, markers through
10,000, b/d/both, and retention. Maximum all-pattern storage was 4.88281 GiB;
maximum capture/retention was 224.911 GiB. Representative scalar engine times
ranged from below timer resolution to 0.14 seconds at C=8,N=5000. These are
synthetic regression signals, not production runtime claims.

## 37. Staged implementation

The earlier proposal for Phase 17X internal Tier 2 and Phase 17Y public Tier 2
is superseded. Phase 18 first aligns and simplifies STBLR and MTBLR across CSR,
block-eigen, and packed-BED operators. Phase 19 adds MTBLR BayesR across
supported operators; Phase 20 adds MTBLR BayesRC and annotation-informed
models; Phase 21 implements extended diagnostics for both STBLR and MTBLR,
using this MT contract as an input. Within Phase 21, retain separate internal
Tier 2, public Tier 2, and later selected-marker Tier 3 checkpoints.

## 38. Existing-route protection

Phase 15/17 owners, public MT/scalar interfaces, packed-BED families, summary MT,
raw schemas, native core, wrappers, registration, fixtures, and public Rd remain
protected and unchanged.

## 39. Mutation sensitivity

All 52 required ordering, semantic, ownership, safety, compatibility, and
production-protection mutations are detected.

## 40. Installed-check behavior

Ordering, deduplication, selection, naming, memory, safety, binary applicability,
control, and warning helpers are portable. Only static checkout/hash assertions
skip narrowly in installed-only contexts.

## 41. Tests and CI

Phase 17W is included in the exact fast filter. Phase 17W itself passed 79
expectations. The protected-owner filter passed 2,165/2,165 expectations. The
exact fast tier passed 2,743 expectations with one established opt-in skip and
no failures, warnings, or errors. The full source suite passed 5,600
expectations with two established opt-in skips and no failures, warnings, or
errors. The installed test suite completed successfully.

## 42. Package check

The built-tarball check completed with zero errors and zero warnings. Its three
notes are unchanged and classified: the established Phase 17C long fixture
path, installed size (5.2 MB; 4.2 MB in `libs`), and legacy scalar-backend
`std::cout` objects. The first check invocation lacked the Rtools path and
failed installation; the required rerun under the recorded GCC 13.2 Rtools
environment passed.

## 43. Diff hygiene

Final ordinary and ignore-EOL statistics agree. Production, wrappers,
registration, DESCRIPTION, NAMESPACE, Rd, fixtures, and the tracked baseline
tarball have no diff. The 24 compiler products produced during validation were
removed; no DLL, object, `.Rcheck`, or temporary product remains in the
workspace.

## 44. Deviations and blockers

No contract deviation or implementation blocker remains. One initial package
check invocation was environment-related (Rtools absent from `PATH`) and was
resolved by the required toolchain-enabled rerun.

## 45. Recommended next phase

> align and simplify the STBLR and MTBLR architectures across CSR, block-eigen, and packed-BED operators, including naming, public controls, chain and trace semantics, convergence diagnostics, outputs, memory accounting, warnings, and operator reductions, before implementing extended MT-only convergence traces.

## 46. Readiness marker

PHASE 17W COMPLETE — EXTENDED MT BED CONVERGENCE CONTRACT FORMALIZED
