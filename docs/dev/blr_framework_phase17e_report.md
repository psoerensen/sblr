# BLR Framework Phase 17E Report

## 1. Executive summary

The corrected authoritative public `mtblr()` execution now uses six explicit
typed, binding-neutral production contracts and one callable
`run_mt_default_core()`. Phase 17C numerical behavior, the legacy 20-position
schema, and public R formatting are unchanged. Native aggregation, a named
converter, a generic operator, and alternative-backend migration were not
started.

## 2. Repository baseline

Phase 17E began from a clean `master` at `31c4973` (`Extract corrected public
multivariate execution`), the committed Phase 17D baseline. No combined hosted
Actions status was visible locally, so hosted CI success is not claimed. The
Windows Rtools44 C++17 toolchain compiled and loaded the baseline. Phase 17D's
full suite had 6,120 passes, zero failures, zero warnings, and eleven opt-in
skips. Its new local matched benchmark recorded updated/fixed/explicit/moderate
means of 0.018/0.008/0.010/0.026 seconds and completed-fit RSS of
133.13/133.14/120.91/123.06 MiB.

## 3. Public call graph

`sblr(algorithm="default")` retains its R validation and dense-summary
preparation, calls unchanged `mtblr()` through the generated wrapper, and the
native adapter constructs typed views/specifications/state. It calls
`run_mt_default_core()` once, reads `MtDefaultCoreResult` in the unchanged
positional finalization, and returns the same 20 positions. R naming,
orientation, and conditional correlations remain in `R/interface_mtblr.R`.

## 4. Phase 17D lexical boundary

Phase 17D included a 240-line body inside `mtblr()` and obtained all inputs and
outputs lexically. Phase 17E moves that include to translation-unit scope,
wraps the same numerical body in an explicit function, qualifies execution
controls, and transfers its outputs into an owning result. No numerical helper
was moved or renamed.

## 5. Dependency inventory

| Symbol(s) | Current type/owner | Mode/lifetime/size | Typed contract | Finalization |
|---|---|---|---|---|
| `wy`, `ww` | nested doubles, native by value | borrowed read-only for fit; marker/trait scale | data view | `wy` yes |
| `yy`, `n` | vectors, native by value | borrowed read-only for fit; trait scale | data view | no |
| `XXvalues`, `XXindices` | nested sparse summaries, caller-owned refs | borrowed read-only; dominant `O(nt*m^2)` representation | data view | no |
| `models`, `sets` | nested integers | borrowed read-only; model/set scale | model spec | no |
| `method` | integer | copied scalar | model spec | no |
| `ssb_prior`, `sse_prior` | nested doubles, native by value | borrowed read-only; trait covariance scale | covariance-prior view | no |
| `nub`, `nue` | doubles | copied scalar | covariance-prior view | no |
| update flags, `nit/nburn/nthin/seed` | booleans/integers | copied fit controls | execution spec | trace sizes use unchanged adapter controls |
| initial `b/B/E/pi` | nested vector/Armadillo/vector, native by value | moved, mutable, core-owned | initial state | final values through result |
| sampling helpers | translation-unit functions | called only by core for fit | existing helpers | no |
| dimensions and five counts | core scalars | core/result lifetime | core result | yes |
| `bm/dm/r/b/d/order` | nested vectors | core workspace then result-owned | core result | yes |
| `vbs/vgs/ves` | nested trace vectors | core workspace then result-owned | core result | yes |
| `cvbm/cvgm/cvem` | nested covariance accumulators | core workspace then result-owned | core result | yes |
| `B/G/E`, `pi/pis` | Armadillo/vector state | core state then result-owned | core result | yes |
| `pistrait/pismarker` | legacy accumulator vectors | core then result-owned | core result | yes |
| beta, probabilities, scores, inverses, temporary states | existing vector/Armadillo workspaces | core-local only | none | no |

## 6. Typed contracts

`MtDefaultDataView` borrows `wy`, `ww`, `yy`, `XXvalues`, `XXindices`, and `n`.
`MtDefaultModelSpec` borrows ordered `models` and zero-based `sets` and copies
`method`. `MtDefaultCovariancePriorView` borrows both prior scales and copies
degrees of freedom. `MtDefaultExecutionSpec` copies three update flags,
iteration controls, and resolved seed. `MtDefaultInitialState` owns mutable
`b/B/E/pi`. `MtDefaultCoreResult` owns dimensions, five counts, summaries,
state, traces, covariances, ordering, final parameters, and legacy probability
accumulators. None contains Rcpp, SEXP, R names, classes, or schema metadata.

## 7. Ownership model

R and the generated wrapper own inputs during conversion/call. The three view
contracts borrow immutable native objects for the complete core call. Execution
controls are small copied values. Existing by-value `b/B/E/pi` parameters move
into initial state, then into core locals; the moved-from adapter parameters are
never read. Workspaces and `std::mt19937` are core-local. Final numerical state
moves into the owning core result. The legacy finalizer reads that result and
the still-live borrowed `wy`; R formatting retains binding ownership.

## 8. Callable core signature

```cpp
MtDefaultCoreResult run_mt_default_core(
    const MtDefaultDataView& data,
    const MtDefaultModelSpec& model,
    const MtDefaultCovariancePriorView& prior,
    const MtDefaultExecutionSpec& execution,
    MtDefaultInitialState initial_state);
```

It is inline, namespace-scoped, binding-neutral, and called exactly once.

## 9. Files changed

- Added `src/blr_mt_default_types.h` and converted
  `src/blr_mt_default_core_impl.h` to the callable core.
- Thinned `src/mtblr.cpp` to contract construction, one call, aliases, and
  positional finalization.
- Added Phase 17E structural/reference/fresh tests and updated Phase 17B-D
  structural assertions for the superseding architecture.
- Added the Phase 17E benchmark/copy audit and CI coverage.
- Updated the implementation plan, capability/reduction matrices, computation
  inventory, naming vocabulary, and this report.

Production R, generated wrappers, fixtures, and protected backends did not
change.

## 10. Numerical operation-order preservation

```text
PHASE17D_NUMERICAL_STATEMENTS=92
PHASE17E_NUMERICAL_STATEMENTS=92
ORDER_EQUIVALENT=TRUE
RNG_CALL_ORDER_EQUIVALENT=TRUE
UPDATE_CALL_ORDER_EQUIVALENT=TRUE
```

Both bodies contain the same 240 numerical lines. Normalizing only
`execution.*` qualification restores Phase 17D MD5
`7e8ea9e4812ce57a701416f8896a97cc`. Nonqualification changes are the function
signature, input aliases, initial-state moves, and result moves after the
retained-sample guard; none changes a numerical statement.

## 11. Correctness-contract preservation

Set-local B updates remain jointly guarded by `execution.updateB`; the global
B update retains its method-and-update guard. Post-burn remains
`it >= execution.nburn`; marker thinning remains relative to burn-in. All five
counts and accumulator-specific denominators remain distinct. Disabled B/E/pi
summaries remain zero, fixed input states remain exact, and the retained-marker
guard prevents zero division.

## 12. RNG ownership

One fit-local `std::mt19937 gen(execution.seed)` remains at the same point after
initialization and before the one Gibbs loop. Draw sites and helper order are
unchanged. The core has no R RNG, Armadillo RNG, worker identity, static state,
or thread-local state.

## 13. Data representation

Summary inputs retain their nested-vector representation, including sparse
per-marker `XXvalues/XXindices`. The typed view adds no dense or sparse
conversion. Armadillo remains limited to dense trait-level B/G/E algebra and
inverses. No Eigen or custom operator was introduced.

## 14. Core result

The result owns `nt/m/nmodels`; marker/B/G/E/pi counts; `bm/dm/r/b/d/order`;
`vbs/vgs/ves`; `cvbm/cvgm/cvem`; final `B/G/E`; final and accumulated
`pi/pis`; and legacy `pistrait/pismarker`. Immutable `wy` is deliberately not
duplicated: the finalizer reads it through the still-live native parameter.

## 15. Legacy finalization

Allocation, resizing, division, copying, and positions 1-20 remain inline in
`src/mtblr.cpp`, adapted only to read `core_result`. Positions 19-20 keep their
zero/unsupported behavior. No native aggregation or conversion object exists.

## 16. Public R interface

`R/interface_mtblr.R`, both generated wrapper files, `NAMESPACE`, native
signature, public arguments/defaults, algorithm names, hidden seed, `.Call`
target, names, dimensions, orientations, states, conditional correlations, and
omission behavior remain byte-identical.

## 17. Phase 17B historical fixtures

All three files remain immutable at SHA-256
`82AF2F...874C9`, `CF32B7...869C`, and `589E3A...13D9`. They remain historical
defect/denominator evidence, not current expected output.

## 18. Phase 17C corrected references

Corrected raw references pass 3/3 and formatted references pass 3/3 with exact
structure and numerical tolerance `1e-12`. Their SHA-256 values remain
`2AD954...C806`, `4F1929...F3A`, and `8F7BD2...3001`.

## 19. Reproducibility

`A;A`, `A;B;A`, OMP environment 1 versus 2, intervening canonical CSR,
intervening packed-BED, and fresh versus reused-process comparisons remain
exact in structure and equal within `1e-12`. Multiple chains remain unsupported.

## 20. Scientific identities

Finite outputs, dimensions/names/order, binary states, bounded `dm`, symmetric
positive-semidefinite covariance matrices, trace lengths/orientation, fixed
B/E/pi, normalized updated `pim`, retained counts, finite covariance means,
disabled zero summaries, and covariance-diagonal trace identities remain valid.

## 21. Alternative multivariate protection

`mt_cpg.cpp`, `mt_cpg_arma.cpp`, `mt_cpg_omp.cpp`, and `mt_cpg_omp_csr.cpp`
remain byte-identical. Protected `mtblr_hybrid()` and `mtblr_eigen()` regions
are unchanged and neither calls/includes the new core. The worker-sensitive
OpenMP risk evidence remains active.

## 22. Scalar, packed-BED, and block-eigen protection

Canonical scalar CSR, canonical packed-BED, shared packed-BED, and block-eigen
source/reference protections pass. Generated wrappers and `NAMESPACE` retain
their baseline hashes.

## 23. Performance, memory, ownership, and I/O

Phase 17E timings were: updated `0.39,0.03,0.04,0.02,0.01` seconds (mean 0.098,
median 0.03, RSS 120.65 MiB); fixed B `0.03,0.01,0.03,0.01,0.01` (0.018, 0.01,
120.71); explicit sets `0.02,0.02,0.02,0.02,0.03` (0.022, 0.02, 119.53); and
moderate dense `0.14,0.18,0.22,0.22,0.24` (0.200, 0.22, 120.69). Short Windows
debug timings are noisy regression signals, not speed claims. Completed-fit RSS
is not peak RSS; peak RSS was not sampled.

Dense XX remains `O(nt*m^2)` and is borrowed by const reference with no typed
copy. Existing wrapper conversions are unchanged. Initial state and result
transfers use moves; final positional allocation is unchanged. No unexplained
`O(nt*m^2)` copy and no MCMC-time I/O were introduced.

## 24. CI coverage

The fast workflow includes ordinary Phase 17E tests. Extended CI adds
`SBLR_RUN_PHASE17E_FRESH: "true"` without removing prior variables. Benchmarks
are not in the fast gate. Hosted CI success is not claimed.

## 25. Tests

- Phase 17B historical evidence and fixture integrity: passed.
- Phase 17C corrected references and contracts: passed.
- Phase 17D extraction/supersession protections: passed.
- Phase 17E ordinary: 160 passed and one opt-in fresh skip.
- Phase 17E fresh enabled: 161 passed, zero failures.
- Raw/formatted corrected references: 3/3 and 3/3.
- Alternative and canonical protections: passed.
- Full suite: 6,282 passed, zero failures, zero warnings, and twelve opt-in skips.

## 26. Deviations and blockers

No correctness or ownership blocker remains. Debug-build timings varied
materially between short runs, including before Phase 17E, while completed-fit
RSS remained stable and operation-order/reference tests were exact. The typed
result necessarily owns the state required by finalization; moves avoid full
duplicate transfers. The route remains legacy and noncanonical until
finalization/aggregation and binding conversion are separated.

## 27. Recommended next phase

Move the legacy posterior finalization into one typed native
aggregate/result-finalization function and leave a thin positional schema
adapter in `mtblr()`, while preserving the Phase 17C corrected references,
public 20-position output, R formatting, and the typed callable-core contract.

## 28. Readiness marker

PHASE 17E COMPLETE — TYPED PUBLIC MULTIVARIATE CORE BOUNDARY ACTIVE
