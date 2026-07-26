# Unified BLR Framework Phase 17S report

## 1. Executive summary

Public `mtblr_bed()` multichain and optional OpenMP chain execution is active.

## 2. Repository baseline

The baseline is clean `master` at
`8aca84822f791e40a16dd6fec336608e0b8228e4`. R 4.4.1 UCRT,
GCC/G++/GFortran 13.2, OpenMP (12 reported maximum threads), R BLAS, and R
LAPACK are used. The Phase 17R baseline recorded 5,032 source expectations
with zero failures or warnings, two established opt-in skips, and a built
package check with zero errors, zero warnings, and three established notes.
Hosted CI is not locally visible.

## 3. Phase 17R verification

All Phase 17R native files, wrappers, registration, version-1 raw schema,
aggregation, seed, failure, and formatting paths remain unchanged.

## 4. Public signature

The four controls `nchains=1L`, `ncores=1L`, `chain_seeds=NULL`, and
`keep_chains=FALSE` appear immediately after `seed` and before memory controls.

## 5. Argument validation

Counts are positive integer-compatible scalars within R's integer range.
Retention is scalar logical. Explicit seeds have length `nchains`, are finite
integer-compatible signed 32-bit values, and preserve duplicates and order.

## 6. Public dispatch

Every public execution calls `mtblr_bed_chains_internal()` exactly once and
never calls `mtblr_bed_internal()`.

## 7. Default single-chain reduction

Default chain controls preserve every pre-existing numerical and metadata field;
only additive chain, stability, diagnostic, and memory fields appear.

## 8. Public seed behavior

`seed` remains the base seed. NULL explicit seeds become `integer()` for native
default resolution; explicit signed values pass unchanged in supplied order.

## 9. OpenMP behavior

Requested workers are capped by chains. Static chain-only OpenMP execution is
owned by Phase 17R and is independent of logical chain seeds.

## 10. OpenMP fallback

Without OpenMP, the native route warns once on the main thread and runs serially.
The R adapter does not duplicate the warning.

## 11. BLAS policy

The package does not mutate BLAS environment variables and records their input
snapshot. One-thread BLAS is recommended during concurrent chain execution.

## 12. Memory decomposition

Memory separates shared immutable components, private state per live worker,
typed results per chain, pooled output, and optional compact output per chain.

## 13. Memory warning

The pre-execution requested-worker upper bound drives one warning. The message
reports data/model/chain/core/worker/retention sizes and disclaims measured RSS
and measured peak RSS.

## 14. Raw enrichment

Existing marker, trait, BED, phenotype, set, and alignment metadata remain.
Data metadata adds chain/core/worker/OpenMP/seed/retention/BLAS/memory fields.

## 15. Input metadata

Input records base/requested/resolved seeds, chain policy, requested/used cores,
aggregation policies, BLAS policy/environment, and analytical memory.

## 16. Fit formatting

`.as_mtblr_fit()` remains the sole formatter. BED diagnostics, phenotype
preprocessing, and memory remain the only backend-specific additions.

## 17. Pooled posterior fields

`bm`, `dm`, `covb`, `covg`, `cove`, and `pim` pool retained samples using actual
counts.

## 18. Primary final-state fields

`b`, `d`, `r`, `vb`, `vg`, `ve`, and `pi` come from primary chain 1 and are not
averaged.

## 19. Trace semantics

`vbs`, `vgs`, and `ves` are iterationwise arithmetic means across chains.

## 20. Stability summaries

Marker stability uses sample SD/min/max across per-chain posterior means. These
are not posterior uncertainty or formal convergence diagnostics.

## 21. Retained chains

Compact records include marker summaries/state, traces, covariance summaries,
probabilities, seeds, timing, and numerical diagnostics while excluding shared
prepared data and marker/sample residuals.

## 22. Public/internal equality

Permanent tests cover one/two/four chains, serial/OpenMP, explicit/default
seeds, retention, covariance modes, updates, patterns, traits, sets, alignment,
initialization, and phenotype centering.

## 23. Pre-17S default protection

Full/diagonal, fixed/all-update, one/two/three-trait, initialization, ID, and
explicit-row default reductions use the Phase 17O serial route as oracle.

## 24. Validation failures

Zero, negative, fractional, nonfinite, missing, vector counts; malformed
retention; and wrong-length/nonfinite/fractional/out-of-range seeds fail closed.

## 25. Reproducibility

Repeated, fresh-process, serial/OpenMP, worker-count, explicit/default seed,
retained/non-retained, and intervening-fit reductions are owned permanently.

## 26. Memory tests

Tests prove one shared packed/phenotype/map/`X'Y`/order allocation, worker-private
scaling, per-chain result scaling, conditional retained scaling, and warnings.

## 27. OpenMP/BLAS tests

Worker counts follow availability and `min(ncores,nchains)`; BLAS environment
values remain byte-for-byte unchanged.

## 28. Public API protection

Only `mtblr_bed()` gains four controls. CSR, block-eigen, scalar BED, and dense
public signatures/defaults are unchanged, with one existing export.

## 29. Native protection

Phase 17R native source hashes are unchanged.

## 30. Mutation sensitivity

The Phase 17S audit guards all forty required activation, memory, semantic,
documentation, formatter, native, wrapper, namespace, API, and CI mutations.

## 31. Benchmark

All 72 small/moderate full/diagonal configurations completed for 1, 2, and 4
chains; 1, 2, and 4 requested cores; and retained-chain output disabled and
enabled. They report preparation, native dispatch, total, chain timing, worker,
analytical-memory, and object-size regression signals without adapter-overhead
or linear-speedup claims. Observed workers were capped by `nchains`.

## 32. Existing-route protection

Phase 17O/R numerics, public Phase 17P behavior, summary MT, scalar BED, schemas,
BED alignment, and backend-consistency owners remain active.

## 33. Installed-check behavior

Numerical, public, validation, memory, BLAS, formatting, and reproducibility
tests are portable. Only static source/hash/API assertions skip without source.

## 34. Generated documentation

Roxygen updates only `man/mtblr_bed.Rd`; NAMESPACE remains unchanged.

## 35. Tests and CI

The exact fast filter includes 17S and passed 2,408 expectations with zero
failures or warnings and one established opt-in skip. Focused Phase 17P/R/S
tests passed 421 expectations. The final full source suite passed 5,265
expectations with zero failures or warnings and two established opt-in skips.
Installed package tests passed.

## 36. Package check

The built-tarball check completed with zero errors and zero warnings. The same
three classified baseline notes remain: the long Phase 17C fixture path,
installed size (5.2 MiB, including 4.2 MiB under `libs`), and legacy scalar
backend `std::cout` symbols. No new note was introduced.

## 37. Diff hygiene

Final audit confirms unchanged fixtures, Phase 17R native sources, wrappers,
registration, NAMESPACE, and unrelated public APIs. Generated documentation is
limited to the existing `mtblr_bed.Rd` page; compiled and check artifacts are
removed before handoff.

## 38. Deviations and blockers

No contract deviation or blocker is currently known. Timing remains
intentionally nondeterministic.

## 39. Recommended next phase

> audit and formalize convergence-diagnostic requirements for public MT BED multichain fits, including split-R-hat, effective sample size, Monte Carlo standard errors, trace-retention requirements, memory scaling, marker-level diagnostic scope, warning thresholds, and behavior when compact chains are not retained, without yet implementing convergence diagnostics.

## 40. Readiness marker

PHASE 17S COMPLETE — PUBLIC MT BED MULTICHAIN AND OPENMP EXECUTION ACTIVE
