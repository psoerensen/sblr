# Unified BLR Framework Phase 17R Report

## 1. Executive summary

The internal MT BED multichain route is active as `mtblr_bed_chains_internal()`.
It dispatches independent unchanged Phase 17O chains over one prepared dataset.

## 2. Repository baseline

The baseline is master at `0efe4a8b267690bf001bcd9a43f786bf56c6b053`
with a clean initial tree. R 4.4.1, GCC/G++/GFortran 13.2.0, OpenMP (maximum
12 threads), R BLAS, and R LAPACK were used. The Phase 17Q baseline source
suite had 4,958 passing expectations. Baseline source and installed checks
passed with the three established package notes and no error or warning.

## 3. Phase 17Q verification

The chain-only topology, +9176 seeds, ownership, failure ordering, pooling,
trace means, primary state, stability, compact output, and schema decisions map
directly to production types, tests, audits, raw fields, and formatter fields.

| Contract | Production symbol | Test/audit owner | Raw/fit field | Documentation |
|---|---|---|---|---|
| one shared preparation | `prepare_mt_bed_adapter()` | Phase 17O/17R exactness; chains audit | shared `wy`/order | sections 4-5 |
| chain tasks and seeds | `MtBedChainTask`, `resolve_mt_bed_chain_seed()` | seed and worker tests; mutation audit | `chain_seeds` | sections 6-8 |
| private execution | `run_mt_bed_chain_task()` | serial/OpenMP tests; chains audit | per-chain diagnostics | sections 9-13 |
| deterministic aggregate | `aggregate_mt_bed_chains()` | pooling/stability oracle | pooled marker/variance/pi | sections 14-18 |
| compact chains | `MtBedChainSummary` and R converter | retained-chain tests | `chains`; `fit$chains` | sections 19, 21-23 |
| one final path | default finalizer/legacy/raw calls | single-chain and architecture tests | version-1 raw | sections 20-23 |

## 4. Shared prepared adapter

`MtBedPreparedAdapter` owns one stationary `PackedBedMatrix` through a
`unique_ptr`, maps, centered Y, marker X'Y, stable order, and trait-major X'Y.
Its method creates a fresh immutable view after final owner placement.

## 5. Serial-route refactor

`mtblr_bed_internal()` calls shared preparation, then one unchanged core,
finalizer, legacy adapter, and raw converter. Serial output is unchanged.

## 6. Chain task and result types

Tasks carry zero-based chain and uint32 seed. Initialized results carry chain,
seed, failure text, elapsed seconds, and the existing `MtBedCoreResult`.

## 7. Seed resolution

Chain `c` uses `(uint32(seed) + 9176*c) mod 2^32`; explicit signed ints convert
modulo 2^32 in supplied order.

## 8. Single-chain reduction

Chain zero preserves Phase 17O. Full/diagonal and fixed/all-update reductions
match all ordinary serial numerical fields exactly.

## 9. Chain-private ownership

Each core owns beta, b, state, B, E, G, pi, residual, workspace, RNG, counts,
accumulators, traces, and numerical diagnostics.

## 10. OpenMP dispatch

The sole new OpenMP loop is a static chain-level parallel-for with one result
slot per chain and `min(ncores,nchains)` workers.

## 11. OpenMP fallback

Without OpenMP, requests above one core warn once, use one worker, and retain
serial numerical behavior.

## 12. BLAS policy

`ncores` governs package workers only. No BLAS environment is changed; BLAS
should normally use one thread during chain parallelism.

## 13. Exception capture

Every task catches standard and unknown C++ exceptions inside the worker.

## 14. Aggregation

Marker, covariance, and pi accumulators and retained counts are summed in chain
order and finalized once.

## 15. Trace aggregation

`vbs`, `vgs`, and `ves` are added in chain order and divided by `nchains` once.

## 16. Primary-chain final state

Final r, b, state, B, G, E, pi, order, and legacy diagnostics come from chain 1.

## 17. Aggregate diagnostics

Jitter attempts and E updates are summed, maximum jitter is maximized, and all
values are retained per chain.

## 18. Stability summaries

Marker posterior means yield sample SD, min, and max for bm and dm. One-chain
SD is zero and min=max=the pooled mean.

## 19. Compact retained chains

Optional records contain statistical states, summaries, traces, covariances,
pi, diagnostics, and seconds; they exclude r and all shared prepared objects.

## 20. Internal signature

All serial arguments are followed by `method`, `nchains`, `ncores`, integer
`chain_seeds`, and `keep_chains`; the route is registered but not exported.

## 21. Raw schema

The extension remains `mtblr_raw` version 1, backend `mt_bed_bayesc`, data level
`individual`, with additive meta, stability, diagnostics, and `chains` fields.

## 22. Raw validation

Serial raw follows the old branch. Extended raw validates counts, dimensions,
seeds, timing, policies, diagnostics, stability, and compact exclusions.

## 23. Fit-formatting support

The sole `.as_mtblr_fit()` names optional stability and compact records and
exposes chain count, seeds, diagnostics, and present-but-NULL chains.

## 24. Public-route protection

Public `mtblr_bed()` retains exact formals and calls only the serial route.

## 25. Serial multichain tests

Tests cover one to four chains, both E modes, update controls, sets, traits,
patterns, initial states, and both retention modes.

## 26. Serial/OpenMP tests

Worker-count reductions compare pooled and retained fields after timing and
worker diagnostics are normalized.

## 27. Aggregation oracle

Independent chain means verify pooled counts, trace means, primary finals,
sample SD/min/max, and stable order.

## 28. Failure tests

Malformed states list all failing chains in ascending one-based order under
serial and parallel dispatch, with no partial output.

## 29. Retained-chain tests

Names, indices, seeds, dimensions, diagnostics, and exclusions pass; disabled
retention is stable NULL.

## 30. Reproducibility

Repeated, worker-count, default/explicit seed, full/diagonal, and existing
fresh-process reductions are deterministic on this toolchain.

## 31. Ownership and threading audit

The audit reports one owner/view/read/preparation, chain-indexed tasks/results,
one static loop, no R/Rcpp/printing in workers, and one finalization path.

## 32. Mutation sensitivity

All 38 required ownership, threading, seed, failure, aggregation, output,
public-route, and protected-core mutations are detected.

## 33. Benchmark

Small full/diagonal runs report chain, dispatch, total, and object-size signals
without a speedup claim.

## 34. Existing-route protection

Phase 17O/P, summary MT, scalar BED, schema, BED, and backend owners pass.

## 35. Installed-check behavior

Numerical chain tests are portable; source assertions skip only without source.

## 36. Generated wrappers

Attributes add exactly one internal wrapper and one native registration.

## 37. Tests and CI

The fast workflow includes 17R. Its exact local filter passed 2,175
expectations with one established opt-in skip. The final full source suite
passed 5,032 expectations with zero failures, zero warnings, and two
established opt-in skips. Phase 17R itself owns 74 passing expectations.

## 38. Package check

The built-tarball check passed its installed tests with zero errors and zero
warnings. Its three notes are unchanged from baseline: the long Phase 17C
fixture path, installed size (5.2 MiB, including 4.2 MiB under `libs`), and
legacy scalar-backend `std::cout` symbols. No Phase 17R note was introduced.

## 39. Diff hygiene

Final audit covers wrappers, namespace, fixtures, EOL, and artifacts.

## 40. Deviations and blockers

No blocker remains. Timing is intentionally nondeterministic.

## 41. Recommended next phase

> activate the validated Phase 17R multichain route in public `mtblr_bed()` by adding `nchains`, `ncores`, `chain_seeds`, and `keep_chains`, extending the analytical memory warning and public metadata, documenting BLAS oversubscription and retained-chain output, and preserving the current default single-chain result exactly.

## 42. Readiness marker

PHASE 17R COMPLETE — INTERNAL MT BED MULTICHAIN AND OPENMP EXECUTION ACTIVE
