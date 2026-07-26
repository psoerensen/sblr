# Phase 17Q report

## 1. Executive summary

The multichain contract is implementation-ready. Phase 17Q adds only audit,
oracle, documentation, benchmark, and CI-filter material.

## 2. Repository baseline

Baseline: master at `5b31edf67fedcd48fd79da6e4adb9b57c5ab6d5a`, initially clean,
R 4.4.1 UCRT, GCC/G++/GFortran 13.2.0, `-fopenmp`, and R BLAS/LAPACK.
OpenMP reports 12 processors/maximum threads; BLAS/OpenMP thread environment
variables were unset. Hosted CI was not locally visible. The baseline source
suite passed 4,917 expectations with two established opt-in skips. Its check
exposed one committed command-output artifact as a portability warning plus two
legacy notes; the artifact was removed and the final check has no warnings.

## 3. Phase 17P verification

`mtblr_bed()` remains the public serial one-chain route through the unchanged
Phase 17O owner, view, core, finalizer, raw converter, and formatter.

## 4. Scalar multichain inventory

BayesC, BayesR, and BayesRC allocate one immutable packed owner, enumerate
trait × chain jobs into deterministic slots, use static OpenMP scheduling and
logical seeds, keep task-local mutable state/RNG, catch worker failures, and
aggregate typed results after workers. BayesC/BayesR own stability summaries;
BayesRC demonstrates optional retained chains. Their task topology is not MT's.

## 5. Phase 17O thread-safety audit

Classification: safe for independent concurrent calls after a narrow adapter
refactor. `run_mt_bed_bayesc_core()` borrows immutable binding-neutral data and
owns all mutable state and `std::mt19937`; no Rcpp, R RNG, BED I/O, file handle,
or mutable static reaches it. The current Rcpp binding combines preparation and
single execution, so preparation must be separated before dispatch.

## 6. Task topology

One zero-based chain is one all-trait task and result slot; count is `nchains`.

## 7. Shared immutable data

Owner/view, Y, maps/xx, wy, order, models, sets, priors, frequencies, and public
metadata are prepared once and shared.

## 8. Chain-private state

Beta, b, state, B/E/G/pi, R, workspace, RNG, counts, accumulators, traces,
diagnostics, timing, and failure payload are private.

## 9. Seed audit

Scalar packed routes use `seed + 1000003*(trait+1) + 9176*(chain+1)` and thus do
not preserve the supplied seed for their first trait/chain. Joint MT must not
invent a trait index.

## 10. Selected seed contract

Chain `c` uses modulo-2^32 `seed + 9176*c`; chain zero exactly preserves Phase
17O and conversion never invokes signed C++ overflow.

## 11. Single-chain reduction

One default chain resolves the current seed and executes the unchanged core
once, preserving all fields; stability is degenerate.

## 12. Internal-route decision

Phase 17R should introduce `mtblr_bed_chains_internal()` over a shared prepared
adapter also used by `mtblr_bed_internal()`, not duplicate preparation or core.

## 13. OpenMP contract

Static chain-level dispatch uses `min(ncores,nchains)` workers and chain-indexed
slots. No nested/marker/set/trait OpenMP or R API is allowed.

## 14. OpenMP-unavailable policy

Warn once on the main thread for `ncores>1`, run serial, and diagnose one worker.

## 15. BLAS contract

`ncores` governs package workers only. Global BLAS settings remain untouched;
users should generally use one BLAS thread to avoid oversubscription.

## 16. Failure contract

Workers catch standard and unknown exceptions. Failures are inspected and
reported in chain order; any failure prevents aggregation and partial output.

## 17. Timing and progress

Per-chain, mean, maximum, and dispatch seconds are diagnostics. Only the main
thread may print preparation/final summaries; worker progress is absent.

## 18. Aggregation contract

Additive posterior accumulators and retained counts are pooled in chain order;
traces are iterationwise chain means; shared prepared data is retained once.

## 19. Final-state contract

Final b/state/r/B/G/E/pi come from one-based primary chain 1. Posterior means are
pooled. Binary states are never averaged.

## 20. Chain-stability summaries

Marker posterior means yield sample SD/min/max across chains; one chain yields
zero SD and min=max=mean.

## 21. Retained-chain contract

Optional compact records retain summaries, final state/covariances, traces,
diagnostics, seed, and timing, but exclude packed/shared data, Y, wy/order,
sample residuals, genetic values, and marker r.

## 22. Raw schema decision

Optional additive fields are safe in `mtblr_raw` version 1. `chains` is stable
present-but-NULL when disabled; backend/data-level remain unchanged.

## 23. Fit-formatting decision

The sole `.as_mtblr_fit()` additively exposes stability, chain counts/seeds,
diagnostics, and present-but-NULL or named compact chains.

## 24. Future public API

Phase 17S adds `nchains=1L`, `ncores=1L`, `chain_seeds=NULL`, and
`keep_chains=FALSE` after `seed`. Phase 17Q changes no public formals.

## 25. Memory scaling

Shared bytes count genotype/Y/preparation once. Concurrent private bytes scale
with workers; typed results and optional retained output scale with chains. The
oracle reports an analytical upper bound separately from RSS.

## 26. Complexity

Preparation is once per fit; CPU work is approximately chain-linear, wall time
is worker-bounded, and aggregation is post-dispatch. No speedup claim is made.

## 27. Reproducibility matrix

The contract specifies exact single/multichain, serial/parallel, explicit/default
seed, fresh-process, intervening-fit, and completion-order reductions on one
toolchain.

## 28. Aggregation oracle

Pure-R fabricated chains verify unequal-count weighting, trace means, primary
final state, non-averaged binary state, stability, order, and compact retention.

## 29. Existing-route protection

Phase 17P/O/N/M/L/J/I/C, Phase 15A/B, all scalar BED families, schemas, BED
interface, and backend consistency remain production-unchanged.

## 30. Mutation sensitivity

The Phase 17Q audit detects the required topology, seed, ownership, threading,
failure, aggregation, schema, API, wrapper, and production-protection mutations.

## 31. Audit benchmark

Synthetic cases report shared/private/result/retained bytes and aggregation time
for requested chain/core grids; they make no OpenMP performance claim.

## 32. Installed-check behavior

Task, seed, aggregation, stability, retained-chain, and memory oracles are
portable. Only repository source assertions skip without a checkout.

## 33. Tests and CI

Phase 17Q is included in the exact fast filter. Phase 17Q alone passed 41
expectations. The protected focused tier and exact fast tier passed; the fast
tier had the established peak-RSS opt-in skip. Both full source runs passed
4,958 expectations, 0 failures, 0 warnings, and the two established opt-in
skips. Built-package tests passed the same portable tests; source assertions
skip narrowly only when source is unavailable.

## 34. Package check

The built-tarball check completed with 0 errors, 0 warnings, and 3 established
notes: the Phase 17C long fixture path, installed size (5.1 MB; libs 4.2 MB),
and legacy scalar `std::cout` symbols. No Phase 17Q note was introduced.

## 35. Diff hygiene

Wrappers, registration, NAMESPACE, man pages, and fixtures are unchanged.
Generated objects/DLL are removed after validation. The committed accidental
root command-output artifact `tatus --short`, which caused the baseline check
warning, is deleted. Diff and EOL audits pass.

## 36. Deviations and blockers

The requested `src/stblr_cpg_omp_bed_bayesrc.cpp` and
`src/blr_scheduled_execution_aggregate_impl.h` names do not exist at this
checkpoint; the actual BayesRC scheduled TU and its model-specific aggregate
header were audited. The initial compilation retry lacked Rtools on `PATH`; the
established Rtools environment compiled successfully. Neither deviation is a
contract blocker.

## 37. Recommended next phase

> implement the internal multichain and optional OpenMP chain-dispatch route defined by the Phase 17Q contract, preparing the packed-BED dataset once, sharing immutable genotype and phenotype data, executing the unchanged Phase 17O core independently per chain, aggregating deterministically, and returning optional compact chain output without changing the public `mtblr_bed()` interface.

## 38. Readiness marker

PHASE 17Q COMPLETE — MT BED MULTICHAIN AND OPENMP CONTRACT FORMALIZED
