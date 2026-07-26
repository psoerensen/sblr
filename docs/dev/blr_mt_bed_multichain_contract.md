# Multichain individual-level MT packed-BED contract

## 1. Purpose

This contract specifies future multichain and chain-parallel MT BED execution.
Phase 17Q is audit-only; Phase 17R will implement the internal route.

## 2. Current status

Phase 17P exposes one serial chain. `mtblr_bed()` and
`mtblr_bed_internal()` have no multichain or OpenMP controls today.

## 3. Scalar packed-BED precedents

| backend | task | owner/view | dispatch | seed | aggregation |
|---|---|---|---|---|---|
| BayesC | trait × chain | one `FastPackedBedMatrix`, common view | static OpenMP | trait and chain offsets | typed slots; mean/stability fields |
| BayesR | trait × chain | one `FastPackedBedMatrixBR`, common view | static OpenMP | same helper | typed slots; mean/stability fields |
| BayesRC | trait × chain | one `FastPackedBedMatrixBR`, ownership-equivalent view | static OpenMP | same helper | typed slots; optional retained chains |

All readers finish before worker dispatch. Mutable residuals, states, workspaces,
and RNGs are task-owned. R objects are assembled after dispatch.

| inventory item | scalar BayesC | scalar BayesR | scalar BayesRC | current MT BED |
|---|---|---|---|---|
| task/result index | `chain*nt+trait` via `BedFamilyTaskIndex` | same | same | one serial all-trait task |
| seed | common trait/chain offset helper; uint64 then engine cast | same | same | supplied seed cast by fit-local engine |
| immutable context | common view, maps, order, phenotype, priors | common view plus components/scheduler | ownership-equivalent view plus annotations | common view, maps, Y, wy, order, models/sets/priors |
| mutable state | effects/state/residual/workspace/RNG per task | component/effect/state/residual/workspace/RNG per task | annotation/component/effect/state/residual/workspace/RNG per task | joint beta/b/state/B/E/G/pi/R/workspace/RNG |
| OpenMP | static trait × chain loop | static trait × chain loop | static trait × chain loop | none |
| exception | typed `failed/error`, checked in job order | typed failure plus deferred progress, checked in job order | typed failures, first aggregate failure | ordinary serial exception propagation |
| timing/progress | per job; main-thread summary, some historical worker output | per job; events replayed on main thread | per job; post-worker binding | MT diagnostics only; no live worker concept |
| aggregation/final state | arithmetic chain aggregation, including averaged scalar final state | arithmetic chain aggregation | arithmetic chain aggregation | no chain aggregation |
| stability | sample SD/min/max | sample SD/min/max | model-specific aggregate fields | none |
| retained chains | no | no | optional complete typed task records | none |
| raw/fit | scalar raw converter/formatter | scalar raw converter/formatter | scalar raw with optional chains | shared MT raw/finalizer/formatter |
| memory | one owner plus task-private state/results | same | same plus optional retained tasks | one owner plus one joint state |

The scalar final-state averaging precedent is explicitly rejected for MT because
joint binary states and correlated B/E/G are coherent states, not estimands to
average fieldwise.

## 4. Critical MT distinction

One MT task is one complete joint chain containing all traits. A trait × chain,
marker × chain, or set × chain topology is invalid.

## 5. Task topology

`task_count = nchains`. The task index and deterministic result slot are the
zero-based chain index. Completion order and worker identity never select slots.
A future binding-neutral `MtBedChainTask` contains `chain` and resolved seed.

## 6. Shared immutable data

One prepared fit shares the `PackedBedMatrix`, `BedPackedGenotypeView`, centered Y,
marker maps and `xx`, `X'Y`, marker order, model patterns, sets, covariance
priors, allele frequencies, and BED/sample/marker metadata. All are immutable
before dispatch.

## 7. Chain-private state

Each chain owns latent beta, effective b, binary state, B, E, G, pi, sample
residual R, decoded-marker workspace, one RNG, model counts, posterior accumulators,
traces, jitter/E diagnostics, timing, and failure payload. Initial values are
copied identically to every chain.

Thus every chain owns one sample residual R and no mutable state is shared.

## 8. Owner and view lifetime

Exactly one stationary packed owner and one immutable borrowed view outlive all
workers. No chain rereads BED, copies packed bytes, or holds a file handle.

## 9. Seed policy

For zero-based chain `c`, default resolution is
`uint32(seed + 9176*c)`, evaluated without signed overflow and reduced modulo
2^32. Thus chain zero preserves Phase 17O exactly. Explicit seeds have length
`nchains`, remain in supplied order, are finite integer-compatible values, and
are converted by the same modulo-2^32 rule. Workers never affect seeds.

## 10. Single-chain reduction

With `nchains=1` and default seeds, the resolved seed equals the current seed,
one unchanged Phase 17O core executes, and every existing numerical field is
exact. Additive metadata and degenerate stability summaries are permitted.

The audited seam is `run_mt_bed_bayesc_core()` in
`src/blr_mt_bed_core_impl.h`, called by `mtblr_bed_internal()` in
`src/mtblr.cpp`. `MtBedDataView` borrows immutable
objects; `MtBedInitialState`, the sample residual, decoded workspace,
accumulators, diagnostics, and local `std::mt19937` are execution-owned.
`sampleBset`, `sampleB_latent`, `sampleB`, `samplePi`, and covariance draws all
receive that engine. No R RNG, Rcpp object, file operation, or mutable static is
reachable from the core. This establishes the narrow-refactor classification.

## 11. Future internal route

Phase 17R should extract a binding-neutral `prepare_mt_bed_adapter()` used by
both current `mtblr_bed_internal()` and a new
`mtblr_bed_chains_internal(..., int nchains, int ncores,
std::vector<int> chain_seeds, bool keep_chains)`. Empty native seeds mean use
the default resolver. Preparation occurs once; the numerical core is unchanged.

## 12. OpenMP dispatch

Dispatch is chain-level only:

```cpp
#pragma omp parallel for schedule(static) num_threads(worker_count)
for (int chain = 0; chain < nchains; ++chain)
  results[chain] = run_one_mt_bed_chain(...);
```

`worker_count=min(ncores,nchains)`. There is no nested, marker, set, or trait
parallelism, worker seed, R API, Rcpp allocation, printing, or callback.
There is no R API use inside workers.

## 13. OpenMP-unavailable behavior

If `ncores>1` is requested without OpenMP, warn once and run serial on the main
thread with `used_workers=1` and `openmp_available=FALSE`. Results remain exact.

## 14. BLAS interaction

`ncores` controls package OpenMP workers only. The package does not change BLAS
settings. Users should normally configure BLAS to one thread during chain
parallelism and inspect `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
`VECLIB_MAXIMUM_THREADS`, `OMP_NUM_THREADS`, and `OMP_THREAD_LIMIT`.
This explicit BLAS oversubscription policy is diagnostic guidance, not a global
environment mutation.

## 15. Exceptions

Every initialized result slot contains `failed`, `error`, and `chain`. Workers
catch `std::exception` and unknown exceptions. After joining, slots are checked
in ascending chain order. The first chain is the primary error and all failures
are included in diagnostic text. There is no partial aggregation or OpenMP
cancellation.

## 16. Timing

Diagnostics contain per-chain seconds, their mean and maximum, and dispatch
elapsed seconds. Timing is nondeterministic and tests normalize it.

## 17. Progress

Workers emit no progress. `verbose` may report preparation and final summaries
on the main thread. Native parallel execution has limited R-interrupt
responsiveness; Phase 17R adds no callbacks or event queue.

## 18. Aggregation

In ascending chain order, sum bm, dm, covB, covG, covE, pi accumulators and their
retained counts, then call the existing finalizer once. This weights unequal
retained counts correctly. `vbs`, `vgs`, and `ves` use an iterationwise arithmetic
mean across chains. `wy`, order, models, and sets remain shared.

## 19. Final-state policy

`primary_chain = 1`, `final_state_policy = "primary_chain"`,
`posterior_summary_policy = "pooled_retained_samples"`, and
`trace_policy = "iterationwise_chain_mean"`. Final b, state, r, B, G, E, and pi
come from chain 1. Binary states are never averaged.

## 20. Chain-stability summaries

`bm_sd`, `bm_min`, `bm_max`, `dm_sd`, `dm_min`, and `dm_max` summarize per-chain
posterior marker means in marker-by-trait orientation. SD is the sample SD with
denominator `nchains-1`. For one chain SD is zero and min=max=the pooled mean.

## 21. Retained chains

`keep_chains=FALSE` retains only pooled output, stability, seeds, and diagnostics.
When true, `chain1`, `chain2`, ... contain index, seed, bm, dm, final b/state,
traces, covariances, final B/G/E, final/mean pi, jitter/E diagnostics, and seconds.
They omit packed bytes, Y, maps, wy, order, models, sets, sample R, genetic U,
and final marker r.

## 22. Raw schema

The additive extension remains `mtblr_raw version 1`, backend `mt_bed_bayesc`,
data level `individual`. Meta adds `nchains` and `keep_chains`; marker adds the
six optional stability fields; MT BED diagnostics add cores, workers, OpenMP,
seeds, timing, and aggregation policies. `chains` is a stable present-but-`NULL`
field when not retained and a compact named list otherwise.

## 23. Fit formatting

`.as_mtblr_fit()` remains the only general formatter. It additively exposes
stability fields plus `fit$nchains`, `fit$chain_seeds`,
`fit$chain_diagnostics`, and a present-but-`NULL` `fit$chains` when disabled.

## 24. Future public API

Phase 17S adds, after `seed=1`, `nchains=1L`, `ncores=1L`,
`chain_seeds=NULL`, and `keep_chains=FALSE`, before memory and verbose controls.
Counts are positive integers, seeds are NULL or length `nchains`, and
`keep_chains` is scalar logical. Phase 17Q changes no formals.

## 25. Memory

Report `shared_bytes`, `private_state_bytes_per_chain`,
`result_bytes_per_chain`, `retained_chain_bytes_per_chain`, `worker_count`,
`nchains`, `estimated_concurrent_bytes`, `estimated_retained_output_bytes`,
`estimated_total_bytes`, and GiB. The packed owner is counted once;
concurrent private state scales with workers, typed results with chains, and
optional public records with chains. This is an analytical upper bound, not RSS
or measured peak RSS.

## 26. Complexity

BED, maps, wy, and order cost occurs once. Total Gibbs CPU work is approximately
linear in `nchains`; concurrent wall time is bounded by workers. Aggregation and
stored results scale with chains. No linear-speedup claim is made.

## 27. Reproducibility

Permanent reductions cover single-chain Phase 17O equality, serial/parallel
equality, repeated and fresh-process runs, explicit/default seeds, worker-count
changes, intervening scalar and summary fits, and completion-order perturbation.
Different seeds must alter a stochastic field on the same toolchain.

## 28. Failure semantics

Failure messages and chain order are deterministic. If any task fails, no pooled
or retained result is returned and no successful partial result is finalized.

## 29. Unsupported scope

Marker/set/trait parallelism, distributed execution, checkpoints, asynchronous
progress, per-chain initials, and cross-toolchain bitwise guarantees are absent.

## 30. Phase 17R implementation plan

Extract shared preparation, define chain task/result bundles, resolve seeds,
preallocate slots, dispatch static chain workers, capture failures, aggregate
accumulators/traces/states, add optional raw fields, and prove all serial and
parallel reductions while leaving the Phase 17O core unchanged.

## 31. Phase 17S public activation plan

Add the four public controls, validate them in R, extend memory warnings and
metadata, dispatch the Phase 17R route, format optional chains through the one
formatter, document examples, and preserve `nchains=1` exactly.

## Phase 17R implementation status

The contract is now implemented internally. `mtblr_bed_chains_internal()`
prepares one stationary owner/view, dispatches chain-only tasks through one
static OpenMP loop when available, captures all errors, pools actual counts,
averages traces, uses chain 1 final state, computes sample-SD stability, and
optionally returns compact chains. Public activation remains Phase 17S.

## Phase 17S activation status

The four chain controls are now public. Every `mtblr_bed()` call executes the
chains route once, including the default one-chain reduction. Public metadata
reports requested/used workers, resolved seeds, aggregation policies, BLAS
environment, stability, optional compact records, and multichain memory.

## Phase 17T diagnostic trace ownership

Formal diagnostics require a separate typed post-burn trace bundle constructed
before aggregation/discard. It is independent of compact-chain retention and
initially contains only Tier 1 B/G/E diagonals.
## Phase 17U diagnostic trace ownership

Tier 1 convergence capture is independent of `keep_chains`. One immutable
prepared dataset still feeds one core call per chain; post-burn diagonal traces
are copied before aggregation, then diagnosed in dependency-free R code.
Compact chains, aggregation, final-state policy, seeds, and OpenMP dispatch are
unchanged.

## Phase 17V public selection

The public adapter selects exactly one ordinary or trace-capable native route.
It never reruns chains and never uses pooled or compact-chain data as formal
diagnostic draws.
