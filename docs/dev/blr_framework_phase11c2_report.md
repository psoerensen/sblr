# Unified BLR Framework Phase 11C2 report

## 1. Executive summary

The corrected public scheduled packed-BED BayesC per-chain execution now runs
behind a typed binding-neutral context and returns a typed per-chain result.
Task dispatch, cross-chain aggregation and R conversion remain in the adapter.

## 2. Repository baseline

The clean baseline was `master` at Phase 11C1 commit `a1c1a99`. Initial status
and `git diff --check` were clean. R 4.4.1 with Rtools44/GCC and OpenMP compiled
the baseline. Its full suite passed 4,747 expectations with two opt-in skips;
the Phase 11A/11B/11C1 focused run passed 124 expectations.

## 3. Lexical dependency inventory

| Symbols | Type/dimensions | Classification, owner, lifetime | Typed destination |
|---|---|---|---|
| `G` packed bytes, `n/m/nbytes/stride` | fit-sized packed marker-major storage | borrowed immutable; adapter fit owns storage through all tasks | `BedPackedGenotypeView` |
| `marker_maps`, `marker_order` | marker-length vectors | borrowed immutable marker metadata; fit lifetime | context references |
| `y_mat` | samples x traits | borrowed immutable phenotype; adapter fit lifetime | `phenotype` |
| `b_init` | traits x markers | borrowed immutable initial effects; adapter fit lifetime | `initial_effects` |
| `B`, `E`, `ssb_prior_mat`, `sse_prior_mat` | traits x traits | borrowed immutable statistical inputs; adapter fit lifetime | typed matrix references |
| `pi` | length two | borrowed immutable model input | `initial_pi` |
| `nub`, `nue`, `adjE`, pi-prior values and update flags | scalars | copied immutable model controls; chain-call lifetime | context scalar fields |
| `nit`, `nburn`, `nthin`, `rebuild_every`, `progress_every` | scalars | copied execution controls | context scalar fields |
| full-sweep controls | frequency/convention | copied scheduler control; identical CSR/BED meaning | `ScheduledSweepControl` |
| null-skip controls | base/max/burn-in/growth | copied scheduler control; identical CSR/BED meaning | `NullSkipControl` |
| candidate controls | threshold/lifetime | copied scheduler control; identical CSR/BED meaning | `CandidateControl` |
| resolved seed, trait and chain indices | scalar | adapter-resolved logical-chain metadata | context fields |
| effects, residuals, indicators, variances and pi | marker/sample/scalar | mutable chain-owned sampler state | core locals |
| active/candidate/due vectors and buckets | marker/iteration sized | mutable chain-owned scheduler state | core locals |
| `BedScheduledBayesCChainRng` | one engine/normal/uniform | mutable logical-chain-owned; one call lifetime | core local |
| posterior accumulators, CPO, timing and failure | chain sized | mutable chain-owned then returned | typed result |
| job index, `job_results`, OpenMP loop, aggregate matrices | task/cross-chain | adapter-only; excluded from context/core | retained adapter storage |
| R names, schema and list construction | binding metadata | adapter/converter-only; excluded | retained adapter conversion |

## 4. Files changed

The selected adapter and implementation header activate the boundary. New
`blr_bed_scheduled_bayesc_types.h` defines the view, context, validation and
result. Phase 11C2 tests, benchmark and this report were added. The plan and
capability matrix were updated. Obsolete protected-source hashes and the Phase
11C1 lexical function-name assertion were replaced with typed-path protection.

## 5. Packed genotype view

`BedPackedGenotypeView<PackedGenotype>` holds a const reference to fit-owned
decoded storage plus a const byte pointer and packed-size, marker-count,
sample-count, bytes-per-marker and stride metadata. It owns nothing, contains no
path or file handle, and must not outlive the adapter's packed matrix.

## 6. Typed per-chain execution context

`BedScheduledBayesCChainExecutionContext<PackedGenotype>` borrows genotype,
marker maps/order, phenotype, initial effects, priors and initial parameters;
copies MCMC/model/update/scheduler controls and the already resolved chain seed;
and carries logical trait/chain indices. It contains no R/Python types, names,
schema metadata, file state, worker identity or mutable cross-chain output.
Field-specific validation covers packed dimensions, phenotype/initial/prior
dimensions, iteration controls, scheduler controls, indices and prior values.

## 7. Scheduler-contract reuse

`ScheduledSweepControl`, `NullSkipControl` and `CandidateControl` are reused.
Phase 11A established identical initialization, update/reset timing and
traversal consequences for these public BED BayesC concepts. The BED adapter
sets the existing probability-adaptive growth rule. `NeighborWakeupControl` and
the wholesale CSR `ScheduledExecutionControl` are deliberately not used.

## 8. Callable per-chain core

`run_bed_scheduled_bayesc_chain(context)` validates once, binds borrowed hot
fields locally, allocates chain state, constructs chain RNG/scheduler state,
runs the unchanged MCMC and marker traversal against decoded storage,
accumulates summaries, and returns one typed chain result. It performs no BED
I/O, dispatch, cross-chain aggregation or R conversion.

## 9. Typed per-chain result

`BedScheduledBayesCChainExecutionResult` is the relocated authoritative former
`ChainResultSTScheduled` vocabulary: marker posterior means/PIPs, final effects
and indicators, variance/pi/VLE/VLD traces and finals, CPO diagnostics,
retained-sample count, elapsed time, failure flag and error text. It contains no
R metadata, RNG state, genotype ownership or mutable scheduler internals.

## 10. RNG ownership

One `BedScheduledBayesCChainRng` is constructed inside the callable core from
the adapter-resolved seed. Its `std::mt19937`, normal and uniform distributions
exist for one logical chain execution. No worker, static, thread-local,
cross-chain or cross-fit distribution state exists.

## 11. Genotype and I/O ownership

The adapter retains BED magic/mode validation, blocked `fseek`/`fread`, decoding
and fit-owned immutable packed storage. The core borrows it and performs only
existing in-memory marker access. No additional decode pass, disk access, file
handle, missing-value change or per-chain genotype copy was introduced.

## 12. Chain-local scheduler and sampler state

Effects, residuals, indicators, variance/pi parameters, active and candidate
lists, due buckets, visit/candidate flags, posterior accumulators, CPO, RNG and
failure/timing state remain local to a single core invocation. No mutable state
is shared among logical chains.

## 13. Native adapter

The adapter retains R decoding/validation, BED decoding, map/order preparation,
trait-chain mapping, unchanged static OpenMP task dispatch, resolved seed
mapping, context construction, task-result storage, cross-chain aggregation,
timing/failure reporting and R conversion.

## 14. Temporary adapter aliases

`MarkerMapSTScheduledChains` and `ChainResultSTScheduled` remain as binding-file
aliases to the stable typed marker-map and result. They preserve the existing
helper and aggregation spelling and are candidates for removal in Phase 11C3.
No numerical compatibility shim or fallback path remains.

## 15. Corrected frozen references

All three Phase 11B corrected raw and all three formatted references match
exactly. Phase 11A BayesC remains historical. Phase 11A BayesR/BayesRC
references remain protected.

## 16. Reproducibility

Phase 11B continues to protect exact `A; A`, `A; B; A`, normalized `1,2,2,1`,
fresh/reused process, intervening BayesR/BayesRC and different-chain-count
results. Single-route versus public multichain-one-chain remains intentionally
nonidentical because their established seed mappings differ.

## 17. Performance, memory, and I/O

`Rscript tools/benchmarks/blr_phase11c2_bed_scheduled_bayesc.R` ran five timed
repetitions after warm-up for dense 1x1/2x1/2x2, aggressive and conservative
scheduling plus 2,000-marker/200-sample workloads. Moderate medians were 0.09 s
(dense 1x1) and 0.08 s (aggressive 2x2), ranges 0.03 and 0.04 s. Phase 11B's
moderate medians were 0.10, 0.17 and 0.16 s, but these short Windows timings do
not establish an improvement. Completed-fit whole-process RSS was about
122--139 MB; it is not peak memory. The moderate BED remained 100,003 bytes,
the blocked reader is unchanged, and OS page-cache effects remain a limitation.
There is no evidence of material runtime, memory or I/O regression.

## 18. Public API and schema

Public arguments/defaults/routing, native exported signature, generated
wrappers, `NAMESPACE`, `stblr_raw_v1`, formatted fields, types, dimensions,
classes and actual `NULL` behavior are unchanged.

## 19. Protected backends

The experimental and sparse BED BayesC routes, packed-BED BayesR/BayesRC,
canonical CSR, block-eigen, multivariate sources and generated wrappers remain
unchanged from `a1c1a99`.

## 20. Tests

The combined Phase 11A/11B/11C1/11C2 focused run passed 169 expectations with
two expected opt-in skips (67, 35, 22 and 45 expectations respectively). The
enabled Phase 11B fresh-process matrix passed 37/37 with no skips. The final
full suite passed 4,792 expectations with zero failures or warnings and two
expected opt-in fresh-process skips. Native compilation/package loading and
`Rcpp::compileAttributes()` succeeded; generated wrappers remained unchanged.

## 21. Deviations and blockers

No known blocker. Conditional progress output uses binding-neutral `std::cout`
inside the existing guarded progress branch instead of `Rcpp::Rcout`; draw
sites, state transitions and returned results are unchanged. Cross-chain work
is intentionally deferred.

## 22. Recommended Phase 11C3

> centralize public scheduled packed-BED BayesC task-result aggregation and result conversion around the typed chain results, remove temporary adapter aliases, preserve public schemas and Phase 11B references, benchmark the final migrated path, and close the public multichain BayesC migration while leaving the experimental route, BayesR, and BayesRC unchanged.

## 23. Readiness marker

PHASE 11C2 COMPLETE — SCHEDULED PACKED-BED BAYESC TYPED CHAIN BOUNDARY ACTIVE
