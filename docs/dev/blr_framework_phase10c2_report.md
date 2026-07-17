# BLR Framework Phase 10C2 Report

## 1. Executive summary

The corrected scheduled ordinary-CSR BayesC execution now runs behind a typed, binding-neutral `CsrScheduledBayesCExecutionContext` and returns a typed `CsrScheduledBayesCExecutionResult` from `run_csr_scheduled_bayesc()`. The Phase 10B chain-owned RNG model, scheduler semantics, numerical trajectory, inline R converter, public route, and schema are preserved.

## 2. Repository baseline

- Branch: `master`
- Starting commit and Phase 10C1 commit: `b6c436f Extract corrected scheduled CSR BayesC execution block`
- Initial status: clean
- Initial `git diff --check`: clean
- Toolchain: R 4.4.1 on Windows with the Rtools44 GNU C++17/OpenMP toolchain
- Baseline Phase 10A–10C1 focused suite: 124 expectations passed
- Baseline Phase 10C1 suite: 44 passed
- Phase 10C1 full suite: 4,486 passed, 0 failures, 0 warnings, 0 skips
- Corrected raw and formatted references, same-process sequences, normalized core ordering, and isolated-process checks were exact

## 3. Lexical dependency inventory

| Symbols | Type / dimensions | Classification | Ownership and lifetime | Context/result treatment |
|---|---|---|---|---|
| `ld` | `STLDCSR`; marker-indexed CSR values, row pointers, indices | immutable borrowed input | adapter-owned, outlives call | `context.ld`; no result field |
| `wy_mat`, `ww_mat`, `b_mat` | `arma::mat`, traits × markers | immutable borrowed input | adapter-owned prepared statistics/initial effects | `wy`, `ww`, `initial_b` |
| `yy_vec` | `arma::vec`, traits | immutable borrowed input | adapter-owned | `context.yy` |
| `n` | `std::vector<int>`, traits | immutable borrowed input | R-decoded adapter storage | `sample_sizes`; binding converter also uses it |
| `ssb_prior_mat`, `sse_prior_mat` | `arma::mat`, traits × traits | immutable borrowed input | adapter-owned | context prior references |
| `B`, `E` | `arma::mat`, traits × traits | immutable initial state | adapter-owned; copied only into chain scalar state | `initial_B`, `initial_E` |
| `pi` | `std::vector<double>`, length 2 | immutable initial state | adapter-owned; copied per chain | `initial_pi` |
| `order` | `std::vector<int>`, markers | immutable borrowed input | adapter-owned and shared read-only | `marker_order`; no per-chain copy |
| `d_init`, `r_init` | vector-of-vectors, traits × markers when active | immutable borrowed input | R-decoded adapter storage | `initial_d`, `initial_r` |
| `m`, `nt` | scalar marker/trait counts | immutable scalar control | copied metadata | `marker_count`, `trait_count` |
| `nit`, `nburn`, `nthin` | scalar iteration controls | immutable scalar control | Phase 10A control | `scheduled.iterations/burnin/thinning` |
| `nchains`, `ncores` | scalar execution controls | immutable scalar control | Phase 10A control | `scheduled.chains/cores` |
| `seed`, `chain_seeds` | base seed and optional chain vector | immutable borrowed/control | adapter vector outlives call | Phase 10A execution control |
| `full_sweep_every` | scalar | scheduler control | copied into typed sweep control | `scheduled.sweep` |
| `null_skip_base/max`, burn-in-only | scalars | scheduler control | typed skip control | `scheduled.skip` |
| candidate threshold/lifetime | scalar | scheduler control | typed candidate control | `scheduled.candidate` |
| wake-up enabled/threshold/max and CSR friends | scalars plus borrowed neighbor storage | scheduler control / borrowed input | CSR storage remains shared read-only | `scheduled.neighbor` plus `context.ld` |
| `nub`, `nue`, `adjE`, pi prior parameters | scalar model/prior controls | immutable scalar control | copied metadata | direct context fields |
| update/init/rebuild flags | booleans | immutable scalar control | copied metadata | direct context fields |
| task descriptor, output, failure, timing arrays | task-count vectors/matrices | mutable task-local/core-owned | allocated once inside core | timing returned; failures translated before return |
| effects, residuals, states, variances, pi | marker vectors/scalars per task | mutable chain-local state | allocated inside logical chain | aggregated result only |
| active/candidate/due buckets, skip and wake-up state | marker/bucket vectors per task | mutable chain-local scheduler state | allocated within chain execution | not publicly returned |
| `ScheduledChainRng` and jitter distribution | engine/distributions per task | chain-owned RNG state | one logical chain execution | not exposed in result |
| posterior task accumulators | task × marker/trace matrices | diagnostic/output accumulators | core-owned | aggregated into typed result |
| `bm_mat` through `nsamples_vec` | aggregate Armadillo outputs | result outputs | core-owned then moved | typed result fields |
| `task_seconds` | task-count vector | diagnostic output | core-owned then moved | typed result |
| R names, dimnames, schema/class strings, `n_trace` | binding-only metadata | excluded | binding source only | not present in context; converter retains access |

No former execution dependency remains lexically resolved from the wrapper.

## 4. Files changed

- `src/blr_csr_scheduled_bayesc_types.h`: added the backend-specific typed context, typed result, and field-specific validator.
- `src/blr_csr_scheduled_bayesc_core_impl.h`: converted the lexical block into the callable templated core, bound original names to context fields, and returned moved native aggregates.
- `src/st_cpg_omp_csr_scheduled.cpp`: constructs Phase 10A controls/context, calls the core, and exposes typed-result aliases to the existing converter.
- `tests/testthat/test-blr-framework-phase10c2.R`: added permanent architecture, ownership, exact-reference, reproducibility, reduction, and protection tests.
- `tools/benchmarks/blr_phase10c2_scheduled_csr.R`: added an identical-workload typed-path benchmark wrapper.
- Framework implementation plan and capability matrix: marked typed execution migration active.
- This report records inventory, ownership, validation, behavior, performance, and test evidence.

## 5. Typed execution context

The templated context borrows the CSR operator; prepared summary statistics; initial effects, residual/state inputs, covariance values, priors, marker order, sample sizes, and seed vector. It copies only small model flags/scalars and references one `ScheduledExecutionControl`. All borrowed resources must outlive the call. It contains no Rcpp, SEXP, schema, names, dimnames, class, or Python types.

## 6. Scheduler contracts

The adapter activates `ScheduledSweepControl`, `NullSkipControl`, `CandidateControl`, `NeighborWakeupControl`, RNG ownership metadata, and the encompassing `ScheduledExecutionControl`. Iteration zero remains a full sweep; non-full traversal remains active markers, candidates, then due buckets; skipped markers still avoid marker-update RNG.

## 7. Callable numerical core

The binding-neutral signature is:

```cpp
template <class Operator>
CsrScheduledBayesCExecutionResult run_csr_scheduled_bayesc(
    const CsrScheduledBayesCExecutionContext<Operator>& context);
```

It validates the prepared boundary, allocates task outputs, creates chain RNG/scheduler state, executes the unchanged OpenMP static task loop and MCMC scheduler, detects task failures, aggregates current native outputs, and returns the typed result. It performs no R decoding or R object construction.

## 8. Typed execution result

The result owns aggregated marker means/states and chain-summary matrices; final effects/residuals; variance, pi, VLE, and VLD traces; final variance/pi values; retained-sample counts; task timings; and marker/trait/chain/task dimensions. No new public scheduler counters or unsupported chain payloads were invented.

## 9. RNG ownership

Each logical trait-chain still constructs exactly one `ScheduledChainRng` after final seed resolution. Its `std::mt19937`, normal distribution, and uniform distribution are chain-owned and fit-bounded. Variable-parameter gamma and chi-square distributions remain at their existing draw sites. There is no static, thread-local, worker-owned, shared, or persistent distribution state.

## 10. Scheduler-state ownership

Active/candidate lists, due buckets, scheduled/last-updated indices, candidate flags/lifetimes, wake-up state, effects, residuals, variances, posterior accumulators, and RNG are chain-local. CSR rows/indices/values, marker order, statistics, priors, seeds, and neighbor data are borrowed immutable. Only preallocated task output slots are shared, with each task writing its own row.

## 11. Native adapter

The binding source retains R argument decoding, public validation, Armadillo conversion, shared-`ww` checks, CSR construction, marker ranking, Phase 10A control construction, context construction, one core call, temporary typed-result aliases, inline conversion, and exception translation.

## 12. Existing converter

The Rcpp converter remains inline and behaviorally unchanged. Temporary aliases map result fields to its established names: `bm_mat`, `dm_mat`, chain-summary matrices, `b_out_mat`, `r_out_mat`, `d_out_mat`, trace matrices, final variance/pi vectors, `nsamples_vec`, and `task_seconds`. Phase 10C3 should remove these aliases while centralizing conversion.

## 13. Corrected frozen references

Phase 10B corrected raw references: 3/3 exact. Corrected formatted references: 3/3 exact. Phase 10A defective references remain unchanged historical artifacts.

## 14. Reproducibility

Repeated `A; A`, intervening scheduled `A; B; A`, normalized `1,2,2,1`, intervening unscheduled work, different chain counts, and explicit chain seeds remain exact. Isolated fresh-process comparison is run separately against reused-process results.

## 15. Dense reduction

Dense scheduled execution remains non-identical to canonical unscheduled BayesC at the previously documented sampled-marker difference. Typed-boundary activation did not alter that classification, and no formulas were changed to force reduction.

## 16. Performance and memory

The Phase 10C2 benchmark reused the exact Phase 10B/10A driver: 2,000 markers, 30 iterations, dense 1-chain/1-core, dense 2-chain/1-core, dense 2-chain/2-core, aggressive and conservative skipping, and neighbor wake-up, with warm-up and five timed repetitions.

Phase 10C2 medians were 0.03, 0.05, 0.05, 0.01, 0.04, and 0.01 seconds respectively. Phase 10B medians were 0.08, 0.18, 0.11, 0.03, 0.10, and 0.04 seconds. These short Windows timings are noisy and do not support a speed claim; they show no regression. Completed-fit RSS after fits was approximately 128.7–142.6 MiB. The method is process RSS after completed fits, not peak RSS; no meaningful typed-boundary memory regression is evident.

## 17. Public API and schema

Public arguments, defaults, route, native exported signature, generated wrappers, `NAMESPACE`, `stblr_raw_v1`, field names/order/types/dimensions, actual `NULL` fields, and formatted fit structure are unchanged. Scheduled BayesR remains unsupported.

## 18. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior BayesC, group BayesC, learned-annotation BayesC, block-eigen, BED/individual scheduled backends, and multivariate sources remain byte-identical to the starting checkpoint.

## 19. Tests

- Phase 10C2 suite: 61 passed.
- Phase 10A–10C2 focused suite: 185 passed.
- Isolated fresh-process Phase 10B matrix: 39 passed.
- Full suite: 4,547 passed.
- Failures: 0; warnings: 0; skips: 0.

## 20. Deviations and blockers

Binding-neutral diagnostic output uses `std::cout` instead of `Rcpp::Rcout`; message content, placement outside the marker hot loop, control flow, RNG, and returned state are unchanged. The initial compile exposed that a namespace-level implementation header cannot be included at a function-local lexical seam; the include was moved to the existing helper/core boundary before the exported wrapper, with the call retained at the approved execution seam. No numerical blocker remains.

## 21. Recommended Phase 10C3

> centralize the typed scheduled result-to-R conversion into one named binding-layer converter, remove temporary result aliases, retain wrapper/output aggregation semantics, replace migration-era structural checks, benchmark the final migrated path, and close the scheduled BayesC migration.

## 22. Readiness marker

PHASE 10C2 COMPLETE — SCHEDULED CSR TYPED EXECUTION BOUNDARY ACTIVE
