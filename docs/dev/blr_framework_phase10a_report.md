# Unified BLR Framework Phase 10A report

## 1. Executive summary

Scheduled ordinary-CSR contracts, deterministic references, scheduler semantics,
and RNG ownership were audited without changing production execution. The only
scheduled ordinary-CSR model is BayesC. A fit-lifetime defect was reproduced:
the production `static thread_local std::normal_distribution<double>` retains a
cached variate across logically independent fits and chains. Fresh-process
references are stable, but same-process and worker-assignment results are not.

## 2. Repository baseline

Branch `master` started clean at `ab18b9c` (`Canonicalize learned-annotation CSR
BayesC`), which is the Phase 9G commit. `git diff --check` passed. R 4.4.1 UCRT,
Rtools44 GCC 13.2.0, C++17, and OpenMP were used. The baseline native compile
and complete suite passed: 4,362 expectations, zero failures, warnings, or skips.

## 3. Scheduled-route inventory

| Public R route | Model | Native entry/source | Operator | Chains/tasks | Converter/schema |
|---|---|---|---|---|---|
| `stblr_csr(method="bayesc", scheduled=TRUE)` | BayesC | `stblr_cpg_omp_csr_scheduled`; `src/st_cpg_omp_csr_scheduled.cpp` | streamed upper-triangle CSR | trait-major/chain-minor tasks; OpenMP static | native named raw list; `stblr_raw_v1` formatter |
| `stblr_csr(method="bayesr", scheduled=TRUE)` | BayesR | rejected in R | none | none | none |

No scheduled ordinary-CSR BayesRC/SBayesRC or conditionally scheduled ordinary
route exists. Scheduled packed-BED BayesR/BayesRC and scheduled individual-level
samplers are distinct implementations and were inspection-only protected scope.

## 4. Scheduler-state inventory

| State/control | Type/dimension/init | Owner/reset | Update timing/effect |
|---|---|---|---|
| `full_sweep_every` | integer | immutable call control | full when skipping is disabled, value <=0, or `it % value == 0` |
| `null_skip_base/max` | integers | immutable call control | adaptive next-visit interval, capped by max |
| candidate threshold/lifetime | double/integer | immutable controls | entry after update; expiry when `it-last_interesting > lifetime` |
| burn-in-only flag | bool | immutable | disables skipping after burn-in |
| neighbor controls | bool/double/integer | immutable; CSR borrowed | wake after effect change; CSR order; max zero means unlimited |
| `scheduled` | vector of marker buckets | chain-local, new per task | due markers indexed by future iteration modulo bucket count |
| `scheduled_at` | marker integer, -1 | chain-local | prevents stale bucket visits |
| `last_updated` | marker integer, -1 | chain-local | prevents duplicate update in an iteration |
| `is_candidate` | marker byte, 0 | chain-local | candidate status |
| candidate/active lists and membership bytes | marker-capacity vectors/bytes | chain-local | ordered traversal; compacted every 50 iterations |
| `last_interesting` | marker integer, -1e9 | chain-local | lifetime/expiry clock |

Production has no returned attempted/skipped/full-sweep/candidate/wakeup
diagnostics. Those names are internal result vocabulary, not new public fields.

## 5. Scheduling semantics

Iteration zero is a full sweep. `full_sweep_every=1` makes every iteration full;
a frequency above the run length leaves only iteration zero full. Burn-in and
thinning do not change the modulo convention. When skipping is burn-in-only,
all post-burn-in iterations are full.

A null update selects an interval from probability bands (4x, 2x, base, half
base, or 1), caps it, and adds task-local uniform integer jitter. Non-null or
high-probability markers enter/refresh candidates. Candidate lifetime zero still
allows the current iteration. Non-full traversal is active list, candidate list,
then the due bucket; duplicate visits are suppressed. Skipped markers consume no
marker-update RNG.

Neighbor wake-up traverses borrowed CSR adjacency order after an effect change
exceeding the threshold. Membership flags suppress duplicate list insertion,
while `last_interesting` is refreshed. Wake-up may affect a later list in the
same iteration or the next eligible traversal.

## 6. Marker traversal

Full sweeps use the fixed descending rank of `(wy/ww)^2`. Scheduled sweeps use
stable active/candidate/bucket insertion order. Traversal is chain-local, but
static OpenMP assigns trait-chain tasks to worker threads. Skips omit all marker
draws, so candidate and wake-up decisions alter the subsequent RNG stream.
Diagnostics in production do not expose attempted versus completed updates.

## 7. RNG ownership inventory

| Distribution | Storage/owner | Cache/reset | Risk |
|---|---|---|---|
| `std::mt19937` | task-local, seeded from trait-chain seed | reconstructed per task | safe engine ownership |
| uniform integer jitter | task-local | reconstructed per task | safe |
| uniform real | `static thread_local` in marker helper | persistent worker-thread object; stateless in practice | safe only under library assumptions |
| normal | `static thread_local` in marker helper | cached second variate; never reset per task/fit | stateful cross-fit and thread-assignment risk; Phase 10B correction required |
| chi-square | local variance update | reconstructed at draw site | fit/chain local |
| gamma | local probability update | reconstructed at draw site | fit/chain local |

Explicit chain seeds seed engines, not distribution-internal state. Reseeding or
reconstructing an engine therefore does not reset the persistent normal cache.

## 8. `static thread_local` investigation

Both persistent distributions are in `src/st_cpg_omp_csr_scheduled.cpp` inside
the marker update helper. The validation-only diagnostic reproduces the exact
normal pattern. On one and two OpenMP workers, a first draw cached its paired
normal; reconstructing the engine with the same seed and drawing through the
same distribution returned the cached paired variate. Calling
`distribution.reset()` restored the fresh seeded result exactly.

The distribution lifetime is worker-thread lifetime, not fit, chain, trait, or
task lifetime. Risk classification: normal = `requires Phase 10B correction`;
uniform real = `safe only under current standard-library/stateless assumptions`.

## 9. Fresh-process reproducibility

Three configurations were generated in separate R processes. Repeated
process-isolated generation matched their complete normalized raw and formatted
objects exactly. One-core, two-core, and explicit-chain-seed fresh baselines are
therefore retained. Timing and the recorded LD prefix are the only normalized
execution metadata.

## 10. Same-process reproducibility

For the compact two-chain fixture, `A; A` was exact. `A; B; A`, where B changed
seed and scheduling/draw consumption, was not exact. The first formatted
difference was `bm[1,1]`: expected `-0.008899249`, observed `0.007956757`.
Scheduled/unscheduled ordering is asymmetric: ordinary unscheduled fits do not
use the scheduled persistent object, whereas an intervening scheduled fit can
change its cached state. Different chain counts likewise change consumption.

## 11. Thread-assignment reproducibility

The same explicit chain seeds under the sequence 1,2,2,1 cores produced:
1-vs-2 unequal, repeated 2-core unequal, and returning to 1-core also unequal
(including sampled values and, in one run, last-bit aggregation differences).
This localizes sensitivity to persistent worker-thread
distribution state/task assignment rather than output ordering. Tests cover the
one- and two-worker diagnostic; four-worker experimentation adds no new route
because the compact fixture has two tasks.

## 12. Frozen references

Scheduled CSR BayesC has 3 raw and 3 formatted references:

1. `dense_one`: one chain/core, every iteration full, no skipping/wakeup;
2. `skip_two_one`: explicit two-chain seeds, one core, adaptive skipping,
   candidates, finite lifetime, and neighbor wake-up;
3. `skip_two_two`: explicit seeds, two cores, burn-in-only skipping, candidates,
   bounds, and neighbor wake-up.

Scheduled CSR BayesR has 0 references because the public route is rejected.
Same-process unstable sequences are diagnostic assertions, not canonical data.
The generator is explicitly marked as a manual maintenance tool.

## 13. Reduction to unscheduled execution

Dense controls (`full_sweep_every=1`, skip base/max 1, threshold/lifetime 0,
wakeup disabled) did not reduce exactly to canonical unscheduled BayesC. The
first stable sampled-output difference occurs in marker summaries. This is
classified as distinct implementation/RNG ownership and update-policy behavior,
not a failed canonical reduction; no production formula was changed.

## 14. Typed scheduled contracts

`src/blr_scheduled_execution_types.h` defines binding-neutral sweep, null-skip,
candidate, neighbor-wakeup, and execution controls. Validation covers positive
dimensions/counts, seed length, skip ordering, thresholds/lifetimes, immutable
borrowed friend ownership, friend dimensions, and state dimensions. Borrowed
resources must outlive future execution; no CSR or friend-list per-chain copy is
represented.

## 15. Typed state and result vocabulary

Configuration is separated from chain-owned mutable marker state. The state
vocabulary mirrors production's scheduled-at, last-updated, candidate/list
membership, and last-interesting fields. Diagnostics/results name attempted,
completed, skipped, full-sweep, candidate, expiry, and wake-up categories while
accurately marking current production scheduler diagnostics and optional chain
payloads unavailable.

## 16. Rcpp validation bridge

`blr_phase10a_validate_scheduled_execution_cpp()` is internal. It constructs and
validates contracts, returns normalized metadata plus `validated=TRUE`,
`invokes_sampler=FALSE`, and `consumes_rng=FALSE`. It has no public export or
NAMESPACE entry. A second internal bridge performs only the standard-library
distribution diagnostic.

## 17. Production behavior statement

`src/st_cpg_omp_csr_scheduled.cpp` remains byte-identical
(`fdac03befb742f4f6fa7c22ccbbbc920`). Routing, scheduler state, formulas, RNG,
conversion, signatures, `stblr_raw_v1`, and formatted schema are unchanged.

## 18. Performance and memory baseline

`Rscript tools/benchmarks/blr_phase10a_scheduled_csr.R` uses a warm-up and five
repetitions. For 2,000 markers/30 iterations, observed elapsed distributions
were: dense 1-chain/1-core .04,.06,.05,.03,.03 s (median .04); dense
2-chain/1-core .05,.06,.17,.25,.21 (median .17); dense 2-chain/2-core
.09,.13,.11,.12,.13 (median .12); aggressive skipping .03,.02,.04,.03,.03 s;
conservative skipping five times .11 s; wake-up configuration
.03,.03,.04,.03,.05 (median .03). Empty benchmark CSR means wake-up controls
are measured without neighbor edges. Production does not return update counts,
so those are recorded unavailable rather than inferred.

Memory is process RSS sampled through `ps` before and after completed fits; it
ranged from 127.8 to 136.3 MiB after the primary completed-fit repetitions. It
is not peak RSS and interval sampling was unavailable. This is a pre-migration
per-backend baseline, not a BayesC/BayesR comparison.

## 19. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior, group, and learned-annotation
sources; block-eigen; scheduled BED/individual; BED BayesRC; and multivariate
sources were hash/Git audited and unchanged. NAMESPACE and public signatures are
unchanged. Generated wrappers changed only for the two internal bridges.

## 20. Tests

Baseline full suite: 4,362 pass, zero failure/warning/skip. The opt-in
fresh-process Phase 10A run passed 34/34 with no failure, warning, or skip. The
ordinary focused run contributed 43 passing expectations; process spawning is
only activated by its environment flag. The final full suite passed 4,405
expectations with zero failures, warnings, or skips. Phase 10A focused
tests cover contract round trips and errors, source ownership, cached normal
state, three process-isolated fixture sets, same-process ordering, core order,
dense reduction, and protected hashes. Process-isolated execution is opt-in in
ordinary sandboxed tests because nested process pipe creation is denied; the
manual generator was run successfully outside that restriction.

## 21. Confirmed defects or risks

| Source | Reproduction | Statistical consequence | Migration consequence | Recommended correction |
|---|---|---|---|---|
| persistent normal in scheduled marker helper | engine reseed diagnostic; `A;B;A`; 1,2,2,1 cores | sampled effects and downstream summaries depend on prior fit/worker state | cannot freeze same-process canonical trajectory or migrate safely | make distribution state explicit and chain-owned, reset per task/fit |
| persistent uniform real | source ownership audit | no cached variate observed/expected, but ownership remains worker rather than chain | ownership contract is misleading | make explicit chain-owned state with normal correction |
| no scheduler counters returned | source/result audit | no statistical effect; update savings cannot be observed directly | benchmark diagnostics limited | retain public schema; consider internal counters only in later bounded work |

## 22. Recommended Phase 10B

Replace scheduled ordinary-CSR stateful static/thread-local distributions with
explicit chain-owned distribution state, preserving current within-fresh-process
numerical trajectories where possible, and establish deterministic fit-local
reproducibility before migrating the scheduled numerical core.

## 23. Readiness marker

PHASE 10A COMPLETE — SCHEDULED CSR CONTRACTS AND RNG AUDIT READY
