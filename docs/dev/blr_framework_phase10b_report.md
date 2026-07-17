# Unified BLR Framework Phase 10B report

## 1. Executive summary

Scheduled ordinary-CSR BayesC RNG ownership was corrected without changing its
scheduler or statistical formulas. Each logical trait-chain task now owns its
engine, normal distribution, and uniform distribution for one chain execution.
Fit-local, call-order, and worker-assignment reproducibility are exact.

## 2. Repository baseline

Branch `master` started clean at Phase 10A commit `fb9fd03` (`Audit scheduled
CSR execution and RNG ownership`). R 4.4.1 UCRT, Rtools44 GCC 13.2.0, C++17,
and OpenMP were used. Baseline compilation and the 4,405-expectation full suite
passed. Phase 10A reproduced the defect before editing.

## 3. Confirmed Phase 10A defect

The marker helper in `src/st_cpg_omp_csr_scheduled.cpp` held worker-owned
`static thread_local` normal and uniform distributions. The normal distribution
cached a second variate across engine reconstruction. `A;B;A` first differed at
formatted `bm[1,1]` (`-0.008899249` versus `0.007956757`), and `1,2,2,1` cores
was not reproducible. Consequently posterior marker and variance summaries could
depend on earlier fits and OpenMP worker assignment despite identical seeds.

## 4. Files changed

- `src/st_cpg_omp_csr_scheduled.cpp`: routes marker uniform/normal draws through
  a task-local `ScheduledChainRng`; no scheduler statement or formula changed.
- `src/blr_scheduled_execution_types.h`: defines the chain RNG and explicit
  ownership/lifetime contract.
- `src/blr_phase1_rcpp.cpp` and generated `R/RcppExports.R`/
  `src/RcppExports.cpp`: ownership round trip and validation-only RNG diagnostic.
- Phase 10A/earlier structural tests: replace intentionally obsolete production
  and generated-wrapper hashes while retaining Phase 10A audit artifacts.
- Phase 10B tests, helper, fixtures, generator, benchmark, report, plan, and
  capability matrix: permanent corrected-path coverage and documentation.

## 5. Corrected ownership model

After the unchanged trait/chain seed formula resolves `task_seed`, the task
constructs one `ScheduledChainRng`. It owns `std::mt19937`,
`std::normal_distribution<double>`, and
`std::uniform_real_distribution<double>`. Its lifetime is the complete logical
chain execution, it is never shared, and destruction occurs before the task
ends. OpenMP thread number is not an ownership or seed identity. Chi-square and
gamma distributions have variable parameters and remain newly constructed at
their unchanged update sites inside the owning chain execution.

## 6. Numerical changes and nonchanges

Intentionally changed: inherited cross-fit, cross-chain, and worker-thread
cached distribution state. Of the Phase 10A fresh fixtures, `dense_one` and
`skip_two_two` remain byte-exact; sequential `skip_two_one` changes as expected
because chain 2 no longer inherits chain 1's cached normal (`bm[1,1]`
`-0.008899249` to `-0.007343946`).

Unchanged: seed formulas, engine algorithm, distribution algorithms/parameters,
logical draw sites, skipped-marker no-draw rule, posterior/inclusion/effect,
residual, variance, global parameter formulas, traversal, task order, OpenMP
static scheduling, accumulation, and conversion.

## 7. Source audit

No `static thread_local`, `thread_local`, fit-persistent normal, worker-owned
distribution, old/new selector, or fallback remains in scheduled ordinary-CSR
production. The sole stateful marker distributions are members of
`ScheduledChainRng`, constructed once per trait-chain task. The Phase 10A
unsafe-pattern diagnostic remains isolated in validation-only bridge code.

## 8. Post-correction references

Three distinct Phase 10B configurations provide 3 raw and 3 formatted exact
references under `blr_phase10b_scheduled_csr`: dense one-chain/one-core;
two-chain/one-core adaptive skipping with explicit seeds, candidates, lifetime,
and wake-up; and two-chain/two-core burn-in-only skipping and wake-up. Phase 10A
references were neither regenerated nor overwritten. `keep_chains=TRUE` remains
unsupported for scheduled CSR and is protected as rejection behavior.

## 9. Same-process reproducibility

`A;A`, `A;B;A`, scheduled A / unscheduled fit / scheduled A, and one-chain A /
multi-chain fit / one-chain A are exact. An intervening chain consuming an odd
number of normals does not affect a reconstructed same-seed chain RNG.

## 10. Core/thread reproducibility

The normalized `1,2,2,1` core sequence is exact. Explicit two-chain seeds yield
identical sampled and aggregated results when chains run sequentially on one
worker or concurrently on two workers. Task order and OpenMP `schedule(static)`
remain unchanged.

## 11. Fresh-process reproducibility

All three fresh-process results compare exactly with reused-process results and
their post-correction fixtures. The opt-in process-isolated suite passed 39
expectations without failures, warnings, or skips.

## 12. Scheduler preservation

Iteration-zero full sweep, full-sweep modulo, null skip base/max and adaptive
growth, candidate entry/lifetime/expiry, burn-in-only skipping, neighbor order,
threshold/max, active/candidate/due traversal, and skipped-marker no-draw
behavior are unchanged. Production still exposes no scheduler event counters;
the unchanged schema therefore reports benchmark attempted/skipped/full-sweep,
candidate, and wake-up counts as unavailable rather than inventing them.

## 13. Dense reduction

Dense scheduled controls still do not reduce exactly to canonical unscheduled
BayesC. The first difference remains in sampled marker summaries. This is
unrelated to ownership correction and reflects distinct scheduled production
implementation/update behavior. No formulas were changed to force reduction.

## 14. Public API and schema

Public arguments, defaults, routing, seed mapping, native public signature,
`NAMESPACE`, `stblr_raw_v1`, formatted fields, chain order, actual NULL behavior,
and `keep_chains` rejection are unchanged. Scheduled BayesR remains unsupported.

## 15. Protected backends

Canonical BayesC, BayesR, SBayesRC, fixed-prior, group, learned-annotation,
block-eigen, BED scheduled, individual scheduled, and multivariate sources are
byte-identical to the Phase 10B starting commit. Generated wrappers changed only
for the internal diagnostic bridge.

## 16. Performance and memory

`Rscript tools/benchmarks/blr_phase10b_scheduled_csr.R` reused Phase 10A's
warm-up and five-repetition workloads. For 2,000 markers/30 iterations, Phase
10B medians were dense 1-chain/1-core .08 s, dense 2-chain/1-core .18 s, dense
2-chain/2-core .11 s, aggressive skipping .03 s, conservative skipping .10 s,
and wake-up .04 s. Phase 10A medians were .04, .17, .12, .03, .11, and .03 s;
these short Windows timings are noisy and show no material systematic regression.

Completed-fit process RSS was 127.8–137.6 MiB versus Phase 10A's 127.8–136.3
MiB. This is completed-fit RSS sampled using `ps`, not peak memory. The small
per-chain distribution objects cause no meaningful memory regression. Scheduler
event counts remain unavailable in the public result.

## 17. Tests

Phase 10A audit behavior and artifacts remain covered. Phase 10B covers typed
ownership validation, reconstruction, identical seeds, odd intervening draws,
source restrictions, 3/3 raw and 3/3 formatted references, `A;A`, `A;B;A`,
intervening scheduled/unscheduled work, different chain counts, normalized
1/2-core execution, fresh/reused processes, dense reduction, public limitations,
and protected hashes. Phase 10A audit tests passed 40/40; Phase 10B tests passed
39/39 in both ordinary and process-isolated modes. The final full suite passed
4,441 expectations with zero failures, warnings, or skips.

## 18. Deviations and blockers

Scheduled production does not return scheduler event counters or retained chain
payloads; Phase 10B preserves that schema and documents those benchmark values
as unavailable. The Phase 10A one-core sequential multichain numerical baseline
changes intentionally because it embodied the confirmed defect. No blocker
remains for deterministic scheduled execution.

## 19. Recommended Phase 10C

Establish typed scheduled execution inputs/results around the corrected
deterministic production implementation and mechanically extract the scheduled
numerical core while preserving scheduler semantics and post-correction
references.

## 20. Readiness marker

PHASE 10B COMPLETE — SCHEDULED CSR RNG OWNERSHIP CORRECTED
