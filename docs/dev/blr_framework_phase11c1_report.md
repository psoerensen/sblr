# Unified BLR Framework Phase 11C1 report

## 1. Executive summary

The corrected public scheduled packed-BED BayesC logical-chain execution body
was mechanically extracted without changing Phase 11B behavior.

## 2. Repository baseline

The clean baseline was `master` at Phase 11B commit `2ce809c`. Initial status
and `git diff --check` were clean. R 4.4.1 with Rtools44/GCC and OpenMP compiled
the baseline, whose full suite passed 4,724 expectations.

## 3. Production call graph

`stblr_bed(method="bayesc")` prepares aligned inputs and calls
`stblr_cpg_omp_bed_marker_scheduled_chains()`. That exported native adapter in
`st_cpg_omp_individual_scheduled_chains.cpp` validates and decodes BED data,
prepares tasks, and calls `run_one_scheduled_bed_chain()` once per trait-chain.
The helper owns the MCMC/marker execution. Native aggregation follows, then the
existing inline `stblr_raw_v1` construction and R formatter.

## 4. Selected extraction target

The public multichain route was selected because it contains its own task
orchestration and calls its own complete logical-chain numerical helper. The
417-line `run_one_scheduled_bed_chain()` definition is the true active
public-path MCMC body.

## 5. Unselected route

`stblr_cpg_omp_bed_marker_scheduled()` is an experimental single-chain route
with a separate numerical implementation. The public route does not delegate
to it, so it remains byte-identical to the starting commit.

## 6. Extraction seam

The moved block begins at the stable declaration
`static ChainResultSTScheduled run_one_scheduled_bed_chain(` and ends at the
closing brace immediately before `// [[Rcpp::export]]` for the public native
entry. The source now includes the guarded implementation header at that seam.

## 7. Files changed

The selected multichain source now includes the new implementation header.
Phase-era source hashes were updated to permanent post-extraction hashes. A
Phase 11C1 structural/reference test and this report were added; the plan and
capability matrix record migration-in-progress status.

## 8. Number of lines moved

Exactly **417 lines** moved. A line-for-line comparison against commit
`2ce809c` reports 417 old and 417 new body lines with no differences. Only the
include guard, implementation-detail comment and replacement include are new.

## 9. Extracted content

The header contains logical-chain seed/RNG construction, effects, residuals,
variance and global-pi state, adaptive scheduler vectors/buckets, MCMC loop,
full/active/candidate/due traversal, marker draws, genotype access, residual
updates, accumulation, final chain summaries, timing and failure capture.
Task allocation, OpenMP task dispatch and native multichain aggregation remain
in the adapter around the extracted per-chain body.

## 10. RNG ownership preservation

One `BedScheduledBayesCChainRng` remains per logical trait-chain, owning its
`std::mt19937`, normal and uniform distributions for one chain execution. No
static, thread-local, worker-owned or fit-persistent state was introduced.

## 11. Scheduler preservation

Full-sweep timing/frequency, null-skip growth/reset, candidate threshold and
lifetime, due buckets, burn-in-only policy, marker order and
active/candidate/due traversal are byte-preserved. Skipped markers still avoid
marker RNG and genotype update work. OpenMP task scheduling remains static.

## 12. Genotype and I/O preservation

BED validation, blocked `fseek`/`fread`, SNP-major decoding, missing handling,
centering/scaling and fit-local immutable packed storage remain outside the
header and unchanged. The header performs only existing in-memory marker
access; it adds no disk read or per-chain genotype copy.

## 13. Existing result conversion

Native aggregation and the complete named raw-list construction remain in the
selected source after the include. Field names, dimensions, classes, actual
`NULL`, schema/version and formatted output are unchanged.

## 14. Corrected frozen references

All three Phase 11B raw and all three formatted references match exactly after
normalizing only declared timing/core metadata. Phase 11A BayesC remains a
historical pre-correction artifact and was not regenerated.

## 15. Reproducibility

Phase 11B continues to protect exact `A; A`, `A; B; A`, normalized `1,2,2,1`,
fresh/reused process, intervening BayesR/BayesRC and different-chain-count
results. Single versus multichain-one-chain remains intentionally nonidentical
because the established multichain seed includes its chain offset.

## 16. Protected backends

The experimental single-chain and sparse BayesC routes, BayesR, BayesRC,
canonical CSR, block-eigen, multivariate sources, generated wrappers and
`NAMESPACE` remain unchanged from `2ce809c`.

## 17. Tests

The focused Phase 11A/11B/11C1 run passed 124 expectations with two expected
opt-in skips (67 + 35 + 22 expectations respectively). The Phase 11B
fresh-process matrix passed 37/37 expectations with no skips. The final full
suite passed 4,747 expectations with zero failures or warnings and two
expected opt-in fresh-process skips; the separately enabled Phase 11B
fresh-process check passed. Native compilation and package loading succeeded.

## 18. Deviations and blockers

No blocker. The extracted unit is the complete per-logical-chain numerical
body; public task dispatch and native aggregation deliberately remain in the
adapter until typed-boundary work.

## 19. Recommended Phase 11C2

> replace the lexically dependent corrected packed-BED BayesC execution include with an explicit typed BED execution context and callable core, preserving logical-chain RNG ownership, scheduler semantics, genotype storage ownership, and existing result conversion until corrected references pass again.

## 20. Readiness marker

PHASE 11C1 COMPLETE — CORRECTED SCHEDULED PACKED-BED BAYESC EXECUTION BLOCK EXTRACTED
