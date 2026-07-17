# Unified BLR Framework Phase 11A report

## 1. Executive summary

Individual-level and packed-BED backend families were audited without changing
production execution. All discovered “individual-level” production samplers
operate on packed SNP-major BED data. BayesC has a confirmed worker-persistent
distribution risk; BayesR and BayesRC use chain-local distributions. Shared
scheduled infrastructure is only partially appropriate.

## 2. Repository baseline

The clean baseline was `master` at Phase 10D commit `2a3cddb`. Initial status
and `git diff --check` were clean. R 4.4.1 with Rtools44/GCC and OpenMP compiled
the package; the baseline full suite passed 4,623 expectations.

## 3. Route inventory

| Family | Public route | Model | Native entry/source | Mode/status |
|---|---|---|---|---|
| packed BED | `stblr_bed(method="bayesc")` | binary BayesC | `stblr_cpg_omp_bed_marker_scheduled_chains`, `st_cpg_omp_individual_scheduled_chains.cpp` | scheduled multichain active |
| packed BED internal | none | binary BayesC | `..._sparse`, `st_cpg_omp_individual.cpp` | sparse unscheduled/experimental |
| packed BED internal | none | binary BayesC | `..._scheduled`, `st_cpg_omp_individual_scheduled.cpp` | scheduled single-chain/experimental |
| packed BED | `stblr_bed(method="bayesr")` | K-component BayesR | `..._scheduled_chains_bayesr`, corresponding C++ source | adaptive scheduled multichain |
| packed BED | `stblr_bed(method="bayesrc")` | annotation-component BayesRC | `..._scheduled_chains_bayesrc`, corresponding C++ source | multichain but unscheduled full sweep |
| legacy public helper | `stblr_bed_marker()` | BayesC | scheduled single/multichain native entries | packed-BED marker route |

No in-memory dense genotype public route was found in this family. Parallel
units are trait-chain tasks; BED decoding may separately use OpenMP static
marker blocks.

## 4. Statistical-model inventory

BayesC uses a binary null/non-null Gaussian marker prior and global pi/vb.
BayesR uses K components, Dirichlet probabilities and component-scaled marker
variance. BED BayesRC uses annotation-dependent probit-stick component
probabilities and coefficient updates, not canonical CSR SBayesRC's complete
execution contract. BayesRC always performs full sweeps; its “scheduled chains”
name describes task execution, not adaptive marker scheduling.

## 5. Genotype representations

All routes require PLINK BED magic bytes and SNP-major mode. Selected marker
rows are read by `fseek`/`fread` in configurable blocks, repacked into a
fit-local marker-major byte matrix, and file handles close after preparation.
Missing code 1 is omitted for AF estimation and mapped to the marker mean during
standardized access. Marker/sample orders follow `cls`, optional `rows`, and R
alignment metadata. There is no memory map and no repeated MCMC disk read.

## 6. Ownership and lifetime

| Object | Owner/lifetime | Mutability/copies |
|---|---|---|
| BED paths, `cls`, rows, AF, phenotype | binding/fit input | copied to native immutable metadata |
| blocked raw buffer | reader, one block | mutable and reused |
| packed genotype matrix | fit | immutable; one copy shared by chains |
| marker maps/order | fit | immutable shared |
| phenotype/residual/effects/components | chain | mutable chain-local |
| scheduler vectors/buckets | chain | mutable chain-local |
| BayesC engine | chain | chain-local |
| BayesC normal/uniform distributions | OpenMP worker/process | persistent and unsafe |
| BayesR/BayesRC engine/distributions | chain | fit-bounded and safe |

## 7. Scheduler-control inventory

BED BayesC and BayesR expose full-sweep frequency, null-skip base/max,
candidate threshold/lifetime and burn-in-only skipping. Defaults are 10, 50,
200, 0.001, 20 and false through `stblr_bed()`. Validation matches field
ranges. There is no neighbor/friend wake-up. BayesRC ignores these controls and
documents full sweeps. The sparse experimental BayesC route instead uses
`null_update_prob`, so its semantics are implementation-specific.

## 8. Mutable scheduler state

Adaptive chains own `scheduled_at`, due buckets, active/candidate flags and
lists, list-membership flags and last-interest/update state. State initializes
per chain, changes after marker updates, and is periodically compacted. It is
internal and no event-counter diagnostics are returned. BayesRC has no adaptive
scheduler state.

## 9. Traversal semantics

Iteration zero is a full sweep because iteration modulo frequency is zero.
Full sweeps use marker order. Non-full BED BayesC/BayesR traversal is active,
candidates, then due buckets. Skipped markers are not decoded or sampled and
consume no marker-update RNG. BED is already packed in memory, so traversal
does not issue disk reads. Trait-chain task order is deterministic; worker
assignment can affect unsafe BayesC distribution state.

## 10. Cross-backend scheduler comparison

| Semantic | CSR scheduled | BED BayesC | BED BayesR | BED BayesRC | Decision |
|---|---|---|---|---|---|
| full sweeps/modulo zero | yes | identical | identical | every iteration | subset only |
| adaptive skip | yes | identical thresholds | identical thresholds | unsupported | shareable C/R subset |
| candidates/due | yes | same ordering | same ordering | unsupported | shareable C/R subset |
| neighbor wake-up | yes | absent | absent | absent | CSR-specific |
| genotype decoding | CSR operator | packed bytes | packed bytes | packed bytes | BED-specific |
| RNG ownership | chain | worker distribution | chain | chain | statistically critical split |

## 11. RNG ownership inventory

All routes use `std::mt19937` and deterministic trait/chain seed formulas.
Unscheduled sparse BayesC constructs marker distributions per call. Scheduled
single/multichain BayesC uses active `static thread_local` uniform and normal
distributions: engine state is chain-local, cached distribution state is
worker/process-persistent. BayesR and BayesRC construct uniform/normal objects
inside each chain; gamma and chi-square objects are local at parameterized draw
sites. BayesC is classified cross-fit and worker-assignment risk; BayesR/RC are
chain-local and safe.

## 12. Persistent-distribution investigation

The retained Phase 10A diagnostic proves that a reused normal distribution's
cached variate survives engine reconstruction/reseeding. Source audit matches
that exact unsafe pattern in both scheduled BED BayesC files. The compact
two-chain fixture is same-process repeatable on one core but differs between
one and two cores in sampled fields, localizing worker assignment sensitivity.

## 13. Fresh-process reproducibility

One fresh-process raw/formatted reference was frozen for each of BayesC,
BayesR and BayesRC. Opt-in child-process comparisons reproduce all three
exactly. BayesC fresh-process stability does not remove its reused-process and
worker-assignment risk.

## 14. Same-process reproducibility

BayesC `A; A` is exact for the compact fixed worker configuration, while its
core-assignment comparison is not. BayesR and BayesRC repeated and intervening
calls are exact after declared core/timing metadata normalization. No unstable
BayesC same-process trajectory is treated as canonical.

## 15. Core/thread reproducibility

BayesR and BayesRC pass normalized `1,2,2,1`. BayesC one-core versus two-core
sampled output differs with identical trait-chain seeds, consistent with
worker-persistent distribution ownership. Explicit chain seeds are not exposed
by `stblr_bed()`; deterministic internal seed mapping is documented instead.

## 16. Frozen references

BayesC: 1 raw + 1 formatted fresh-process audit reference. BayesR: 1 raw + 1
formatted deterministic reference. BayesRC: 1 raw + 1 formatted deterministic
reference. Fixtures cover binary, mixture, annotation-component, multichain,
core, seed, scaling and packed decoding behavior. Generation is a manual tool.

## 17. Reduction and parity results

Adaptive scheduled and sparse BayesC do not have an exact supported public
reduction because the sparse route uses probability-based visits and different
distribution construction. BayesC/BayesR share decoding and scheduler shape
but not statistical formulas. BayesRC reduces to fixed-pi BayesR in the
existing fixed-annotation test, but adaptive scheduling is absent. No exact
CSR/BED parity is claimed because likelihood, operator and RNG order differ.

## 18. Typed audit contracts

`blr_genotype_backend_audit_types.h` defines packed source, scaling, scheduler
kind, chain/RNG ownership, execution audit and result vocabularies plus
field-specific validation. It is binding-neutral and deliberately does not
activate production execution.

## 19. Shared scheduler-contract decision

`ScheduledSweepControl`, `NullSkipControl` and `CandidateControl` are reusable
as a semantic subset for BED BayesC/BayesR after RNG correction.
`NeighborWakeupControl` is not applicable. `ScheduledExecutionControl` requires
a backend-specific genotype/I/O and RNG-ownership wrapper rather than reuse
unchanged. BayesRC requires a separate full-sweep contract.

## 20. Validation bridges

No Rcpp bridge was needed: the audit header is source-validated and production
R validation already covers BED paths, magic/mode, dimensions, rows, seeds,
annotations and controls. Consequently no new bridge can invoke a sampler,
decode a file or consume RNG, and generated wrappers remain unchanged.

## 21. Production behavior statement

All five audited production sources are byte-identical to the starting commit.
No sampler, decoder, signature, route, schema or formatter changed.

## 22. Runtime, memory, and I/O baseline

`Rscript tools/benchmarks/blr_phase11a_individual_bed_audit.R` records tiny
BayesC/BayesR/BayesRC chain/core cases and 2,000-marker/200-sample BayesC/BayesR
cases, five repetitions after warm-up. The tiny BED was 7 bytes and the
moderate BED was 100,003 bytes. Moderate BayesC times were
`0.11, 0.11, 0.09, 0.08, 0.09` seconds (median 0.09, range 0.03); moderate
BayesR times were `0.10, 0.10, 0.11, 0.10, 0.09` seconds (median 0.10, range
0.02). Tiny medians were 0.01--0.02 seconds and were dominated by startup,
clock resolution and cache effects, so no comparative speed claim is made.

The benchmark reports individual times, mean/median/range, file sizes and
whole-process RSS after completed fits. The latter was about 122--127 MB for
tiny cases, 135.5 MB after moderate BayesC, and 127.2 MB after moderate BayesR;
it is explicitly **not peak memory**. The first access includes open/read/pack
work; repetitions may benefit from the operating-system page cache. Process
memory, BED file size and I/O/cache effects are kept distinct, and BayesC,
BayesR and BayesRC are not treated as equivalent workloads.

## 23. Protected backends

Canonical CSR, scheduled CSR, block-eigen, multivariate, public route,
generated-wrapper and namespace hashes remain unchanged.

## 24. Tests

The final full suite passed 4,690 expectations with zero failures, zero
warnings, zero errors and one deliberate opt-in fresh-process skip. The Phase
11A focused file contributed 68 expectations: 67 pass and one opt-in test skips
in the ordinary run. With `SBLR_RUN_PHASE11A_FRESH=true`, the fresh-process raw
and formatted comparisons for all three models also pass (73 expectations,
no skips). Contract, source/RNG, same-process, core-order, scheduler, reduction,
individual/BED focused and canonical-protection checks all pass.

## 25. Confirmed defects and risks

| Backend | Source/reproduction | Consequence | Recommendation |
|---|---|---|---|
| scheduled BED BayesC single/multichain | active worker `static thread_local` normal/uniform; one/two-core fixture differs | sampled trajectory depends on worker history/assignment | correct chain/fit-owned distributions before migration |
| BED page-cache timing | blocked reader and repeated benchmark fits | warm/cold I/O timings differ, not statistical | report separately; do not infer sampler speed |

## 26. Recommended Phase 11B

> correct fit-local and chain-local RNG/distribution ownership in the affected individual-level or packed-BED backend before any numerical-core migration.

## 27. Readiness marker

PHASE 11A COMPLETE — INDIVIDUAL AND PACKED-BED BACKEND AUDIT READY
