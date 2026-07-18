# Unified BLR Framework Phase 15A Report

## 1. Executive summary

Canonical packed-BED BayesC, BayesR, and BayesRC consolidation opportunities
were audited without modifying production numerical implementations. Exact task
mapping and seed resolution are the strongest future sharing candidates;
model-specific numerical and schema boundaries remain separate.

## 2. Repository baseline

Branch `master` started clean at `95a7524` (Phase 14E). R 4.4.1 with Rtools44,
Rcpp 1.1.1, and sblr 0.1.0 are the local toolchain. The Phase 14E full baseline
was 5,428 passes, zero failures/warnings, and eight opt-in skips. Canonical
benchmark baselines are Phase 11D BayesC, Phase 13E BayesR, and Phase 14E
BayesRC; each records warm-up timings, completed-fit RSS and I/O limitations.

## 3. Canonical family inventory

| Model | Public/native path | View/context/core | Result/aggregation/converter | Fixtures/benchmark |
|---|---|---|---|---|
| BayesC | `stblr_bed(bayesc)` / scheduled-chains adapter | `BedPackedGenotypeView`, `BedScheduledBayesCChainExecutionContext`, `run_bed_scheduled_bayesc_chain` | typed chain/aggregate, `aggregate_bed_scheduled_bayesc_results`, `stblr_bed_scheduled_bayesc_result_to_raw` | Phase 11B / 11D |
| BayesR | `stblr_bed(bayesr)` / BayesR adapter | `BedBayesRPackedGenotypeView`, `BedBayesRChainExecutionContext`, `run_bed_bayesr_chain` | typed chain/aggregate, `aggregate_bed_bayesr_results`, `stblr_bed_bayesr_result_to_raw` | Phase 13A / 13E |
| BayesRC | `stblr_bed(bayesrc)` / BayesRC adapter | `BedBayesRCPackedGenotypeView`, `BedBayesRCChainExecutionContext`, `run_bed_bayesrc_chain` | typed chain/aggregate, `aggregate_bed_bayesrc_results`, diagnostics, `stblr_bed_bayesrc_result_to_raw` | Phase 14A / 14E |

Each has one route, dispatch, core, aggregate path, converter, and RNG path.

## 4. Files changed

Only this report, the plan/matrix, Phase 15A audit tests, and the audit benchmark
are changed. No production source is edited.

## 5. Packed-genotype representations

BayesC and BayesR views have identical fields: storage reference, byte pointer,
packed size, marker/sample counts, bytes per marker, and stride. Both expose
decoded SNP-major marker storage; lookup/scaling remains in model marker maps.
BayesRC borrows the same lifetime class but exposes only storage, counts, and
bytes per marker; access is through its storage API. All preserve missing-value,
centering, scaling, and marker order established by their adapters. Decision:
BayesC/R are representation-equivalent; BayesRC is ownership-equivalent only.

## 6. Storage ownership and lifetimes

All decode once per fit, retain fit-owned immutable storage, borrow it in each
chain, make no full per-chain copy, and perform no MCMC-time disk access. Marker
maps and phenotypes are fit-owned; chain state and RNG are private; R names and
schema metadata are binding-owned. A production ownership type would add little;
reuse the convention and tests.

## 7. Task mapping

All use `task_count = trait_count * chain_count`, `chain = job / trait_count`,
`trait = job % trait_count`, job-major result slots, static scheduling, and
trait-then-chain aggregation derived from that layout. Decision: reuse unchanged
as a small constexpr/index helper, protected by ordering tests.

## 8. Seed policies

All resolve before context construction using an unsigned result of
`seed + 1000003*(trait+1) + 9176*(chain+1)`, then store the resolved logical-chain
seed in the context. Worker identity contributes nothing. Decision: reuse
unchanged only with exactly preserved casts and overflow behavior.

## 9. OpenMP dispatch

All preallocate one typed result per job and use `schedule(static)`. BayesC/R
render detailed completion information after workers; BayesR also replays typed
progress events. BayesRC throws failed task errors before aggregation. A generic
dispatch template risks obscuring failure and progress behavior; reuse through a
narrow common task-index wrapper only.

## 10. Failure vocabulary

All chain results carry `failed` and `error`; BayesC/R count and render task,
trait, and chain before translating the first failure, while BayesRC stores
aggregate messages and translates its first failure. Public failure schemas
differ. A `task/trait/chain/message` envelope is reusable through a narrow common
wrapper, but conversion remains model-specific.

## 11. Timing vocabulary

BayesC/R retain per-chain seconds and aggregate mean/max by trait. BayesRC does
not expose the same chain timing vocabulary in its current result. Benchmark
elapsed time and completed-fit RSS are presentation measurements. Reuse
convention only; numerical timing fields are not interchangeable.

## 12. Progress handling

BayesC currently renders progress from its chain implementation, BayesR captures
binding-neutral events and deterministically renders after dispatch, and BayesRC
has no equivalent progress-event payload. Only an event envelope is conceptually
shareable; payload and presence remain model-specific. Introducing progress into
BayesRC or changing BayesC rendering is outside this audit.

## 13. Optional genotype diagnostics

BayesC computes requested `wy`/residual scores through its aggregation context;
BayesR and BayesRC fill analogous diagnostics after typed aggregation. All
preserve actual `NULL` when unrequested and require genotype/phenotype access.
Reuse the post-aggregation boundary convention only; do not move genotype into
aggregate results or force one implementation.

## 14. Converter metadata

Marker/trait names, model/backend identifiers, MCMC controls, input metadata,
schema version/classes, dimnames, and actual-NULL policy are common concepts.
BayesC scheduler/prior metadata, BayesR component labels/scales, and BayesRC
annotation names/alpha priors are model-specific. A narrow common metadata
envelope is plausible; converters and schemas remain separate.

## 15. Orientation and helper functions

Trait/chain scalar matrices and marker-by-trait matrices use common conventions.
Component arrays differ: BayesR has global ordered mixture probabilities;
BayesRC has marker priors and annotation/stick axes. Dimname attachment is
binding-specific. Reuse helper only for exact two-dimensional common shapes;
component/annotation orientation is unsafe to consolidate.

## 16. Aggregate-result conventions

Dimensions, marker means/PIP, variance traces, CPO, retained counts, timing,
failures, and optional diagnostics are vocabulary-level overlaps. BayesC binary
probability and scheduler semantics, BayesR global components, and BayesRC
marker priors/alpha/step variances are distinct. Retain model-specific aggregate
types; a base class or inheritance is harmful.

## 17. Reference harnesses

All fixture families separate raw/formatted objects, normalize declared timing
metadata only, support fresh processes, and report exact structural differences.
Reuse helper only for recursive first-difference reporting and process launching;
retain model-specific capture and normalization functions and never normalize
sampled values.

## 18. Architecture-test helpers

Zero-safe occurrence counting, active-source filtering, forbidden-token scans,
include-user checks, protected hashes, and file-I/O/RNG scans are exactly reusable
test helpers. Semantic scheduler, component, stick, and annotation checks remain
model-specific.

## 19. Benchmark reporting

R/package/compiler versions, samples, markers, iterations, burn-in, thinning,
chains, cores, individual/summary times, completed-fit RSS, optional sampled
peak RSS, BED size, and page-cache caveats form the common convention. BayesC
scheduler, BayesR mixture, and BayesRC annotation/update controls remain required
model-specific metadata. Workloads must not be speed-ranked across models.

## 20. Shared-infrastructure decision matrix

| Candidate | BayesC | BayesR | BayesRC | Equivalence | Decision | Risk | Phase |
|---|---|---|---|---|---|---|---|
| packed view | full byte view | full byte view | reduced storage view | C/R exact; RC ownership only | reuse through narrow common wrapper | access/stride | 15B optional |
| ownership | fit immutable | fit immutable | fit immutable | exact contract | reuse convention only | lifetime | tests/docs |
| task mapping | chain-major | chain-major | chain-major | exact | reuse unchanged | ordering | 15B |
| seed resolution | same constants/cast | same | same | exact | reuse unchanged | overflow/RNG | 15B |
| OpenMP dispatch | static + model reporting | static + events | static + early throw | structural only | reuse helper only | failure/order | later |
| failure vocabulary | failed/error | failed/error | failed/error + aggregate strings | envelope equivalent | reuse through narrow common wrapper | schema | 15B optional |
| timing | chain mean/max | chain mean/max | different exposure | partial | reuse convention only | meaning | none |
| progress | core rendering | typed events | none | different | retain model-specific | thread/side effects | none |
| diagnostics | aggregation-aware | post aggregate | post aggregate | boundary only | reuse convention only | genotype/schema | none |
| converter metadata | common names/schema concepts | components | annotations | partial | reuse through narrow common wrapper | field order | later |
| orientation | common 2-D shapes | component axes | annotation/stick axes | partial | reuse helper only | transposition | later |
| reference harness | exact recursive compare | exact | exact | helper equivalent | reuse helper only | over-normalization | tests |
| architecture tests | source scans/hashes | same | same | exact helper layer | reuse unchanged | shallow checks | tests |
| benchmark metadata | common envelope + scheduler | mixture | annotation | convention equivalent | reuse convention only | false ranking | tooling |
| contexts/results/aggregates | binary | mixture | annotated mixture | statistically distinct | unsafe to consolidate | numerical/schema | never |

## 21. Candidate audit contracts

No audit-only C++ header was added: it would duplicate documentation without
proving production reuse. The report and source-matched tests are sufficient.

## 22. Direct equivalence tests

Tests directly protect identical task formulas, seed constants, static dispatch,
immutable borrowing, absence of core file I/O and worker-derived seeds, common
schema class/NULL conventions, distinct converters, and starting-commit hashes.

## 23. Model-specific contracts protected

BayesC retains binary inclusion, adaptive scheduling, null skipping and candidate
state. BayesR retains ordered global mixtures, adaptive scheduling, and custom
inverse-CDF sampling. BayesRC retains marker-specific annotation probit sticks,
latent/sequential alpha updates, and full sweeps without scheduler state.

## 24. Exact reference results

BayesC Phase 11B, BayesR Phase 13A, and BayesRC Phase 14A each retain three raw
and three formatted exact configurations: 9/9 raw and 9/9 formatted in total.

## 25. Reproducibility

Permanent same-process, fresh-process, core-order, worker, chain-count, and
intervening-fit matrices remain exact for all families. The fixed-alpha
BayesRC-to-fixed-pi BayesR reduction remains exact. The explicitly enabled
BayesC/BayesR/BayesRC fresh-process selection passed 157 assertions with zero
failures, warnings, or skips.

## 26. Public API and schema

Arguments, defaults, routes, native signatures, raw/formatted schemas, ordering,
actual `NULL`, wrappers, and `NAMESPACE` are unchanged.

## 27. Protected backends

Canonical and experimental packed-BED production, CSR, block-eigen,
multivariate, wrapper, and namespace files are unchanged from `95a7524`.

## 28. Performance, memory, and I/O baselines

Phase 11D, 13E, and 14E remain the independent canonical baselines. Completed-fit
RSS is not peak RSS; sampled peaks are opt-in; page cache affects BED reads. No
cross-model speed claim is made and no production I/O path changed. In this run,
BayesC tiny medians were 0.01--0.02 seconds (moderate opt-in workloads 0.24--0.25
seconds) with completed-fit RSS roughly 122--129 MB; BayesR tiny medians were
0.01--0.02 seconds with RSS roughly 121--122 MB; BayesRC tiny medians were
0.02--0.03 seconds with RSS roughly 122--131 MB. First-configuration warm/cache
outliers were retained rather than interpreted as speed. The Phase
15A audit benchmark successfully resolved all three canonical scripts and their
model-specific fixture/control metadata. The new focused audit passed 69
assertions. Compilation/loading succeeded; the final full suite passed 5,497
assertions with zero failures/warnings and eight documented opt-in skips.

## 29. Risks

| Source | Semantic/numerical/RNG/ordering/schema risk | Safeguard |
|---|---|---|
| common view | BayesRC access differs | preserve wrapper and reference hashes |
| task/seed helper | cast or order changes trajectories | exact formula/type tests and all fixtures |
| dispatch/failure | reporting and exception order differ | keep adapter policy explicit |
| metadata/orientation | field or axis reordering | schema snapshots and converter-specific tests |
| reference helpers | excessive normalization | normalize declared metadata only |
| benchmark convention | invalid speed ranking | model-specific workloads and labels |

## 30. Recommended Phase 15B

Introduce only the proven representation-equivalent packed-BED common
infrastructure—immutable genotype-view vocabulary for compatible views,
task-index mapping, and narrowly compatible failure envelopes—behind permanent
references while leaving all model-specific contexts, numerical cores,
schedulers, probability policies, aggregators, and converters unchanged.

## 31. Readiness marker

PHASE 15A COMPLETE — PACKED-BED FAMILY CONSOLIDATION BOUNDARIES AUDITED
