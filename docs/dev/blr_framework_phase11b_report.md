# Unified BLR Framework Phase 11B report

## 1. Executive summary

Scheduled packed-BED BayesC RNG/distribution ownership was corrected without
changing genotype decoding, scheduler rules or statistical formulas.

## 2. Repository baseline

The clean baseline was `master` at Phase 11A commit `e2035ab`. Initial status
and `git diff --check` were clean. R 4.4.1, Rtools44/GCC and OpenMP were used;
the Phase 11A suite passed 4,690 expectations.

## 3. Confirmed Phase 11A defect

Both scheduled BayesC sources owned normal/uniform distributions as
`static thread_local`. One/two-core execution first differed in numeric
`raw$marker$bm`, followed by marker summaries, variance traces and log-CPO.
Engines were chain-seeded but cached distribution state belonged to workers.

## 4. Files changed

The two scheduled BayesC sources received ownership-only edits. A binding-neutral
chain RNG header, corrected fixtures/tests, benchmark and manual fixture tool
were added. Audit contracts, Phase 11A structural checks, implementation plan
and capability matrix were updated. BayesR/BayesRC and decoder sources were not.

## 5. Corrected ownership model

Each logical trait-chain constructs one `BedScheduledBayesCChainRng` after the
existing seed is computed. It owns one `std::mt19937`, standard normal and
`[0,1)` uniform distribution for exactly one chain execution. Worker identity
does not participate in ownership.

## 6. Distribution inventory

Normal/uniform marker distributions are chain-persistent only. Uniform-integer
jitter remains chain-local beside the engine. Variable-parameter chi-square,
gamma and Beta-equivalent draws remain locally constructed at unchanged update
sites. No discrete distribution or worker-indexed RNG container was introduced.

## 7. Intentional numerical changes and preserved behavior

Intentionally changed: inherited worker-owned cached distribution state and
cross-chain/cross-fit cache sharing. Unchanged: engine, seed formulas, logical
draw sites, branch order, scheduler-dependent draws, decoding, formulas,
aggregation and public output structure.

## 8. Single-chain route correction

The existing trait seed constructs one chain RNG inside the OpenMP trait task.
The marker helper receives it by reference; other helpers retain an engine
reference to the same chain-owned engine.

## 9. Multichain route correction

Each trait-chain task computes the existing seed, constructs its chain RNG and
passes it to marker updates. Sequential or reassigned tasks cannot share caches.

## 10. Genotype and I/O preservation

BED magic/mode checks, SNP-major blocked reads, marker/sample order, missing
mapping, scaling, fit-local packed storage and immutable sharing are unchanged.
No per-chain genotype copy, new buffer or MCMC-time disk read was added.

## 11. Scheduler preservation

Full-sweep timing, null-skip growth/reset, candidate threshold/lifetime, due
buckets, burn-in-only behavior and active/candidate/due traversal are unchanged.
Skipped markers still perform neither marker RNG nor genotype update work.

## 12. Post-correction references

Three raw and three formatted references cover single 1x1, multichain 2x1 and
multichain 2x2. Metadata records `bed_scheduled_bayesc_chain_rng_v1`, baseline,
dimensions, seed, controls, iterations and schema. Phase 11A BayesC remains
unchanged as historical pre-correction evidence.

## 13. Same-process reproducibility

`A; A`, `A; B; A`, intervening BayesR, intervening BayesRC and one/two/one
chain-count sequences are exact after elapsed/core metadata normalization.

## 14. Core/thread reproducibility

Normalized `1,2,2,1` is exact. Identical logical seeds are worker-independent.

## 15. Fresh-process reproducibility

The isolated-process 2x2 raw and formatted comparison is exact.

## 16. Single-chain versus multichain parity

The routes are intentionally nonidentical: the multichain seed adds
`9176 * (chain + 1)` and aggregation metadata differs. Equality was not forced.

## 17. BayesR and BayesRC protection

Both sources are byte-identical to `e2035ab`; Phase 11A references remain exact.

## 18. Public API and schema

Arguments, defaults, routing, signatures, wrappers, `NAMESPACE`, raw/formatted
schemas and actual `NULL` behavior are unchanged.

## 19. Performance, memory, and I/O

`Rscript tools/benchmarks/blr_phase11b_bed_bayesc_rng.R` used tiny and
2,000-marker/200-sample workloads with five timed repetitions after warm-up.
Moderate medians were 0.10 seconds (single 1x1), 0.17 seconds (2x1) and 0.16
seconds (2x2), with ranges 0.03, 0.03 and 0.08 seconds. Phase 11A's comparable
public 1x1 median was 0.09 seconds; these short Windows timings do not establish
a material regression or improvement. The 100,003-byte BED and blocked-read
strategy are unchanged. Whole-process completed-fit RSS was roughly 125--140
MB and is not peak memory; page-cache effects remain a stated limitation.

## 20. Tests

The Phase 11B ordinary run passed 34 expectations with one opt-in fresh-process
skip; the opt-in run passed all 36. The final suite passed 4,724 expectations
with zero failures, warnings or errors and two deliberate fresh-process skips
(Phase 11A and Phase 11B). Corrected references, same-process/core-order,
protected BED, canonical and schema checks all pass.

## 21. Deviations and blockers

No blocker. The old worker-cache-dependent Phase 11A trajectory is intentionally
not reproduced.

## 22. Recommended Phase 11C

> establish typed execution contracts and mechanically extract the corrected scheduled packed-BED BayesC implementation, beginning with the route whose execution boundary and deterministic reference coverage are clearest, while retaining BayesR and BayesRC production paths unchanged.

## 23. Readiness marker

PHASE 11B COMPLETE — SCHEDULED PACKED-BED BAYESC RNG OWNERSHIP CORRECTED
