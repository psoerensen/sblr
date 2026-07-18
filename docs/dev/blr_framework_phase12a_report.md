# Unified BLR Framework Phase 12A Report

## 1. Executive summary

Cross-cutting architecture, documentation, validation, CI, and measurement
hardening was completed without changing posterior targets, MCMC traversal,
RNG draw sites, genotype decoding, public APIs, or schemas.

## 2. Repository baseline

Branch: `master`. Starting and Phase 11D commit: `f27ace3` (`Canonicalize
scheduled packed-BED BayesC`). Initial status was clean and `git diff --check`
passed. The toolchain is R 4.4.1, Rtools44, GNU C++17 and OpenMP. The Phase 11D
suite baseline was 4,891 passes, zero failures/errors/warnings and two opt-in
skips. No independently visible workflow existed initially.

The recovered Phase 10D scheduled-CSR baseline completed all seven workloads;
the representative moderate dense workload recorded 0.24, 0.17, 0.18, 0.20,
and 0.19 seconds (mean 0.196, median 0.19, range 0.07) with completed-fit RSS
approximately 128--139 MiB. The recovered Phase 11D packed-BED baseline recorded
0.19, 0.47, 0.76, 0.33, and 0.18 seconds for moderate dense 1x1 (mean 0.386,
median 0.33) and 0.19, 0.17, 0.25, 0.19, and 0.22 seconds for moderate
aggressive 2x2 (mean 0.204, median 0.19), with completed-fit RSS approximately
123--138 MiB. These short Windows debug timings are noisy behavior baselines,
not speed claims, and completed-fit RSS is not peak RSS.

## 3. Hardening inventory

| Finding | Location/current behavior | Risk and decision | Tests |
|---|---|---|---|
| stale plan header | implementation plan named Phase 9C | misleading; synchronize now | status text |
| stale matrix header | matrix named Phase 4 | misleading; synchronize now | status text |
| stale summaries | historical “ready/migration” text | retain history, label point-in-time | historical label |
| control bytes | plan/reduction matrix form-feed, backspace and tab-corrupted notation | malformed math; repair | byte scanner/LaTeX tokens |
| brittle counts | framework `length(gregexpr())` | absent pattern counted as one; shared helper | zero/one/multiple |
| sweep zero | shared validator rejected zero while cores/public adapters accepted it | contract disagreement; accept zero | -1/0/1/large |
| CSR `pi` mean | converter traversed retained `pis` | posterior aggregation misplaced; move unchanged arithmetic to typed core result | exact fixtures/layer test |
| CSR core output | unconditional `std::cout` | reusable-core side effect; remove | source boundary |
| BED core output | conditional progress output | callback refactor risks public progress semantics; narrowly defer | remains guarded by `progress_every` |
| commented snapshots | canonical CSR/BED binding sources contained 1,034/7,437 inactive lines | maintenance/source-test risk; remove, Git retains history | snapshot guard |
| CI absent | no workflows | regressions not independently visible; add fast/manual tiers | workflow structure |
| peak RSS absent | completed-fit RSS only | cannot observe peak; add sampled child-process tool | opt-in smoke |

## 4. Files changed

Documentation/status files, three generated Rd files, shared and phase-specific
tests, scheduled-control/type/core/binding files at approved hardening seams,
two canonical binding sources for inactive-tail deletion, two GitHub Actions
workflows, and peak-RSS tooling changed. No backend migration was started.

## 5. Current-status documentation

The plan and matrix now identify Phase 11D as the latest completed model phase,
list all canonical ordinary and scheduled routes, and list audited/protected
noncanonical BED, block-eigen and multivariate routes. Earlier phase statements
are explicitly historical snapshots.

## 6. Markdown/control-character repair

Four form-feed-corrupted `\frac` sequences, one backspace-corrupted `\bmod`,
and tab-corrupted `\theta` text were repaired in the plan and reduction matrix.
Maintained `docs/dev/*.md` now permits only tab, LF and CR controls.

## 7. Structural-test helper

`helper-source-architecture.R` returns an empty integer vector for the
`gregexpr()` `-1` sentinel and consequently counts absent, single and duplicate
matches as zero, one and their true count. Fixed and regex modes are tested;
affected framework structural assertions use the helper.

## 8. Scheduled validation consistency

For scheduled controls already implementing the `<= 0` full-sweep branch,
zero is accepted and means a full sweep every iteration, positive values retain
periodic behavior, and negative values reject. Shared, CSR, and BED typed/public
validation now agree. Numerical branch behavior was not changed.

## 9. Aggregation and conversion responsibilities

Scheduled-CSR retained-sample posterior mean inclusion probability is now
computed in the binding-neutral typed numerical result using the identical
loop order, indices, thinning rule, arithmetic and final-value fallback. The
converter only shapes the two-column R representation. Timing mean/max remain
converter-side execution-metadata presentation. Packed-BED posterior and CPO
aggregation was already native and required no change.

## 10. Diagnostic-output boundary

Unconditional scheduled-CSR reusable-core configuration/task `std::cout` was
removed; it carried no numerical state and was not part of returned schemas.
No marker loop logs. Packed-BED per-chain progress output is retained only under
the existing explicit `progress_every > 0` guard: moving it safely requires a
future binding-neutral callback contract and was narrowly deferred to avoid
changing progress timing or worker synchronization.

## 11. Historical source cleanup

The inactive duplicate implementation beginning after the active scheduled-CSR
binding and five successive inactive packed-BED snapshots after the active
public binding were deleted. They were uncompiled, unused by runtime tooling,
and remain available at `f27ace3`. Short semantic comments remain allowed.

## 12. Fast CI

`.github/workflows/blr-framework.yml` runs on master pushes, pull requests and
manual dispatch on Ubuntu. It sets up cached R dependencies, verifies generated
attributes by compiling/loading, runs fast deterministic framework/schema and
core-order coverage, and runs a scoped package check.

## 13. Extended validation workflow

`.github/workflows/blr-framework-extended.yml` is manual. It enables fresh-
process matrices, compiles, runs the full suite and package check, runs peak-RSS
smoke measurement, and uploads check artifacts. Costly benchmarks are excluded
from the fast gate.

## 14. Peak-memory measurement

`tools/benchmarks/measure_peak_rss.R` starts a child `Rscript` with `processx`,
samples root plus recursive descendants through `ps`, and records sampled peak
RSS, final sampled RSS, interval, count, method, platform, elapsed time, output
and exit status. Linux uses RSS; Windows uses the equivalent working-set field.
The smoke workload allocates 16 MiB and sleeps long enough for interval
sampling. Completed-fit RSS remains separately labelled and is not called peak.
Canonical benchmark commands can be supplied as the child command/arguments.
The final Windows smoke exited successfully with sampled peak/final working set
76,533,760 bytes, two samples at 0.02-second intervals, and method `ps root plus
recursive descendants`.

## 15. Numerical reference protection

All permanent exact-reference families passed: ordinary CSR BayesC, BayesR,
SBayesRC, fixed-prior BayesC, group BayesC, learned-annotation BayesC,
scheduled CSR BayesC, public scheduled packed-BED BayesC, and packed-BED BayesR
and BayesRC audit references. The corrected scheduled raw/formatted families
remain 3/3 exact for CSR and 3/3 exact for packed BED. No fixture was regenerated.

## 16. Reproducibility

Repeated-call, intervening-fit, one/two-core, reversed-core-order, explicit-seed
and fresh/reused-process results remain exact. The enabled Phase 10A/10B and
11A/11B fresh-process matrix passed 193 expectations without failures, warnings,
or skips.

## 17. Public API and schema

Arguments, defaults, routes, native public signatures, generated Rcpp exports,
`NAMESPACE`, `stblr_raw_v1`, formatted fits, and actual R `NULL` semantics are
unchanged. Rd text only documents the already-active zero policy.

## 18. Protected backends

The Git/formula audit found native changes only at approved hardening seams:
shared validation, scheduled-CSR posterior-summary placement and console
removal, plus inactive canonical snapshot deletion. MCMC loops, formulas, RNG
sites, seed mapping and genotype access remain unchanged. Ordinary CSR,
packed-BED BayesR/BayesRC and experimental BayesC, block-eigen, multivariate,
generated wrappers, public routes and `NAMESPACE` remain protected.

## 19. Tests

The focused affected Phase 10C3/10D/11D/12A run passed 175 expectations with
one opt-in memory skip. Phase 12A with peak measurement enabled passed 53/53.
The enabled fresh-process matrix passed 193/193. The final full suite passed
4,941 expectations with zero failures or warnings and three explicitly opt-in
skips (Phase 11A fresh, Phase 11B fresh, and peak RSS); those opt-in paths were
also run separately and passed. The only warning outside testthat was the benign
notice that `testthat` was built under R 4.4.3.

## 20. CI validation

Workflow presence, triggers and command structure are tested locally without
network access. Actual GitHub-hosted action execution remains external
verification after publication.

## 21. Deviations and blockers

The packed-BED progress callback boundary is narrowly deferred as documented.
Peak RSS is sampled and may miss sub-interval spikes; the interval is always
reported. No statistical or completion blocker remains.

## 22. Recommended next phase

> begin a typed-boundary and deterministic-reference migration of packed-BED BayesR, reusing only the scheduled controls and genotype ownership contracts proven semantically compatible during Phase 11A.

## 23. Readiness marker

PHASE 12A COMPLETE — CROSS-CUTTING ARCHITECTURE HARDENING ESTABLISHED
