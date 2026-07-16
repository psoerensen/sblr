# Unified BLR Framework Phase 9D3 Report

## 1. Executive summary

Group CSR BayesC migration is complete with one typed numerical core and one
named binding-layer result converter. Numerical behavior and public schemas are
unchanged.

## 2. Repository baseline

- Branch: `master`.
- Starting commit and Phase 9D2 commit: `4b5d3c1` (`Activate typed group BayesC execution boundary`).
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1 UCRT, Rtools44 GCC, C++17 and OpenMP on x86_64 Windows.
- Baseline: native compile/load passed; full suite passed 4,149 tests with no failures or skips.
- Contemporaneous Phase 9A group benchmark medians were 0.11 s (one chain/one
  core), 0.19 s (two chains/one core), and 0.19 s (retained two chains/two
  cores); completed-fit RSS was 128.7--137.4 MiB. These are noisy whole-call
  measurements rather than isolated core timings.

## 3. Phase 9D2 structure inventory

| Item | Classification | Resolution |
|---|---|---|
| `CsrGroupBayesCExecutionContext` and policy view | retain permanently | Borrowed execution boundary remains active. |
| `run_csr_group_bayesc()` and typed result | retain permanently | Sole numerical path remains unchanged. |
| `cpg_group_raw_v1` inline-style converter | centralize now | Renamed and typed as `stblr_csr_group_bayesc_result_to_raw()`. |
| field construction and conversion helpers | retain permanently | Kept inside the single binding converter/helper boundary. |
| aliases used by wrapper aggregation | retain permanently | Required to preserve the existing chain schema and order. |
| two converter call branches | retain permanently | They are the mutually exclusive retained-chain and ordinary aggregate exits. |
| inline-converter assertions in Phase 9D1/9D2 | remove now | Replaced with named-converter structural assertions. |
| numerical context/result fields | defer with justification | All remain consumed by execution or conversion; no proven dead fields. |
| `std::cout` diagnostics | retain permanently | Binding-neutral and outside marker updates. |

No duplicate conversion fragment, route selector, fallback, compatibility
logging, unused context field, or lexical execution path was found.

## 4. Files changed

- `src/st_cpg_omp_csr_group.cpp`: named typed-result converter and its two
  mutually exclusive aggregate exits.
- Phase 9D1 and Phase 9D2 tests: permanent named-converter protection.
- `tests/testthat/test-blr-framework-phase9d3.R`: closure architecture, exact
  reference, reproducibility, policy, and restriction checks.
- `tools/benchmarks/blr_phase9d3_csr_group_bayesc.R`: reproducible timing and
  completed-fit RSS matrix.
- Implementation plan and capability matrix: migrated/ready status.
- This report.

## 5. Final execution path

The public group route performs existing R decoding, alignment and validation,
prepares group policy and CSR state, builds the typed context, calls
`run_csr_group_bayesc()`, aggregates chains once, converts the typed aggregate
once, and returns unchanged `stblr_raw_v1` for the unchanged formatter.

## 6. Centralized result converter

`stblr_csr_group_bayesc_result_to_raw(const
CsrGroupBayesCExecutionResult&, ...)` is the sole named group result converter.
It owns binding-only construction of names, dimensions, metadata, diagnostics,
optional chain payloads, LD-swap and selection fields. Storage types, field
order, classes, actual `NULL`, marker/trait/group/chain order, and schema are
unchanged.

## 7. Wrapper-level multichain aggregation

Each core invocation returns native raw numerical state. The existing wrapper
combines chain marker summaries, group means/standard deviations, traces and
optional chain payloads in input chain order, then calls the converter once on
the selected aggregate branch. `keep_chains` remains unchanged. Aggregation is
not duplicated, chains are not converted twice, and normalization is not
reapplied. Moving this boundary could alter chain ordering and schema.

## 8. Native adapter boundary

The adapter is limited to binding decode/validation, alignment, group and CSR
preparation, typed-context construction, core invocation, the single converter,
wrapper aggregation, and exception translation. MCMC, marker/RNG draws, group
updates, normalization, global updates, residual/variance updates and posterior
accumulation remain absent from the adapter.

## 9. Numerical core

One `run_csr_group_bayesc()` implementation and one active marker loop remain
in `src/blr_csr_group_bayesc_core_impl.h`. Phase 9D3 did not modify that header
or any numerical statement.

## 10. Group-policy preservation

The borrowed immutable zero-based mapping, mapping length, nonempty group count,
explicit group order, probabilities, variance multipliers, priors, update flags
and timing are unchanged. Both `normalize_group_vb` modes retain their existing
formula/timing, as do global `pi`/`B` interactions. Invalid mapping, empty-group,
dimension, equal-sample-size/shared-`ww`, normalization-control and unsupported
trait behavior remain unchanged.

## 11. Logging decision

Existing `std::cout` messages are retained. They are outside the marker hot
loop, consume no RNG, do not change control flow or returned state, preserve
message content, and do not add worker-loop writes. They keep the core free of
Rcpp binding calls.

## 12. Exact frozen references

Phase 9D3 tests match all group raw references 3/3 exactly and all formatted
references 3/3 exactly, including types, dimensions, names, classes, actual
`NULL`, order, schema, group state, variance outputs, diagnostics, LD-swap and
selection fields. Fixtures were not regenerated.

## 13. Reproducibility

Exact checks pass for repeated calls, one/two cores, reversed `1,2,2,1` core
order, intervening learned-annotation fits, explicit seeds covered by the Phase
9A/9D suites, retained/dropped chains, fixed/updated probabilities and
multipliers, both normalization modes, and disabled LD-swap.

## 14. Public API and schema

Public arguments, defaults, native signatures, routing, `NAMESPACE`, generated
wrappers, `stblr_raw_v1`, formatted fields and present-but-`NULL` behavior are
unchanged.

## 15. Protected backends

Source audits and framework tests protect canonical BayesC, BayesR, SBayesRC,
fixed-prior BayesC, learned-annotation BayesC, block-eigen, BED, scheduled and
multivariate implementations. No protected native file changed.

## 16. Performance and memory

Command: `Rscript tools/benchmarks/blr_phase9d3_csr_group_bayesc.R`. Each row
has one warm-up and five timed repetitions. The representative workload has
2,000 markers, one trait, two ordered nonempty groups (20/80 mapping), 100
iterations plus 25 burn-in iterations.

| Configuration | Times (s) | Mean / median / min / max (s) | Completed-fit RSS (MiB) |
|---|---|---|---:|
| fixed probabilities/multipliers, 1x1 | .09,.10,.09,.10,.09 | .094 / .09 / .09 / .10 | 136.1 |
| updated multipliers, normalized, 1x1 | .09,.09,.12,.09,.11 | .100 / .09 / .09 / .12 | 129.4 |
| updated multipliers, unnormalized, 1x1 | .09,.10,.10,.10,.09 | .096 / .10 / .09 / .10 | 128.8 |
| updated probabilities/multipliers, 1x1 | .08,.08,.08,.08,.06 | .076 / .08 / .06 / .08 | 129.4 |
| two chains/one core | .18,.21,.17,.19,.21 | .192 / .19 / .17 / .21 | 129.7 |
| retained two chains/two cores | .17,.15,.17,.17,.19 | .170 / .17 / .15 / .19 | 129.2 |

The tiny correctness row was timer-resolution dominated (median 0 s) and is
not used for regression claims. The directly comparable Phase 9A medians were
0.11, 0.19 and 0.19 s; post-migration values are 0.09, 0.19 and 0.17 s. RSS
ranges overlap. No material regression is evident, but timings are short and
order/cache effects remain. RSS is sampled after completed fits for the whole R
process, not interval-sampled peak memory.

## 17. Test results

- Phase 9A + Phase 9D1 + Phase 9D2 + Phase 9D3 focused migration suite: 171
  passed, 0 failed, 0 warned, 0 skipped (Phase 9D3 alone: 33 passed).
- Baseline full suite: 4,149 passed, 0 failed, 0 skipped.
- Final full suite, including group-focused, canonical-model,
  learned-annotation, schema and backend protection: 4,182 passed, 0 failed,
  0 warned, 0 skipped.

## 18. Deviations and blockers

The first benchmark tiny mapping assigned no marker to one declared group and
correctly failed existing validation; the benchmark mapping was corrected to
exercise two nonempty groups. Memory is completed-fit RSS rather than sampled
peak RSS. There are no implementation or validation blockers.

## 19. Recommended next phase

> canonicalize and stabilize group CSR BayesC, remove remaining migration-only wording or aliases, retain permanent exact fixtures, establish the post-migration baseline as canonical, and then begin the learned-annotation BayesC migration.

## 20. Readiness marker

PHASE 9D3 COMPLETE — GROUP BAYESC MIGRATED WITH BEHAVIOR PRESERVED
