# Unified BLR Framework Phase 9E Report

## 1. Executive summary

Group CSR BayesC is canonicalized and stabilized with one typed numerical core,
one active marker loop, one binding converter and one wrapper aggregation path.
No native numerical code changed.

## 2. Repository baseline

- Branch: `master`.
- Starting and Phase 9D3 commit: `14db955` (`Migrate group CSR BayesC to typed core`).
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1 UCRT, Rtools44 GCC, C++17 and OpenMP on x86_64 Windows.
- Baseline focused Phase 9A/9D1/9D2/9D3 suite: 171 passed.
- Baseline full suite: 4,182 passed, 0 failed, 0 warned, 0 skipped.
- Phase 9D3 pre-cleanup representative medians: 0.09 s fixed 1x1,
  0.09 s updated-multiplier 1x1, 0.09 s unnormalized 1x1, 0.08 s
  updated-probability 1x1, 0.18 s 2x1 and 0.19 s retained 2x2;
  completed-fit RSS was 128.4--136.1 MiB.

## 3. Phase 9D3 structure inventory

| Item | Classification | Phase 9E resolution |
|---|---|---|
| Typed context and `GroupBayesCPolicyView` | retain permanently | Canonical borrowed execution contract. |
| `run_csr_group_bayesc()` | retain permanently | Sole canonical numerical core. |
| Typed result | retain permanently | Sole native handoff to aggregation/conversion. |
| `stblr_csr_group_bayesc_result_to_raw()` | retain permanently | Sole binding converter. |
| Wrapper chain aggregation | retain permanently | Sole aggregation path; schema/order sensitive. |
| Three Phase 9A fixture pairs | retain permanently | Canonical exact regression fixtures. |
| Phase-numbered historical tests/reports | defer with explicit justification | They document and protect completed boundaries; ordinary architecture no longer depends on their wording. |
| “migrated/ready” current-status wording | rename for canonical use | Current plan and capability matrix now say canonical. |
| Phase 9D3 benchmark “post-migration” wording | rename for canonical use | Phase 9E script is the canonical baseline. |
| Duplicate structural assertions | consolidate | Phase 9E provides one permanent canonical invariant set while historical tests remain regression records. |
| Context/result fields and aggregation aliases | retain permanently | All are consumed by execution, diagnostics, conversion or chain aggregation. |

No dead native field, duplicate converter, duplicate aggregation path,
temporary route alias, source-MD5 assumption for the group source, fallback or
fixture-generation dependency was found.

## 4. Files changed

- `tests/testthat/test-blr-framework-phase9e.R`: permanent canonical
  architecture, exact-reference, reproducibility, policy and protected-file tests.
- `tools/benchmarks/blr_phase9e_csr_group_bayesc.R`: canonical timing and
  completed-fit RSS baseline.
- Implementation plan and capability matrix: canonical status.
- This report.

No native production file was modified.

## 5. Canonical execution path

Public R validation and alignment feed native group/CSR preparation, the typed
context, canonical `run_csr_group_bayesc()`, typed result, sole wrapper chain
aggregation, sole binding converter, unchanged `stblr_raw_v1`, and unchanged
formatted fit.

## 6. Cleanup and naming

Stable architecture names were retained unchanged. Current documentation and
benchmark wording were renamed from migration/ready language to canonical
language. Historical phase reports and test filenames remain as immutable
provenance. No native alias was removed because none was proven dead.

## 7. Numerical core and implementation header

`blr_csr_group_bayesc_core_impl.h` retains its conventional guard, explicit
translation-unit guard and implementation-detail comment. It is included only
by `st_cpg_omp_csr_group.cpp`, defines the sole inline core and marker loop,
does not redefine compiler/Armadillo configuration, contains no Rcpp/Python
binding type, and adds no marker-loop allocation. The numerical body is unchanged.

## 8. Native adapter

The adapter remains responsible for R decoding, alignment, mapping/order and
policy preparation, CSR/native state, typed-context construction, core call,
wrapper aggregation, converter call and exception translation. Numerical MCMC,
RNG, normalization, residual/variance and posterior formulas remain in the core.

## 9. Result converter

`stblr_csr_group_bayesc_result_to_raw()` is the sole group converter. It
preserves field order/names, storage types, dimensions, dimnames, classes,
actual `NULL`, all marker/trait/group/chain orders, schema metadata, group/global
outputs, diagnostics, timing, LD-swap, selection and optional chain payloads.

## 10. Wrapper-level multichain aggregation

The core returns native per-trait/chain numerical state. The one wrapper loop
preserves chain order and forms marker and group means/standard deviations,
traces and optional chain payloads once. Conversion occurs only after the
selected aggregate branch is complete; no chain is converted twice.
`keep_chains` remains unchanged and normalization is not repeated. This stays
outside the core because it is chain-schema and ordering sensitive.

## 11. Group policy

The mapping remains borrowed, immutable and zero-based with explicit unchanged
group order/count. Probabilities, multipliers, priors, update flags/timing, both
normalization modes, group counts and global `pi`/`B` interactions are exact.
Invalid/missing mapping, empty group, dimension, equal-sample-size/shared-`ww`,
normalization-control and unsupported-trait behavior is unchanged.

## 12. Ownership

CSR storage, statistics, diagonal, mapping, group order/size, priors, fixed
initial values and LD-swap inputs are borrowed immutable resources owned by the
adapter and outlive the core call. Effects, inclusion states, residuals,
variances, mutable group/global state, RNG engines/distributions, accumulators,
diagnostics and workspaces remain chain/trait owned within execution.

## 13. Logging

Existing `std::cout` diagnostics remain outside the marker loop. They consume no
RNG, do not alter control flow or returned objects, preserve diagnostic content,
and do not introduce worker marker-loop writes. This remains the safest
binding-neutral choice.

## 14. Permanent regression fixtures

The three Phase 9A configurations remain permanent: a one-chain normalized
case, a multichain unnormalized case, and an explicit-chain-seed case. Together
they cover nontrivial mapping/order, fixed/updated group policy, both
normalization modes, retained/dropped chains, one/two cores, explicit seeds,
disabled LD-swap, selection and stable schemas. Fixture generation is not part
of ordinary tests.

## 15. Exact reference results

Group raw references: 3/3 exact. Group formatted references: 3/3 exact. No
fixture or expected numerical value was regenerated or changed.

## 16. Reproducibility

Exact behavior passes for repeated calls, one/two cores, reversed 1/2/2/1
ordering, intervening learned-annotation fits, explicit chain seeds,
retained/dropped chains, fixed/updated probabilities and multipliers, both
normalization modes and disabled LD-swap.

## 17. Public API and schema

Arguments, defaults, native signatures, routing, `NAMESPACE`, generated
wrappers, `stblr_raw_v1`, formatted fit and actual-`NULL` conventions remain
unchanged.

## 18. Protected backends

Starting-commit hashes and the full framework suite protect canonical BayesC,
BayesR, SBayesRC and fixed-prior BayesC; learned-annotation contracts, source,
fixtures, proposals/bounds/frequency and diagnostics; block-eigen; BED;
scheduled; and multivariate implementations. No protected file changed.

## 19. Performance and memory baseline

Command: `Rscript tools/benchmarks/blr_phase9e_csr_group_bayesc.R`. Each row has
one warm-up and five repetitions. The primary workload has 2,000 markers, one
trait, two nonempty ordered groups in a 20/80 mapping, 100 iterations and 25
burn-in iterations.

| Configuration | Times (s) | Mean / median / min / max (s) | Completed-fit RSS (MiB) |
|---|---|---|---:|
| fixed probabilities/multipliers, 1x1 | .09,.11,.09,.10,.08 | .094 / .09 / .08 / .11 | 130.8 |
| updated multipliers, normalized, 1x1 | .10,.09,.09,.09,.09 | .092 / .09 / .09 / .10 | 131.6 |
| updated multipliers, unnormalized, 1x1 | .09,.11,.07,.09,.10 | .092 / .09 / .07 / .11 | 129.4 |
| updated probabilities/multipliers, 1x1 | .06,.10,.08,.09,.08 | .082 / .08 / .06 / .10 | 129.3 |
| two chains/one core | .19,.17,.20,.19,.20 | .190 / .19 / .17 / .20 | 129.8 |
| retained two chains/two cores | .20,.19,.19,.20,.19 | .194 / .19 / .19 / .20 | 130.2 |

The tiny row is timer-resolution dominated and excluded from regression claims.
Phase 9E medians differ from Phase 9D3 by at most 0.01 s; RSS ranges overlap.
No material regression is evident. Warm-up and row order are explicit. Memory
is whole-process RSS sampled after completed fits, not interval peak memory.

## 20. Test results

- Phase 9A/9D1/9D2/9D3 baseline focused suite: 171 passed.
- New Phase 9E suite: 36 passed, 0 failed, 0 warned, 0 skipped.
- Baseline full suite: 4,182 passed.
- Final Phase 9A/9D1/9D2/9D3/9E focused suite: 207 passed, 0 failed,
  0 warned, 0 skipped.
- Final full suite, including group, canonical-model, learned-annotation,
  schema and protected-backend coverage: 4,218 passed, 0 failed, 0 warned,
  0 skipped.

## 21. Deviations and blockers

Interval peak-memory sampling was unavailable, so completed-fit process RSS is
reported accurately. The first Phase 9E test run exposed two stale historical
expected hashes; they were corrected to the clean starting-commit values. No
production defect or blocker remains.

## 22. Recommended next phase

> begin Phase 9F1 by mechanically extracting the learned-annotation BayesC execution block while preserving centered logistic inclusion probabilities, exponential variance multipliers, coefficient proposals, bounds, update frequency, acceptance diagnostics, RNG ordering, public schema, and all Phase 9A learned-annotation references.

## 23. Readiness marker

PHASE 9E COMPLETE — GROUP BAYESC CANONICALIZED AND STABILIZED
