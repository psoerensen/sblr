# Unified BLR Framework: Phase 4 Report

## 1. Executive summary

Phase 4 extracted only proven shared scalar execution infrastructure: exact
trait-major task identity/order, established seed resolution, lightweight
chain execution status and RNG-ownership vocabulary, and the common retained-
iteration predicate. Canonical unscheduled CSR BayesC uses these utilities
outside its marker mathematics. All exact BayesC references remain unchanged.

The production unscheduled CSR BayesR source, public route, marker loop, result
construction, and output schema remain unchanged. Phase 4 documents its future
migration seam but does not migrate or reroute it.

## 2. Repository baseline

- Branch: `master`.
- Starting/latest Phase 3 commit:
  `ee56b73` (`Canonicalize and stabilize CSR BayesC core`).
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1, Rtools44, GCC 13.2.0, C++17, and OpenMP.
  `pkgbuild::check_build_tools(debug = TRUE)` compiled and linked a test DLL.
- Initial native load: successful, exit 0, 459.7 seconds. This was a
  successful long-running compilation, not a timeout.
- Full baseline: 3,451 passed; 0 failed, warned, or skipped.
- Focused Phase 1--3/BayesC/BayesR/schema suite: 1,647 passed; 0 failed,
  warned, or skipped.
- Current CSR BayesR and component-summary files: 222 and 35 assertions,
  respectively; both passed in the baseline/full suite.
- Phase 3 moderate five-run medians were 0.61, 1.13, and 0.58 seconds for
  minimal output and 0.59, 1.09, and 0.64 seconds for ordinary output at
  1/1, 2/1, and 2/2 chains/cores. Sampled moderate peak RSS was
  145.07--147.78 MB.

Baseline commands were:

```text
Rscript -e "Rcpp::compileAttributes('.')"
Rscript -e "pkgload::load_all('.', compile = TRUE)"
Rscript -e "devtools::test('.')"
Rscript tools/benchmarks/blr_phase3_csr_bayesc.R --fixture=tiny --nrep=3
Rscript tools/benchmarks/blr_phase3_csr_bayesc.R --fixture=moderate --nrep=5 --peak=FALSE
Rscript tools/benchmarks/blr_phase3_csr_bayesc.R --fixture=moderate --nrep=1
```

## 3. BayesC/BayesR comparative inventory

| concern | BayesC implementation | BayesR implementation | identical | structurally similar | model-specific | safe to extract now | defer | reason/classification |
|---|---|---|---:|---:|---:|---:|---:|---|
| immutable CSR access | borrowed `CsrBayesCDataView` | shared `CsrOperator` | no | yes | no | no | yes | defer until BayesR typed view exists |
| marker diagonal | borrowed const Armadillo row | shared operator diagonal | no | yes | no | no | yes | same data, different boundary representation |
| trait indexing | `task / nchains` | `task / nchains` | yes | yes | no | yes | no | extract task identity |
| chain indexing | `task % nchains` | `task % nchains` | yes | yes | no | yes | no | extract task identity |
| task construction | trait-major `nt * nchains` | trait-major `nt * nchains` | yes | yes | no | yes | no | extract ordered task vector |
| seed resolution | three established seed helpers | same three helpers/branches | yes | yes | no | yes | no | extract exact resolver |
| explicit chain seeds | chain base plus trait offset | same | yes | yes | no | yes | no | resolver validates length |
| RNG engine ownership | chain-local `std::mt19937` | chain-local `std::mt19937` | yes | yes | no | contract | no | document physical chain ownership |
| RNG distributions | constructed at current operation sites | constructed at current operation sites | no | yes | yes | contract only | yes | moving construction could alter cached state |
| OpenMP execution | one static task loop | one static task loop | yes | yes | no | no | yes | keep schedule physically model-local in Phase 4 |
| timing | per-chain elapsed seconds | output placeholders/no matching timer | no | partial | no | status field only | yes | do not change BayesR logging/output |
| failure status | task flag and message | task flag and message | yes | yes | no | yes | no | lightweight shared status vocabulary |
| marker traversal | sorted BayesC order | BayesR order/component loop | no | partial | yes | no | yes | retain model-specific |
| log normalization | binary log-odds thresholds | max-shift categorical weights | no | no | yes | no | yes | reject false commonality |
| binary state draw | one Bernoulli-style uniform comparison | not applicable | no | no | yes | no | yes | BayesC-specific |
| categorical draw | not applicable | component cumulative draw | no | no | yes | no | yes | no second concrete use yet |
| effect sampling | binary active Gaussian | component-scaled Gaussian | no | partial | yes | no | yes | posterior mathematics |
| residual updates | CSR effect-delta update | CSR effect-delta update | no | yes | partly | no | yes | migrate through typed data view later |
| variance updates | active binary scale | component/scale weighted | no | partial | yes | no | yes | model-specific mathematics |
| retained timing | same burn-in/thinning predicate | same predicate | yes | yes | no | yes | no | extract exact predicate |
| posterior accumulation | Armadillo BayesC summaries | component-rich Armadillo summaries | no | yes | yes | predicate only | yes | containers/finalization differ |
| chain aggregation | same trait-major order | same order plus components | no | yes | partly | no | yes | defer until typed BayesR result |
| diagnostics | BayesC/LD-swap/selection | BayesR/updateE/component/LD-swap | no | partial | yes | identity/status only | yes | payloads model-specific |
| result dimensions/metadata | typed `CsrBayesCResult` | Rcpp-owned matrices/lists | no | yes | partly | task/status only | yes | BayesR typed result not yet present |
| Rcpp conversion | one centralized converter | embedded BayesR construction | no | no | yes | no | yes | Phase 5 migration seam |

Classifications are therefore: extract now for task identity/order, seed
resolution, execution status/RNG ownership vocabulary, and retained timing;
retain model-specific for posterior mathematics and state draws; defer typed
views/aggregation/results until BayesR migration; reject binary/categorical
normalization as false commonality.

## 4. Files changed

- `src/blr_scalar_execution.h`: new binding-neutral scalar task, seed,
  retained-timing, status, and RNG ownership contract.
- `src/blr_csr_bayesc_core_impl.h`: mechanically adopts those utilities
  outside the marker-update implementation.
- `src/blr_phase1_rcpp.cpp`: adds three internal test-only bridges for direct
  task, seed, and retained-iteration unit tests.
- `R/RcppExports.R` and `src/RcppExports.cpp`: generated internal registrations
  for those test bridges. No public export was added.
- `tests/testthat/test-blr-framework-phase4.R`: adds unit, source, ownership,
  BayesR-protection, exact-reference, reproducibility, and schema tests.
- `tools/benchmarks/blr_phase4_csr_bayesc.R`: reuses the Phase 3 benchmark
  fixture and measurement implementation.
- `docs/dev/blr_framework_implementation_plan.md`: records Phase 4 completion
  and the bounded Phase 5 task.
- `docs/dev/blr_model_capability_matrix.md`: records shared BayesC use and
  explicitly leaves BayesR not migrated.
- `docs/dev/blr_framework_phase4_report.md`: this report.

`src/st_cpg_omp_csr_bayesr.cpp`, `R/sparse_ld_bed_helper.R`, `NAMESPACE`, and
all protected backend sources are unchanged.

## 5. Shared infrastructure extracted

`ScalarChainTask` and `make_scalar_chain_tasks()` encode the existing
trait-major order. BayesC now iterates the prebuilt task identities; BayesR's
unchanged source has the exact matching division/modulo mapping. The vector is
allocated once before OpenMP, never in an execution or marker loop.

`resolve_scalar_chain_seed()` delegates to the established
`stblr_trait_seed()`, `stblr_chain_seed()`, and
`stblr_seed_with_chain_base()` helpers. BayesC uses it once per task. BayesR's
unchanged branches call the same helpers in the same precedence/order.

`ScalarChainExecutionStatus` contains identity, resolved seed, failure/message,
elapsed seconds, and retained samples without operator or sampler state.
BayesC replaces parallel flag/message vectors with this preallocated status
vector. BayesR has matching task failure/message and retained-sample concepts
that Phase 5 can adopt.

`scalar_iteration_is_retained()` replaces BayesC's identical predicate.
BayesR currently contains the same predicate. It is allocation-free and does
not alter retained-sample timing or summation order.

## 6. Infrastructure deliberately not extracted

- Stable log normalization: BayesC uses a binary thresholded log-odds
  calculation; BayesR uses max-shifted categorical normalization. Semantics
  and boundary handling differ, so extraction is deferred/rejected as false
  commonality.
- Categorical sampling: BayesR has a valid model-independent-looking draw,
  but BayesC has no second categorical use. It remains BayesR-specific until a
  second exact use exists.
- Binary sampling remains BayesC-specific; it is not forced through a general
  categorical helper.
- Effect, variance, residual, component, `selection_s`, and LD-swap formulas
  remain model-specific.
- Full accumulator helpers are deferred because BayesR carries component
  probabilities and uses different containers/finalization. Only the exactly
  shared retained predicate was extracted.
- Chain aggregation and typed result metadata are deferred until BayesR has a
  typed result boundary.

## 7. Canonical BayesC impact

Changes are limited to: include the shared header; construct the same task
order before OpenMP; read trait/chain from the task identity; call the exact
seed resolver; store failures/timing in the status structure; and call the
identical retained predicate.

No marker-loop line, posterior formula, arithmetic operation, marker order,
RNG draw, distribution construction site, OpenMP pragma/schedule, residual or
variance update, chain aggregation, typed result, or R conversion changed.
All frozen values remain exact.

## 8. BayesR impact

The production file's MD5 remains
`99650ce47ff4633c62c4bbc280d595ce`, exactly matching the Phase 4 starting
commit. It does not include `blr_scalar_execution.h`, invoke a Phase 4 helper,
construct a typed input/result, or change result conversion. The R route,
native signature, public API, and formatted/raw schemas are unchanged.

## 9. Task and seed infrastructure

For `T` traits and `C` chains, tasks are exactly:

```text
task = 0 .. T*C-1
trait = task / C
chain = task % C
```

Zero dimensions and task-count overflow are rejected before allocation.
Default single-chain, default multichain, and explicit-chain-base mappings
delegate to the existing helpers. Tests cover 1x1, 1x2, 2x2, explicit bases,
invalid length, repetition, and a large valid seed.

## 10. RNG ownership

Each physical chain continues to own its `std::mt19937` and all stateful
normal, uniform, chi-square, gamma, categorical/discrete, and proposal
distributions. The shared status records only the resolved seed; it owns no
engine or distribution. Distribution construction was not moved, so cached
state and RNG consumption are unchanged.

## 11. Numerical helpers

No stable-normalization or categorical-sampling helper was extracted. BayesC
binary log-odds thresholds and BayesR categorical max-shift/cumulative
semantics are intentionally isolated. This explicit deferral satisfies the
Phase 4 two-current-use rule and avoids changing floating-point/RNG boundaries.

## 12. Accumulator and result metadata

The exact retained-iteration predicate is shared. Marker accumulation,
component accumulation, final division, chain means/SD/min/max, and result
containers remain model-local so summation order and types do not change.
Task identity, resolved seed, failure, timing, and retained-count vocabulary
are now available for future typed results without changing `stblr_raw_v1`.

## 13. Ownership model

CSR values, indices, row pointers, marker diagonals, order, aligned statistics,
fixed priors/scales, and LD-friend views remain borrowed immutable storage
shared across chains. Each chain owns effects, states/components, residual,
variances, probabilities, RNG/distributions, accumulators, diagnostics, and
workspace. `ScalarChainExecutionStatus` owns no CSR payload. No per-chain CSR
copy or marker-loop allocation was introduced.

## 14. Phase 5 migration seam

Phase 5 should leave Rcpp decoding in the current entry, construct a borrowed
immutable CSR view, and add typed BayesR priors containing mixture scales and
component probabilities. Chain-local state should own effects, component
labels, residual, variance/probability state, RNG/distributions, accumulators,
diagnostics, and workspace.

The existing BayesR marker loop must move operation-for-operation behind that
boundary. It may adopt `ScalarChainTask`, `resolve_scalar_chain_seed()`,
`ScalarChainExecutionStatus`, and `scalar_iteration_is_retained()`. Its
component log weights, categorical draw, mixture scale handling, component
accumulators, diagnostics, aggregation, and typed-result conversion must
remain isolated BayesR responsibilities. The final Rcpp helper should perform
one typed-result-to-unchanged-`stblr_raw_v1` conversion.

## 15. Exact BayesC regressions

All seven normalized raw hashes and all seven normalized formatted-fit hashes
match their frozen pre-refactor values exactly. One/multiple chains,
one/two cores, explicit chain seeds, multiple traits, fixed/disabled
`selection_s`, disabled LD-swap, retained/dropped chains, repeated calls,
one/two/two/one ordering, and an intervening fit remain exact.

## 16. BayesR validation

The CSR BayesR production source hash and public routing source are unchanged.
The source contains its original task/seed calls and no Phase 4 typed route.
CSR BayesR's 222 assertions and the 35 component-summary assertions pass; the
full suite also exercises public schemas, multichain behavior, explicit seeds,
`selection_s`, and LD-swap.

## 17. Performance and memory

Commands:

```text
Rscript tools/benchmarks/blr_phase4_csr_bayesc.R --fixture=tiny --nrep=3
Rscript tools/benchmarks/blr_phase4_csr_bayesc.R --fixture=moderate --nrep=5 --peak=FALSE
Rscript tools/benchmarks/blr_phase4_csr_bayesc.R --fixture=moderate --nrep=5 --peak=FALSE
Rscript tools/benchmarks/blr_phase4_csr_bayesc.R --fixture=moderate --nrep=1
```

Phase 4 moderate medians in the first run were 0.63, 1.25, and 0.65 seconds
for minimal output and 0.66, 1.18, and 0.64 seconds for ordinary output. Means
were 1.110, 1.250, 0.686, 0.644, 1.180, and 0.642 seconds. The required repeat
gave medians 0.69, 1.17, 0.67, 0.69, 1.11, and 0.66 seconds and means 1.008,
1.174, 0.670, 0.668, 1.142, and 0.662 seconds.

Relative to Phase 3 medians (0.61, 1.13, 0.58, 0.59, 1.09, 0.64), several
short configurations vary by slightly more than 10%, while others are nearly
unchanged. Repetition changed both magnitude and direction. Inspection shows
no marker-loop change, allocation, scheduling change, RNG change, or new work
proportional to markers; task/status allocation is once per public call and
the only iteration helper is inline. The observed difference is therefore
documented as short Windows/debug-build variability, not an unexplained
algorithmic regression. No speed improvement is claimed.

Moderate peak RSS was 145.17--147.80 MB versus Phase 3's
145.07--147.78 MB, effectively identical. Tiny RSS was 141.95--144.17 MB.
The same optional `processx`/`ps` whole-child sampler polls every 10 ms; no
dependency was added. The performance and memory gates pass with the stated
measurement limitations.

## 18. Test results

- Baseline full suite: 3,451 passed; 0 failed, warned, or skipped.
- Phase 1: 102 passed; 0 failed, warned, or skipped.
- Phase 2: 63 passed; 0 failed, warned, or skipped.
- Phase 3: 105 passed; 0 failed, warned, or skipped.
- New Phase 4 file: 69 passed; 0 failed, warned, or skipped.
- Baseline focused BayesC/BayesR set: 1,647 passed; 0 failed, warned, or skipped.
- CSR BayesR: 222 passed; component summary: 35 passed.
- Combined final focused suite: 1,716 passed; 0 failed, warned, or skipped.
- Full post-change suite: 3,520 passed; 0 failed, warned, or skipped.

The final native load completed successfully in 304.6 seconds with exit status
0 and only existing compiler warnings. Testthat emitted its external
package-built-under-R startup notice; the test suite itself reported `WARN 0`.

## 19. Deviations and blockers

Three internal Rcpp functions were generated solely to execute the
binding-neutral task, seed, and retained-timing helpers independently. They
are not exported through `NAMESPACE`, do not invoke a sampler, and do not alter
public signatures.

The first sandboxed tiny RSS benchmark exited nonzero with complete
`processx` system-error-5 pipe diagnostics. The approved rerun succeeded; this
was not a timeout, compiler error, or package failure. Short moderate timings
remain variable as documented. No acceptance blocker remains.

## 20. Recommended Phase 5 boundary

Migrate the existing unscheduled CSR BayesR implementation behind typed
execution and result boundaries, adopting only the shared scalar
infrastructure validated in Phase 4 while preserving BayesR mathematics, RNG
ordering, speed, memory use, public API, and output schema.

## 21. Readiness marker

PHASE 4 COMPLETE — SHARED SCALAR INFRASTRUCTURE VALIDATED
