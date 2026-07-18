# Unified BLR Framework Phase 14B Report

## 1. Executive summary
The deterministic packed-BED BayesRC per-chain numerical execution was mechanically extracted into one guarded implementation header without changing Phase 14A behavior.

## 2. Repository baseline
Branch `master`; starting and Phase 14A commit `8f0a391`; clean initial status. R 4.4.1 with Rtools44 GNU C++17/OpenMP. Phase 14A focused 52, enabled fresh-process 54, and full-suite 5,221 assertions passed; its tiny benchmark medians were 0.00--0.02 seconds and completed-fit RSS 123--134 MiB.

## 3. Production call graph
`stblr_bed("bayesrc")` -> public BED dispatch -> annotation alignment/prior initialization -> internal R native wrapper -> native export -> BED decode -> static trait-chain tasks -> guarded `run_one_bayesrc_chain()` -> inline aggregation/raw conversion -> common formatter.

## 4. Extraction target
`run_one_bayesrc_chain()` was the sole active per-logical-chain sampler. Its inputs were already prepared and fit-owned, its RNG/state were chain-local, and dispatch/aggregation were separable.

## 5. Extraction seam
The moved block began at `struct ChainResultBayesRC` (former source line 6) and ended at the closing brace of `run_one_bayesrc_chain()` (former line 179). It now follows the common BED and annotation-helper includes. The adapter resumes at `// [[Rcpp::export]]`.

## 6. Declaration classification
| declaration | classification/action |
|---|---|
| `ChainResultBayesRC` | binding-neutral inseparable result; moved |
| `sample_marker_bayesrc()` | chain-only marker update; moved |
| `run_one_bayesrc_chain()` | extraction target; moved |
| shared probit/truncated-normal/alpha helpers | shared helper; retained unchanged |
| BED reader/maps/order and annotation preparation | adapter-owned; retained |
| job vector/task mapping | adapter-owned dispatch; retained |
| aggregate matrices/Rcpp lists | aggregation/binding state; retained |

## 7. Files changed
The production source replaces the block with one include; `blr_bed_bayesrc_core_impl.h` contains the identical body plus guard/comment. Phase 14A source-location tests were updated; Phase 14B tests/report and plan/matrix status were added. No wrappers or schemas changed.

## 8. Lines mechanically moved
`MECHANICAL_LINES=174`  
`HEADER_BODY_LINES=174`  
`IDENTICAL=TRUE`

Comparison was against commit `8f0a391`; no numerical line differs.

## 9. Extracted execution content
The header owns the temporary result vocabulary, marker component/effect update, chain-local engine/distributions, effects/residual/assignments, marker probabilities, alpha/sigma state, MCMC/full sweep, variance/annotation updates, retained posterior accumulation, CPO and finalization.

## 10. Component semantics preservation
Component zero remains the point-mass null with scale zero; active component `k` uses variance `vb*gamma[k]`; assignment and output order remain zero-based; the last component remains the residual stick and `dm=P(component>0)`.

## 11. Probit-stick preservation
The unchanged shared helper computes `p=Phi(A alpha)`, ordered failure masses and the residual final stick, applies the same clipping/floor, then renormalizes each row. `R::pnorm/qnorm` remain at their original boundary.

## 12. Annotation alignment boundary
ID detection/matching/reordering, extra-row counting, duplicate/missing/nonfinite validation, factor expansion, constant-column validation, annotation names and intercept insertion/movement remain entirely in R/adapter preparation.

## 13. Latent-variable preservation
Step indicators remain `I(component_i>j)`. The inverse-CDF truncated-normal branch, uniform draw, CDF clipping, stick/marker order and update flag/frequency are unchanged in the shared helper.

## 14. Annotation coefficient preservation
The `P x (K-1)` orientation, first-column flat-prior option, sequential conditional residualization, mean/variance formulas, parameterized normal draws, alpha/sigma accumulation and final storage are unchanged.

## 15. Variance-update preservation
Step-specific scaled inverse-chi-square, base marker variance, residual variance, genetic variance, `vle/vld`, update flags/order and CPO calculations are identical.

## 16. Full-sweep traversal
Every iteration visits every marker once in `marker_order`. No BayesR scheduler control, candidate/active/due state, skipping or jitter was introduced.

## 17. RNG ownership preservation
Each logical trait-chain constructs one `mt19937`, uniform and standard-normal distribution using `seed + 1000003*(trait+1) + 9176*(chain+1)`. Parameterized coefficient normals and chi-square distributions remain local. No static/thread-local or worker-owned state exists.

## 18. Genotype and I/O preservation
Blocked BED validation/decoding and marker-map construction remain adapter-side. The core borrows fit-owned packed storage, owns no genotype bytes/path/handle, performs no decoding or disk I/O, and introduces no chain genotype copy.

## 19. Existing dispatch, aggregation, and conversion
Static OpenMP task enumeration, failure capture, chain-result vector, all cross-chain formulas, final-prior recomputation, component counts, optional `wy/r`, compact chains and raw-list construction remain inline and statement-identical.

## 20. Exact references
All three Phase 14A raw and all three formatted references match exactly; fixtures were not regenerated.

## 21. Reproducibility
Repeated calls, A-B-A, normalized 1-2-2-1 core order, fresh/reused process, intervening BED BayesR/BayesC, one/multiple/one chains and worker assignment remain exact.

## 22. Annotation and probability identities
Stick/marker probabilities remain finite/nonnegative and row-normalized; final residual stick, `dm=1-P(null)`, assignment range, `P x (K-1)` alpha, alignment/intercept/fixed-alpha/final-prior identities remain exact.

## 23. Reductions and nonreductions
The intercept-only fixed-alpha configuration remains exactly reducible to fixed-pi packed-BED BayesR. Zero alpha retains ordered half-stick probabilities; constant non-intercepts and missing annotations remain rejected; CSR SBayesRC remains conceptually related but execution/policy-different.

## 24. Protected backends
Canonical BED BayesR/BayesC, experimental/sparse BayesC, CSR SBayesRC/BayesR, other CSR, block-eigen, multivariate, generated wrappers and NAMESPACE remain unchanged.

## 25. Performance, memory, and I/O
The Phase 14A benchmark was rerun on the identical path. Configuration medians were 0.00, 0.01, 0.00, and 0.00 seconds; observed ranges were 0.39, 0.02, 0.02, and 0.01 seconds. Completed-fit RSS was 134.3984, 122.9688, 122.9766, and 123.0039 MiB. These short Windows timings are noisy and page-cache affected, and completed-fit RSS is not sampled peak RSS. No material runtime, completed-fit-memory, or I/O change is attributable to the include-only extraction.

## 26. Tests
Ordinary focused validation passed 103 assertions with one opt-in fresh-process skip: Phase 14A passed 52 and Phase 14B passed 51. With fresh-process validation enabled, Phase 14A and Phase 14B passed 105/105 assertions with no skip. The complete suite passed 5,272 assertions with zero failures, zero errors, zero test warnings, and five intentional opt-in skips. The only emitted warning was that `testthat` was built under R 4.4.3 while the runtime is R 4.4.1; it is non-failing.

## 27. Deviations and blockers
The shared probit/latent/alpha helper was deliberately not moved because it is shared and independently unchanged. Rmath binding neutrality is deferred to Phase 14C. No blocker remains.

## 28. Recommended Phase 14C
Replace the lexically dependent packed-BED BayesRC chain include with explicit typed component, annotation, prior, and per-chain execution contracts; a callable numerical core; a typed per-chain result; and a binding-neutral normal-CDF/inverse-CDF boundary while retaining annotation alignment, task dispatch, aggregation, and R conversion until all Phase 14A references pass again.

## 29. Readiness marker
PHASE 14B COMPLETE — PACKED-BED BAYESRC PER-CHAIN EXECUTION BLOCK EXTRACTED
