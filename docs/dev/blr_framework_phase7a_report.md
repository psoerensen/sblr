# Unified BLR Framework: Phase 7A Report

## 1. Executive summary

Binding-neutral ordinary-CSR SBayesRC contracts and deterministic references were established without changing production execution.

## 2. Repository baseline

`master` at `f6785f7`, initially clean. R 4.4.1 with Rtools44/GCC 13.2, C++17 and OpenMP. The full baseline was 3,648 passed, zero failures, warnings, or skips.

## 3. Existing implementation inventory

| Concern | Current source | Dimensions | Ownership/state | Phase 7A representation | Phase 7B requirement |
|---|---|---|---|---|---|
| Public route/alignment | `R/stblr-csr-annot.R`, `R/stblr-csr-sbayesrc.R` | markers × traits | R-owned/shared | IDs and typed dimensions | retain before seam |
| CSR/operator | `src/st_sbayesrc_omp_csr.cpp` | markers, nonzeros | binding-owned immutable | `CsrSBayesRCDataView` | borrow after construction |
| Annotation design | R preparation then native `A` | markers × annotations, column-major | shared immutable | `SBayesRCAnnotationDesignView` | borrow prepared matrix |
| Components | native `gamma` | K, null at index zero | shared immutable | `SBayesRCComponentSpec` | preserve order/scales |
| Probability policy | `st_bayesrc_annotation_prior.h` | markers × K | chain-derived | `SBayesRCProbabilityPolicy` | move unchanged probit stick logic |
| Alpha | native sampler | annotations × (K-1) | chain-local mutable | `SBayesRCAlphaSpec` | preserve update timing/RNG |
| MCMC/RNG | native trait-chain loop | traits × chains | chain-local | `SBayesRCControls` | move operation-for-operation |
| Accumulation/result | native sampler/converter | existing raw dimensions | task-owned then aggregated | `CsrSBayesRCResult` | populate typed result |

The public spelling is SBayesRC for ordinary CSR. BED BayesRC is a distinct protected route.

## 4. Files changed

Added typed contracts, an internal bridge, six fixture definitions/outputs, permanent tests, a benchmark driver, and this report. Generated internal Rcpp wrappers changed only for the bridge. The plan and capability matrix record Phase 7A status.

## 5. Frozen references

Six configurations cover fixed one chain, fixed two chains, learned alpha/two cores, learned alpha/explicit seeds/retained chains, multiple traits, and explicit four-component scales. All use three ordered annotations including an explicit intercept and disabled LD-swap.

## 6. Typed data and ownership contracts

CSR row pointers, indices, values, diagonal and sample sizes are borrowed read-only and must outlive execution. The column-major marker-by-annotation matrix is likewise borrowed and shared; per-chain CSR or annotation payloads are rejected.

## 7. Component contract

Components retain input order, variance-multiplier interpretation, and a null scale of exactly zero at index zero. Non-null scales must be positive finite.

## 8. Alpha and probability-policy contract

Alpha is annotations × `(K-1)` in column-major order. An explicit first intercept column is never silently inserted by the native contract; `intercept_flat` controls its prior. Each ordered step uses `Phi(A alpha_j)`. Component probabilities are successive remaining-mass products, floored by `pi_floor`, then row-normalized. Stick order is 1 through `K-1`; component zero is the null/reference boundary. Alpha variance has one value per step and inverse-chi-square-equivalent controls `sigmaSqAlpha_a/b`; updates occur every `alpha_update_every` iterations.

## 9. Typed controls

Contracts cover iterations, burn-in, thinning, chains, cores, base/explicit seeds, retention, variance and residual flags, alpha updates, residual rebuilding, output diagnostics, and LD-swap controls.

## 10. Typed result vocabulary

The result contract covers marker means/PIP, marker component probabilities and assignments, counts/proportion traces/finals, alpha traces/means/variances/diagnostics, all variance traces including VLE/VLD, timing/failure diagnostics, optional chain payloads, IDs, names, and dimensions.

## 11. Validation

Standard C++ exceptions reject invalid data ownership/lifetime, annotation layout/dimensions/names/values, component/null/scale contracts, alpha dimensions/values/prior/update frequency, probability transformation/order/floor, MCMC/seeds, output consistency, and LD-swap controls.

## 12. Rcpp validation bridge

`blr_phase7a_validate_sbayesrc_contract_cpp()` creates borrowed views over small bridge-owned buffers, validates, and round-trips metadata. It invokes no sampler, reader, alpha update, or probability execution.

## 13. Frozen-reference results

Raw: 6/6 exact. Formatted: 6/6 exact.

## 14. Reproducibility

The permanent configurations protect repeated calls, one/two cores, explicit seeds, multiple traits, and retained/dropped chains. Existing focused suites continue to protect reversed order and intervening fits.

## 15. Annotation and alpha behavior

Annotation order/intercept treatment, fixed and learned alpha, update frequencies 1/2/3/10, alpha and variance summaries, marker probabilities, counts, and component summaries match frozen production behavior exactly. The current CSR output has posterior alpha summaries rather than a retained per-iteration alpha trace/acceptance counter; this absence is preserved.

## 16. Public behavior and protected backends

Production SBayesRC, canonical BayesC/BayesR, block-eigen, BED BayesRC, prior/group/learned-annotation backends, public routing, arguments/defaults, `stblr_raw_v1`, formatted fit, and `NAMESPACE` are unchanged.

## 17. Pre-migration performance and memory baseline

Command: `Rscript tools/benchmarks/blr_phase7a_csr_sbayesrc.R`. For the representative 500-marker, one-trait, four-annotation, four-component, 40-iteration workload, three-run medians were 0.030 s (fixed alpha, 1 chain/1 core), 0.070 s (learned alpha, 2 chains/1 core), and 0.030 s (learned alpha, 2 chains/2 cores). The compact configurations sampled 136.6--136.9 MiB whole-process RSS. RSS is sampled after each run rather than being a sampler-only interval peak; short Windows timings support no speed claim.

## 18. Approved Phase 7B execution seam

The seam is after existing Rcpp decoding, annotation/marker alignment, validation, Armadillo conversion, CSR/operator construction, and prepared annotation/component state; below it are trait-chain execution, alpha/probability updates, marker loop, accumulation, aggregation, and R conversion. Dependencies map to the Phase 7A data, annotation, component, alpha, probability, prior, control, and output contracts. The native implementation is operator-templated for ordinary CSR and block-eigen and uses prepared Armadillo state; Phase 7B must preserve that compile-time sharing without migrating block-eigen independently.

## 19. Tests

Baseline full suite: 3,648 passed. The new Phase 7A suite has 80 passing assertions. The combined SBayesRC/helper, annotation interface/backend/VLE-VLD, component-summary, and block-eigen filter has 857 passing assertions. The final full suite has 3,728 passed, zero failures, warnings, or skips; it also includes BayesC/BayesR, schema, consistency, field-inventory, LD-swap, and public-interface coverage.

## 20. Deviations and blockers

Multiple traits are supported and covered. Alpha acceptance/proposal counts and retained alpha traces are not currently public CSR outputs, so their absence is recorded rather than invented. No migration blocker exists.

## 21. Recommended Phase 7B task

> mechanically move the existing ordinary-CSR SBayesRC/BayesRC trait-chain execution, alpha and marker-probability updates, marker loop, posterior accumulation, and aggregation behind the Phase 7A typed boundary while preserving its mathematics, RNG ordering, annotation behavior, public API, schema, speed, and memory use.

## 22. Readiness marker

PHASE 7A COMPLETE — SBAYESRC CONTRACTS AND REFERENCES READY
