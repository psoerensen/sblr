# Unified BLR Framework Phase 9D2 Report

## 1. Executive summary

The lexically dependent group CSR BayesC execution block is now an explicit typed callable core. `CsrGroupBayesCExecutionContext` activates the Phase 9A group contract, `run_csr_group_bayesc()` executes the unchanged numerical body, and `CsrGroupBayesCExecutionResult` carries the existing native payload to unchanged wrapper aggregation and inline R conversion.

## 2. Repository baseline

Phase 9D2 started clean on branch `master` at commit `4414362` (`Extract group BayesC execution block`), the committed Phase 9D1 baseline. `git diff --check` passed. R 4.4.1 UCRT on x86_64 Windows used Rtools44/GCC, C++17, and OpenMP. `Rcpp::compileAttributes()` left generated wrappers unchanged, native compile/load succeeded, and the focused baseline passed 1,947 assertions with zero failures, warnings, or skips.

## 3. Dependency inventory

| Dependency | Original/prepared type and dimensions | Access/ownership/lifetime | Phase 9A/context/result treatment |
|---|---|---|---|
| marker/trait/group counts and MCMC controls | scalar `int` | copied immutable metadata | explicit count, iteration, burn-in, thinning, core, and seed fields |
| scores, diagonals, effects, residuals, inclusions | `arma::mat`/`arma::Mat<int>`, trait × marker | prepared adapter storage; borrowed for call; effect/residual/inclusion mutable | explicit pointers; mutable chain state remains invocation-owned |
| phenotype sums, variance priors, initial B/E | `arma::vec` and trait × trait matrices | borrowed immutable | explicit typed pointers |
| sample sizes and initialization vectors | native vectors by trait/marker | borrowed immutable | explicit pointers; current dimension checks retained |
| CSR operator and LD friends | `STLDCSR`, `LDLDFriendsGroup` | adapter-owned, borrowed immutable, outlives call | opaque binding-neutral storage pointers with explicit row-pointer count |
| marker order | `std::vector<int>`, marker length | borrowed immutable | explicit pointer; traversal unchanged |
| marker-to-group and sizes | `arma::Row<int>`/`arma::rowvec`, marker/group length | borrowed immutable, shared, zero-based | explicit const pointers and active `GroupBayesCPolicyView` |
| group order | small ordered string vector, group length | copied metadata | explicit Phase 9A policy order; native numeric order `0..G-1` unchanged |
| initial group probability/multiplier | native trait × group vectors | borrowed immutable | explicit pointers; validated trait-by-trait without flattening/copying |
| group probability priors | group-length Armadillo rows | borrowed immutable | explicit pointers and Phase 9A policy fields |
| group multiplier prior | scalar df/scale | copied immutable metadata | explicit fields and Phase 9A policy fields |
| update and normalization policies | booleans | copied immutable metadata | explicit probability, multiplier, B/E, and normalization fields |
| global probability/B interactions | native global probability and B matrix | borrowed immutable initial inputs | explicit pointers; global output remains marker-weighted group probability |
| LD-swap controls | probability and move count; prepared friends | copied controls/borrowed friends | explicit context fields; behavior unchanged |
| timing/failures and posterior/group/variance accumulators | native vectors and Armadillo matrices | allocated and owned inside each core invocation | returned through the existing 28-slot typed native payload |
| chain seeds, chain summaries, aggregate summaries, output controls | wrapper-native values | wrapper-owned outside single-chain core | intentionally remain in sole wrapper aggregation path |
| selection inputs/outputs | no group-core selection policy; disabled conversion payload | outside core | no field invented; existing disabled selection conversion unchanged |

Every name consumed by the former lexical block is now either an explicit context field or a local/chain-owned object. All borrowed resources must outlive `run_csr_group_bayesc()`.

## 4. Files changed

- `src/blr_csr_group_bayesc_types.h`: adds the typed context, typed result, Phase 9A policy activation, and boundary validation.
- `src/blr_csr_group_bayesc_core_impl.h`: wraps the numerical body in the callable core, aliases explicit dependencies, and removes Rcpp use.
- `src/st_cpg_omp_csr_group.cpp`: constructs the context and passes the typed payload to unchanged aggregation/conversion.
- `tests/testthat/test-blr-framework-phase9d2.R`: adds structural, binding-neutrality, exact-reference, reproducibility, policy, and restriction protections.
- Historical Phase 7A/7B2/7B3/8/9A/9B3/9C tests: remove only temporary group-source MD5 entries now superseded by permanent structural and frozen-reference protections.
- This report records the implementation and validation.

## 5. Execution context

The context explicitly carries prepared score/diagonal/state matrices, phenotype and variance inputs, sample sizes, initialization vectors, CSR and LD-friend storage, marker ordering, marker-to-group mapping, group sizes, group initial values and priors, update/normalization flags, MCMC controls, seed, threads, and LD-swap controls. Large storage is borrowed; scalar controls and small order metadata are copied. Adapter-owned resources outlive the synchronous core call.

## 6. Group-contract activation

`GroupBayesCPolicyView` is populated with the borrowed zero-based mapping, explicit group count/order, borrowed priors, update policies, normalization policy, ownership flags, and lifetime flag. Each trait's borrowed initial probability and multiplier vector is validated directly through the Phase 9A contract. Groups are neither reordered nor recoded, no one-hot representation is created, and nonempty-group and current dimension restrictions remain active.

## 7. Callable core

```cpp
CsrGroupBayesCExecutionResult run_csr_group_bayesc(
    const CsrGroupBayesCExecutionContext& context);
```

The adapter constructs the context at the established seam and invokes this sole core. The guarded implementation header is included once after the native helper and concrete CSR definitions.

## 8. Typed result

The typed result owns the existing 28-slot native payload. It contains marker posterior means/PIP/state, traces for marker/genetic/residual variance and global probability, VLE/VLD, final variance/global outputs, group probability/multiplier/inclusion summaries and group sizes, timing/sample diagnostics, LD-swap diagnostics, and the native fields consumed by outer chain summaries. No new trace, diagnostic, public field, or unsupported trait output was invented.

## 9. Numerical preservation

Marker/group lookup and order, group initialization and updates, prior interpretation, normalization formula/timing, global interactions, marker traversal, arithmetic/branch order, effect draws, residual/variance updates, RNG construction/calls, seed mapping, OpenMP static scheduling, retained-sample timing, accumulation, LD swap, diagnostics, and native result population are operation-for-operation unchanged. Only binding logging changed from `Rcpp::Rcout` to `std::cout`; it does not affect numerical state or RNG.

## 10. Wrapper-level aggregation

The outer wrapper still resolves explicit chain seeds, invokes the typed single-chain core once per chain, consumes `execution_result.raw` once, preserves chain order, accumulates means/standard deviations/minima/maxima once, and handles `keep_chains` once. It remains outside the core because it is a distinct established orchestration boundary.

## 11. Inline R conversion

`cpg_raw_*` helpers and `cpg_group_raw_v1()` remain inline in `st_cpg_omp_csr_group.cpp` and unchanged. The typed payload is adapted to the same native raw reference before aggregation and conversion. R field names/order/types, dimensions, classes, actual `NULL`, schema, order metadata, diagnostics, selection payload, and optional chains are unchanged.

## 12. Exact references

All Phase 9A group fixtures remain exact: raw 3/3 and formatted 3/3. No reference was regenerated and no comparison was weakened.

## 13. Reproducibility

Exact checks passed for repeated calls, one/two cores, reversed `1,2,2,1` core order, an intervening learned-annotation fit, explicit chain seeds, retained/dropped chains, nontrivial group mapping, supported fixed/updated group policies, normalization enabled/disabled, and disabled LD swap. Only declared timing, core-count, and fixture-path metadata was normalized.

## 14. Unsupported cases

Zero-based indices remain required, invalid mappings and empty groups remain rejected, group initial/prior dimensions must match group/trait counts, shared-LD multiple-trait inputs still require equal sample sizes and matching `ww`, and no broader trait support was added.

## 15. Protected backends

MD5 structural audits and focused/full tests confirm canonical BayesC, BayesR, SBayesRC, canonical fixed-prior BayesC, learned-annotation BayesC, block-eigen, BED, scheduled, and multivariate sources are unchanged. Generated wrappers, `NAMESPACE`, public signatures, and schemas are unchanged.

## 16. Tests

- Phase 9A: 67 passed.
- Phase 9D1: 33 passed.
- Phase 9D2: 38 passed.
- Final focused suite: 2,802 passed.
- Final full suite: 4,149 passed.
- Failures: zero.
- Test warnings: zero.
- Skips: zero.

Native compile/load succeeded. R reported only that installed `testthat` was built under R 4.4.3; existing compiler unused-function warnings were unchanged.

## 17. Deviations and blockers

The requested backend-specific types header did not exist at baseline and was added in the declared scope. The first Phase 9D2 run found one test-only reference to a helper local to another test file; it was replaced with direct assertions on the shared contract validator. Numerical/reference checks already passed in that run. Binding-neutral logging uses `std::cout` because the callable core cannot contain Rcpp. No blocker remains.

## 18. Recommended Phase 9D3 task

> centralize the typed group BayesC result-to-R converter, remove obsolete migration assertions and aliases, establish post-migration runtime and memory baselines, update framework documentation, and close the group migration.

## 19. Readiness marker

PHASE 9D2 COMPLETE — GROUP BAYESC TYPED EXECUTION BOUNDARY ACTIVE
