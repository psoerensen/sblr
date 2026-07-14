# Unified BLR Framework: Phase 5A Report

## 1. Executive summary

Phase 5A established binding-neutral typed BayesR input/result contracts and
compact frozen ordinary-CSR BayesR references without changing production
execution. The production trait-chain loop, marker loop, aggregation, public
route, and raw conversion remain in their original implementation.

## 2. Repository baseline

- Branch: `master`.
- Starting commit: `b307db3` (`Extract shared scalar BLR infrastructure`).
- Initial status: clean; `git diff --check` passed.
- Toolchain: R 4.4.1, Rtools44/GCC 13.2.0, C++17, OpenMP.
- Baseline native load: passed in 245.1 seconds.
- Baseline full suite: 3,520 passed; 0 failed, warned, or skipped.
- Baseline CSR BayesR/component-summary tests: 222/35 passed.

## 3. Existing BayesR structure

`stblr_cpg_omp_csr_bayesr()` decodes R inputs and constructs Armadillo and CSR
operator storage. The shared template then validates component inputs, creates
trait-major chain tasks, executes the existing shuffled marker loop, accumulates
posterior/component summaries, aggregates chains, and constructs `stblr_raw_v1`.
The block-eigen entry instantiates the same numerical template with its existing
operator construction and diagnostics.

## 4. Frozen references

Six compact four-marker configurations protect complete normalized raw and
formatted outputs: one chain; two chains/one core; two chains/two cores; two
traits; explicit seeds with retained chains; and explicit component scales with
fixed probabilities. Together they cover default and explicit scales,
fixed/updated component probabilities, null/non-null assignments, component
summaries, variance output, chain order, actual `NULL`, and disabled LD-swap.
Each RDS includes commit, R/compiler, fixture, dimensions, scales, seed/control,
chain/core, output, and schema metadata.

## 5. Typed input contracts

`CsrBayesRDataView` describes borrowed immutable CSR pointers, diagonal and
sample-size storage, marker/trait dimensions, shared-read-only ownership, no
per-chain payload, and the lifetime requirement. `BayesRComponentSpec`
preserves scale and probability order, zero-based null index, variance-
multiplier interpretation, Dirichlet prior, and update flag.
`CsrBayesRControls`, priors, output control, and `CsrBayesRExecutionInput`
represent current execution metadata plus marker/trait order.

## 6. Typed result contract

`CsrBayesRResult` and `CsrBayesRChainResult` define dimensioned native payload
categories for marker means/PIP/effects/components, component probabilities and
counts, final/mean proportions, variance/VLE/VLD traces, timing, failures,
diagnostics, optional chains, and input order. This is vocabulary only in
Phase 5A; the sampler does not return it yet. No BayesRC or annotation fields
were added.

## 7. Validation

Standard C++ validation rejects zero dimensions, invalid CSR/diagonal/sample-
size views, mutable or per-chain storage, insufficient lifetime, invalid
component counts/null index/scales/probability dimensions or values, invalid
Dirichlet parameters, MCMC/chain/core/seed errors, inconsistent output control,
and invalid residual-update or LD-swap controls.

## 8. Rcpp validation bridge

`blr_phase5a_validate_bayesr_contract_cpp()` converts compact internal R
metadata into typed structures, validates them, and returns the same normalized
sections. It reads no CSR file, constructs no production operator, references
no sampler, and reports `invokes_sampler = FALSE`. Its generated wrapper is
internal and `NAMESPACE` is unchanged.

## 9. Approved Phase 5B seam

The seam is immediately after existing Rcpp validation, Armadillo conversion,
and CSR/operator construction, and immediately before mixture-vector/order,
trait-chain task allocation and execution, after `operator_context.op` and
`operator_context.ld_swap_friends` are bound. Phase 5B will move execution
through the typed input and result while
preserving the current operational order.

## 10. Ordinary CSR versus block-eigen

Phase 5B can type the ordinary-CSR adapter at the marked seam while retaining a
single shared numerical implementation. The block-eigen entry can remain on
its existing construction/instantiation route. Phase 5A did not modify either
block-eigen source, its adapter, or the shared numerical body.

## 11. Test results

- New Phase 5A file: 72 passed; 0 failed, warned, or skipped.
- Frozen raw comparisons: 6/6 exact.
- Frozen formatted comparisons: 6/6 exact.
- Existing CSR BayesR: 222 passed; component summary: 35 passed.
- Phase 1--4, block-eigen, BayesC, schema, consistency, field, and public
  interface protections passed in the full suite.
- Full post-change suite: 3,592 passed; 0 failed, warned, or skipped.

## 12. Public behavior statement

Public arguments, defaults, routing, native signatures, BayesR mathematics,
component semantics, marker traversal, RNG calls, seed mapping, OpenMP
scheduling, aggregation, `stblr_raw_v1`, formatted fit, and `NAMESPACE` are
unchanged. The production BayesR translation unit is byte-identical to the
starting commit; the seam is encoded here to preserve the Phase 4 source hash.

## 13. Deviations and blockers

The first post-change compilation found a duplicate C++ result member name;
renaming the dimension field resolved the compiler diagnostic, and a fresh
build passed. No production code was implicated. There are no Phase 5A
blockers.

## 14. Recommended Phase 5B task

> move the existing ordinary-CSR BayesR trait-chain execution, marker loop, accumulation, and aggregation operation-for-operation behind the Phase 5A typed boundary, then convert the typed result through one centralized Rcpp converter while preserving the shared block-eigen instantiation.

## 15. Readiness marker

PHASE 5A COMPLETE — BAYESR CONTRACTS AND REFERENCES READY
