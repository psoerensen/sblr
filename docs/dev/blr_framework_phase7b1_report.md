# Unified BLR Framework: Phase 7B1 Report

## 1. Executive summary

The shared SBayesRC execution block was mechanically extracted into one guarded implementation header without changing numerical or public behavior.

## 2. Repository baseline

`master` at Phase 7A commit `0783581`, initially clean. R 4.4.1, Rtools44/GCC 13.2, C++17/OpenMP. The focused Phase 7A/SBayesRC/annotation/block-eigen baseline passed 768 assertions with zero failures, warnings, or skips.

## 3. Extraction seam

The relocated block starts immediately after operator construction, LD-swap friends, marker ranking, and prepared annotation/component state. It begins with `ntasks` and task/output allocation and ends after trait-chain aggregation and chain-summary standard deviations. Existing Rcpp conversion begins with `return_chain_summaries`, `n_trace`, and binding conversion lambdas.

## 4. Files changed

`src/st_sbayesrc_omp_csr.cpp` replaces the execution block with one include. `src/blr_csr_sbayesrc_core_impl.h` contains the 703 mechanically moved lines. Phase 7A source protection was updated structurally; Phase 7B1 adds permanent structural/reference tests and this report.

## 5. Extracted execution content

The block retains task allocation, OpenMP trait-chain execution, RNG and chain state, alpha initialization/update, marker probabilities, the MCMC and marker loops, component/effect/residual/variance updates, posterior accumulation, LD-swap and selection diagnostics, failures/timing, and chain aggregation.

## 6. Ordered-probit preservation

Annotation order, annotations × steps alpha orientation, stick order, `Phi(A alpha_j)`, successive remaining mass, probability floor and normalization are unchanged. Fixed/learned alpha branches and `alpha_update_every` timing are byte-for-byte relocated statements.

## 7. Operator-template sharing

The include remains inside `stblr_cpg_omp_csr_sbayesrc_impl<MakeOperator>`. Existing `BayesROperatorContext<CsrOperator>` and `BayesROperatorContext<BlockEigenOperator>` factories therefore instantiate the same single block without type erasure or duplicated loops.

## 8. R conversion

Existing inline Rcpp result construction remains in `src/st_sbayesrc_omp_csr.cpp` immediately after the include and was not redesigned.

## 9. Frozen references

Raw references: 6/6 exact. Formatted references: 6/6 exact.

## 10. Reproducibility

Repeated calls, one/two cores (normalizing only declared `input$ncores`), reversed 1-2-2-1 ordering, an intervening canonical BayesR fit, explicit seeds, multiple traits, retained/dropped chains, fixed/learned alpha, and explicit update frequencies remain exact.

## 11. Protected backends

Canonical BayesC/BayesR, block-eigen source files, BED BayesRC, fixed-prior/group/learned-annotation backends, generated wrappers, public signatures, and `NAMESPACE` are unchanged.

## 12. Tests

Baseline focused suite: 768 passed. The combined Phase 7A/Phase 7B1 suite passed 104 assertions. The final full suite passed 3,752 assertions with zero failures, warnings, or skips.

## 13. Deviations and blockers

The implementation header intentionally depends on surrounding lexical variables during Phase 7B1. No callable typed execution function or typed result population was introduced. No blocker remains.

## 14. Recommended Phase 7B2 task

> wrap the mechanically extracted shared block in an explicit operator-templated SBayesRC execution context and callable function using the Phase 7A typed contracts, while retaining the existing R result conversion until exact references pass again.

## 15. Readiness marker

PHASE 7B1 COMPLETE — SBAYESRC EXECUTION BLOCK EXTRACTED WITH BEHAVIOR PRESERVED
