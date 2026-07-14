# Unified BLR Framework: Phase 5B Report

### 1. Executive summary
Ordinary CSR BayesR execution was migrated behind typed, operator-aware boundaries with exact behavior preserved.

### 2. Repository baseline
Branch `master`, starting commit `20bd5e9`, clean status, R 4.4.1/Rtools44 GCC 13.2/C++17. Baseline: 3,592 passing tests.

### 3. Migration seam
The boundary is after existing validation, Armadillo preparation, and operator construction and before trait-chain execution.

### 4. Files changed
The BayesR binding, new implementation header, Phase 4/5 tests, benchmark, plan, capability matrix, and this report changed.

### 5. Templated execution context
The context borrows the concrete operator, prepared matrices/vectors, priors, component inputs, initialization, and immutable metadata; it owns small scalar controls. Owners outlive execution.

### 6. Execution core
The 589-line block was mechanically relocated. Arithmetic, traversal, branches, RNG/distribution construction, OpenMP scheduling, accumulation, and aggregation order did not change.

### 7. Typed result
The result carries task/chain marker summaries, components, probability and variance traces, VLE/VLD, final parameters, diagnostics, and aggregate matrices in their existing native orientations.

### 8. Centralized Rcpp converter
`csr_bayesr_result_to_raw` is the single binding-layer typed-result-to-`stblr_raw_v1` conversion callable. The core constructs no R list.

### 9. Ordinary CSR and block-eigen sharing
Both concrete operators instantiate one templated numerical body. Block-eigen source files and public behavior are unchanged.

### 10. Exact frozen references
Raw: 6/6 exact. Formatted: 6/6 exact.

### 11. Reproducibility
Repeated calls, one/two cores, reversed order, intervening BayesC, explicit seeds, multiple traits, and retained/dropped chains are exact.

### 12. Component behavior
Scales, null position, fixed/updated probabilities, traces, counts, assignments, and summaries are exact.

### 13. Public API and schema
Arguments, defaults, native signatures, routing, `NAMESPACE`, `stblr_raw_v1`, and formatted fits are unchanged.

### 14. BayesC protection
Protected files are unchanged; permanent Phase 1--3 protection tests pass 270/270.

### 15. Block-eigen protection
Protected sources are unchanged and focused block-eigen tests pass.

### 16. Performance and memory
The compact post-migration baseline used four markers, one trait, four components, 14 iterations, and 1/1, 2/1, and 2/2 chains/cores. Three-run elapsed ranges were 0--0.58, 0--0.02, and 0.02--0.02 seconds. These clock-resolution-limited runs support no speed claim. No comparable Phase 5A RSS series exists; exact references and operation preservation are primary, and whole-process RSS remains to be strengthened in Phase 6.

### 17. Tests
Phase 5A passes after replacing its obsolete MD5 assertion. New Phase 5B references/structure tests pass. CSR BayesR 222, component summary 35, combined BayesR/component/block-eigen 450, and BayesC protection 270 pass. Full suite: 3,620 passed, 0 failed, warned, or skipped.

### 18. Deviations and blockers
The shared CSR/block-eigen body required a templated context. The converter remains a named binding-layer closure pending possible Phase 6 extraction to a standalone helper. No numerical blocker remains.

### 19. Recommended Phase 6 boundary
> canonicalize and stabilize ordinary CSR BayesR, remove remaining migration-only aliases or scaffolding, verify one production numerical core and one result converter, and retain permanent exact regression fixtures.

### 20. Readiness marker
PHASE 5B COMPLETE — CSR BAYESR EXECUTION MIGRATED WITH BEHAVIOR PRESERVED
