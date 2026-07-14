# Unified BLR Framework: Phase 9A Report

## 1. Executive summary
Comparative binding-neutral contracts and deterministic references were established for fixed-prior, group, and learned-annotation CSR BayesC without changing production execution.

## 2. Repository baseline
Branch `master`, commit `4e7d6b4`, initially clean. R 4.4.1, Rtools44/GCC 13.2, C++17/OpenMP. Baseline: 3,917 passed, zero failures/warnings/skips.

## 3. Comparative implementation inventory
| Concern | Fixed prior | Group | Learned annotation | Classification / implication |
|---|---|---|---|---|
| CSR/order/sample size | ordinary CSR | same | same | shared contract candidate |
| Probability | fixed marker×trait vector | mutable group×trait Beta state | `logit(pi_global)+center(A eta_pi)`, bounded | mathematically different policies |
| Scale | fixed marker multiplier | mutable group multiplier | `exp(center(A eta_vb))`, bounded | mathematically different policies |
| Mapping/layout | marker order | zero-based marker→group | column-major marker×annotation | shared vocabulary only |
| Normalization | none | optional marker-weighted normalization | centered predictors then bounds | false commonality |
| Updates | immutable | group probability/variance updates | coordinate RW Metropolis | backend-specific kernels |
| Results | marker policy metadata | group means/counts/order | eta means/diagnostics | shared base plus policy payload |
| Tasks/seeds/RNG | trait-chain structure | similar | similar | defer equivalence proof to migration |
| Conversion | backend-specific | backend-specific | backend-specific | centralize separately later |

## 4. Shared-contract conclusion
Shared concepts are borrowed CSR ownership/lifetime, marker/trait ordering, scalar MCMC, seed, output and LD-swap controls, explicit policy tags, and marker/variance/timing/failure result vocabulary.

## 5. False commonality and backend-specific behavior
Marker vectors, group state, and annotation coefficients are not interchangeable. Groups are not one-hot annotations; group normalization is not annotation centering. Learned annotation links are centered logistic/exponential, not SBayesRC ordered probit. Update kernels and diagnostics remain separate.

## 6. Files changed
Added one cohesive typed header, three internal validation bridges and generated wrappers, nine fixtures and manual generator, tests, benchmark, framework/matrix updates, and this report. Production sources were untouched.

## 7. Fixed-prior contract
Borrowed marker×trait probability and multiplier arrays preserve marker order. Probabilities lie in `(0,1)` and replace global pi in marker draws; positive multipliers scale global marker variance. Inputs are fixed shared immutable data.

## 8. Group contract
The boundary uses zero-based marker→group indices, explicit count/order, group×trait initial probabilities/multipliers, Beta priors, update flags, and optional multiplier normalization. Every group must be represented; future mutable group state is chain-owned.

## 9. Learned annotation contract
Annotations are column-major marker×annotation with explicit intercept metadata. Coefficients are annotation×trait. Exact link tags are `centered_logit_offset` and `centered_exponential`; priors, RW proposal scales, update frequency, bounds, flags, and diagnostic vocabulary are explicit.

## 10. Typed result vocabularies
Common results cover marker mean/PIP, variance traces, VLE/VLD, timing, failures, diagnostics, chains, and metadata. Policy results add effective fixed inputs; group probabilities/multipliers/counts/order; or learned eta/effective policies/acceptance. No absent trace was invented.

## 11. Validation
Standard exceptions validate lifetime/ownership, dimensions, MCMC/seeds/LD-swap, fixed ranges, group mapping/nonempty groups/order, annotation layout/finiteness/coefficient dimensions, priors, proposals, frequency, bounds, and exact links.

## 12. Rcpp validation bridges
Three internal bridges build typed views over bridge-owned buffers, validate, and round-trip metadata with `invokes_sampler=FALSE`. They read no CSR and perform no sampler or policy update. `NAMESPACE` is unchanged.

## 13. Frozen fixtures
Nine RDS fixtures (three/backend) cover one chain, retained two-chain/two-core, explicit seeds/dropped chains, fixed marker inputs, group mapping/order/update/normalization, and fixed/learned annotation coefficients, bounds, frequency, and diagnostics.

## 14. Exact reference results
Fixed prior: raw 3/3, formatted 3/3 exact. Group: raw 3/3, formatted 3/3 exact. Learned annotation: raw 3/3, formatted 3/3 exact.

## 15. Reproducibility
Repeated seeded calls are exact for all policies; one/two-core and explicit-seed configurations are frozen exactly. An attempted multiple-trait fixed-prior fixture was rejected by current production dimensions and was not forced.

## 16. Policy behavior
Fixed probabilities/multipliers; group mapping/order/probability/multiplier/normalization; learned eta behavior, centered links, proposals, bounds and diagnostics; chain retention; and disabled LD-swap match exactly.

## 17. Ownership
Shared immutable: CSR/operator data, order, marker vectors, group mapping, annotations, priors, proposals. Future chain-owned: effects, states, residuals, variances, mutable group/eta state, RNG/distributions, accumulators, diagnostics, workspace. Large policy payloads are never per-chain copies.

## 18. Approved migration seams
Fixed: decode/alignment/prior/operator preparation → typed boundary → tasks/marker loop/accumulation/aggregation/conversion. Group adds mapping/group-state preparation before and group updates after the seam. Learned adds annotation/coefficient preparation before and coefficient proposals/policy updates after it. All remain ordinary-CSR, Rcpp-conversion and prepared-Armadillo coupled; none is operator-templated.

## 19. Pre-migration performance and memory
For 2,000 markers and 100 iterations, median seconds at 1×1 / 2×1 / 2×2-retained were fixed `0.11/0.68/0.68`, group `0.41/0.79/0.75`, learned `0.51/0.98/0.93`. Completed-fit RSS ranges were 129.2–138.0, 129.3–137.4, and 129.0–133.7 MiB. RSS is not peak; these are not cross-backend comparisons.

## 20. Protected backends and public behavior
Target production files and canonical BayesC/BayesR/SBayesRC, block-eigen, BED, scheduled, and multivariate sources are unchanged. Routes, schemas, and `NAMESPACE` are unchanged; wrappers changed only for internal bridges.

## 21. Tests
New Phase 9A after final reproducibility coverage: 67 passed. The focused protection suite passed 2,129 before the final 12 Phase 9A assertions; the final full suite passed 3,984. All had zero failures, warnings, or skips. Existing policy, annotation, LD-swap, canonical, schema, and routing suites are included.

## 22. Recommended migration sequence
Phase 9B fixed-prior migration; Phase 9C fixed-prior canonicalization; Phase 9D group migration; Phase 9E learned-annotation migration. Immutable policy inputs precede mutable group and Metropolis policies.

## 23. Deviations and blockers
Multiple-trait fixed-prior construction is documented unsupported. Learned acceptance counters are logged but current raw conversion exposes `NULL`; the vocabulary does not invent them. No blocker remains for fixed-prior migration.

## 24. Readiness marker
PHASE 9A COMPLETE — ANNOTATION BAYESC CONTRACTS AND REFERENCES READY
