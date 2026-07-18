# Unified BLR Framework Phase 14C Report

## 1. Executive summary
Packed-BED BayesRC now executes through a typed binding-neutral per-chain boundary while annotation alignment, OpenMP dispatch, aggregation, and R conversion remain adapter-owned. Phase 14A numerical behavior is unchanged.

## 2. Repository baseline
Branch `master`; starting and Phase 14B commit `9d95b94`; clean initial status. Toolchain: R 4.4.1, Rtools44 GNU C++17 and OpenMP. The Phase 14B full suite passed 5,272 assertions with five opt-in skips. Phase 14A/14B tiny baselines had 0.00--0.01 second medians and approximately 123--134 MiB completed-fit RSS.

## 3. Lexical dependency inventory
| Symbol/category | Type/dimensions | Owner and mutability | Lifetime/copy | Context/result/exclusion |
|---|---|---|---|---|
| `G` | packed markers, `m × n` | fit-owned immutable | borrowed for task; no copy | genotype view |
| `maps`, `marker_order` | marker metadata, length `m` | fit-owned immutable | borrowed | context |
| `y` | `n × nt` phenotype | fit-owned immutable | borrowed | context |
| `A` | `m × P` annotation | fit-owned immutable | borrowed; no chain copy | annotation view |
| `gamma` | length `K` scales | fit-owned immutable | borrowed | component spec |
| initial alpha/sigma | `P × (K-1)`, `K-1` | fit-owned immutable | borrowed then chain-copied as mutable state | coefficient-prior spec |
| `b_init`, `B`, `E`, `ssb`, `sse` | marker/trait initial and priors | fit-owned immutable | borrowed | context |
| update/MCMC controls | scalar | adapter-resolved immutable | copied per context | context |
| logical seed | scalar unsigned | adapter-resolved immutable | copied per context | context |
| effects/residual/assignments | `m`, `n`, `m` | chain-owned mutable | one chain call | internal/final result |
| marker probabilities | `m × K` | chain-owned mutable | recomputed at established points | internal/posterior result |
| alpha/latent/workspaces | `P × (K-1)` and per-stick | chain-owned mutable | one chain call | alpha results; latent internal |
| variances and RNG | scalar/traces; engine/distributions | chain-owned mutable | one chain call | traces/finals; RNG excluded |
| task index/worker | scalar | adapter-owned | task dispatch | excluded |
| aggregation containers | marker/trace matrices | adapter-owned mutable | fit | excluded from context |
| R names/schema | binding metadata | adapter-owned | conversion | excluded |

## 4. Files changed
Added production BayesRC types and normal-probability headers, converted the core and adapter to the typed call, updated the shared annotation helper to use the probability interface, added Phase 14C tests/benchmark/report, synchronized architecture documents, and updated obsolete Phase 14A/14B structural assumptions.

## 5. Component specification
`BedBayesRCComponentSpec` borrows scales and records null zero and residual component `K-1`. Validation enforces `K>=2`, exact null scale zero, ordered positive finite active scales, and unchanged `vb*gamma[k]` semantics.

## 6. Annotation specification
`BedBayesRCAnnotationSpec` borrows the prepared marker-by-annotation matrix, records `m`, `P`, and intercept column zero. Marker matching, factor expansion, validation, and intercept preparation remain outside it.

## 7. Coefficient-prior specification
`BedBayesRCCoefficientPriorSpec` borrows initial `P × (K-1)` alpha and step variances and carries the flat-intercept policy, inverse-chi-square hyperparameters, update flag, and interval without transforming them.

## 8. Packed-genotype and annotation views
Both views are immutable references to fit-owned storage. They carry dimensions needed for native validation and introduce neither full genotype nor annotation copies.

## 9. Typed per-chain context
The context contains borrowed genotype, annotation, marker, phenotype, initialization, and prior inputs; copied MCMC/update controls; the resolved logical-chain seed; and trait/chain indices. It contains no R objects, paths, workers, aggregation state, or schema metadata.

## 10. Typed per-chain result
`BedBayesRCChainExecutionResult` owns all marker summaries, component probabilities, final effects/states/residual, variance traces/finals, alpha and step-variance means/finals, prior means, retained count, CPO, and failure information required by existing aggregation.

## 11. Callable numerical core
`run_bed_bayesrc_chain(context)` validates and binds the context, owns all mutable chain state and RNG, executes the unchanged full-sweep loop and updates, and returns the typed result. It performs no R decoding, alignment, BED I/O, task dispatch, aggregation, or conversion.

## 12. Normal probability boundary
`StandardNormalProbability` exposes stateless `cdf(double)` and `quantile(double)` C++ methods. Its dedicated implementation delegates to the exact prior `R::pnorm/qnorm` lower-tail, non-log calls. Existing caller-side finite handling and `1e-12` clipping remain in their original order; the core has no direct R namespace call.

## 13. Probit-stick preservation
The established matrix product, CDF evaluation, `remaining*(1-p)`, remaining-stick multiplication, residual final stick, probability floor, and row renormalization remain statement-ordered and component-ordered.

## 14. Latent-variable preservation
Step indicators, positive/negative branch order, CDF bounds, uniform draw, clipping, inverse CDF, marker order, and stick order are unchanged.

## 15. Annotation coefficient preservation
Alpha remains `P × (K-1)`, intercept first, with the same flat prior, sequential conditional residual updates, means/variances, normal draws, update frequency, accumulation, and final storage.

## 16. Variance and marker updates
Component log weights, inverse-CDF assignment, conditional effect draw, residual update, marker/residual/genetic variance, step inverse-chi-square draws, flags, and CPO ordering are unchanged.

## 17. Full-sweep traversal
Every iteration visits every marker once in existing order. No scheduler, skip state, candidate/active list, due bucket, or jitter was introduced.

## 18. RNG ownership
Each logical chain constructs one `std::mt19937`, one uniform and one standard-normal distribution; parameterized alpha and chi-square distributions remain local at their existing sites. The adapter resolves `seed+1000003*(trait+1)+9176*(chain+1)`; workers own no stochastic state.

## 19. Genotype, annotation, and I/O ownership
Packed genotype and prepared annotation storage are decoded/prepared once per fit and borrowed immutably. Effects, residuals, assignments, alpha, latent values, and small workspaces are chain-private. The core performs no disk access or alignment.

## 20. Native adapter
The adapter retains R validation, annotation preparation, BED decoding, marker/task preparation, static OpenMP dispatch, seed resolution, context construction, inline aggregation, optional genotype diagnostics, and inline raw conversion.

## 21. Temporary aliases
No old result alias remains. Direct aggregate matrices, Rcpp output containers, and inline conversion fragments remain intentionally adapter-local for Phase 14D.

## 22. Exact references
Final validation requires and reports raw 3/3 and formatted 3/3 exact without fixture regeneration.

## 23. Reproducibility
Final validation covers repeated calls, intervening fixtures/backends, 1/2/2/1 core order, chain-count changes, fresh/reused processes, and worker independence with only timing/core metadata normalized.

## 24. Annotation and probability identities
Finite bounded sticks, nonnegative normalized component rows, residual-stick placement, `dm=1-P(null)`, valid states, coefficient orientation, intercept order, alignment, fixed alpha, zero-alpha order, factor expansion, and validation policies are protected.

## 25. Reductions and nonreductions
The intercept-only fixed-alpha reduction to matched fixed-pi packed-BED BayesR remains exact. General BayesR and CSR SBayesRC remain documented policy/execution nonreductions.

## 26. Public API and schema
Public arguments, defaults, routing, native signature, raw/formatted schemas, actual `NULL`, wrappers, and `NAMESPACE` are unchanged.

## 27. Protected backends
Canonical packed-BED BayesR/BayesC, experimental BayesC, CSR SBayesRC/BayesR and other CSR, block-eigen, multivariate, and generated interfaces are protected by hashes and Git comparison.

## 28. Performance, memory, and I/O
The Phase 14C benchmark repeated four Phase 14A tiny workloads after warm-up for five repetitions. Medians were 0.02, 0.02, 0.01, and 0.00 seconds; ranges were 0.38, 0.02, 0.02, and 0.03 seconds. Completed-fit RSS was 123.2266, 123.3477, 123.3555, and 122.0000 MiB. Completed-fit RSS is not sampled peak RSS; Windows timer resolution and page-cache effects dominate these tiny runs. No material runtime, memory, or I/O regression is evident.

## 29. Tests
Ordinary focused Phase 14A--14C validation passed 159 assertions with two opt-in fresh-process skips. With fresh-process validation enabled it passed 163/163. The complete suite passed 5,328 assertions with zero failures, zero errors, zero test warnings, and six intentional opt-in skips. The R 4.4.3-built `testthat` warning under R 4.4.1 is a non-failing build-version warning.

## 30. Deviations and blockers
The probability interface is binding-neutral to the core while its dedicated implementation intentionally remains Rmath-backed for exact behavior. No blocker is currently known.

## 31. Recommended Phase 14D
Introduce one typed aggregate BayesRC result, centralize cross-chain aggregation and final marker-prior recomputation around typed chain results, add one named R converter, remove temporary adapter aliases, and preserve all Phase 14A references before canonicalization.

## 32. Readiness marker
PHASE 14C COMPLETE — PACKED-BED BAYESRC TYPED CHAIN BOUNDARY ACTIVE
