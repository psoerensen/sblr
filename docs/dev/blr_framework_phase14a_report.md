# Unified BLR Framework Phase 14A Report

## 1. Executive summary
Public packed-BED BayesRC was audited without migrating or changing production execution. Its route, full-sweep sampler, marker-specific probit-stick prior, annotation alignment and coefficient updates, ownership, aggregation, schema, deterministic references and future extraction seam are now protected.

## 2. Repository baseline
Branch `master`; starting/Phase 13E commit `c446d1c`; clean initial status. R 4.4.1 with Rtools44 GNU C++17/OpenMP. Compile/load and the 5,169-test Phase 13E suite passed with four opt-in skips. Fast and extended CI workflows were already active.

## 3. Route and call graph
| route element | implementation | classification |
|---|---|---|
| public call | `stblr_bed(method="bayesrc")` | public production |
| method dispatch/preparation | `R/sparse_ld_bed_helper.R` | public production |
| annotation alignment/prior initialization | `R/stblr-bed-bayesrc-internal.R` | shared internal policy |
| native wrapper | `.stblr_bed_bayesrc_native()` | internal binding helper |
| native export/source | `stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc()` / matching C++ file | sole active production sampler |
| component/annotation helper | `st_bayesrc_annotation_prior.h` | shared with CSR SBayesRC conceptually and operationally |
| per-chain helper | `run_one_bayesrc_chain()` | future migration target |
| dispatch/aggregation/conversion | inline after BED decoding in the same C++ source | production adapter |
| formatting | `.as_stblr_fit()` | canonical common formatter |

No experimental second sampler, fallback, selector, adaptive route or historical active BayesRC route was found. The filename says scheduled chains, but execution is unscheduled full-sweep.

## 4. Statistical model
There are `K>=2` ordered components. Component zero is the point-mass null (`gamma[0]=0`); component `k>0` has `b_j|z_j=k ~ N(0,vb*gamma[k])`. Marker-specific probabilities arise from annotation coefficients. Marker effects use the packed-BED BayesR likelihood update; `vb` and `ve` use the existing scaled inverse-chi-square updates, `vg=Var(y-r)`, `vle` uses the packed-marker map, and `vld=vg-vle`. `dm_j=P(z_j>0)` and component probabilities/counts are retained. No MAF-dependent component scaling, selection-s, LD swap, or annotation-dependent effect variance is active.

## 5. Component ordering and scale policy
| index | name | status | multiplier | probability/output |
|---:|---|---|---:|---|
| 0 | `gamma_0.00` | null | 0 | first stick failure; first probability column |
| 1..K-2 | formatted gamma value | Gaussian | input `gamma[k]` | ordered intermediate stick failure columns |
| K-1 | formatted gamma value | Gaussian | input `gamma[K-1]` | residual final stick |

Input order, zero-based internal states, probability columns and formatted labels are preserved. Final probabilities average final marker priors; posterior means average retained marker priors. `dm=1-P(z=0)`.

## 6. Annotation input and alignment
The public input is a numeric/logical matrix or data frame with markers in rows and annotations in columns. IDs come from one recognized ID column or non-default row names; selected BIM/Glist IDs are matched and reordered exactly. Duplicate/missing IDs, nonfinite values, duplicate/empty column names and non-intercept zero-variance columns fail; extra rows are permitted and counted. Factors expand to indicator columns. An all-ones intercept is moved first or added; multiple intercepts fail. Continuous columns may be standardized and binary centering is separately controlled. The resulting dense `m x P` matrix is fit-local and passed read-only to each chain.

## 7. Probit-stick probability construction
For marker `i`, step `j`, `p_ij=Phi(A_i alpha_j)` clipped to `[1e-12,1-1e-12]`. With remaining mass `R_i0=1`, component `k<K-1` receives `pi_ik=R_ik(1-p_ik)` and `R_i,k+1=R_ik p_ik`; component `K-1` receives the residual `R_i,K-1`. Each value is floored at `pi_floor` and the row is renormalized. Thus the first stick failure is null and the last component is the residual stick.

## 8. Annotation coefficients and latent variables
`alpha` is `P x (K-1)` and `sigmaSqAlpha` has length `K-1`. Step indicator `I(z_i>j)` defines eligible markers (all at step zero, only survivors later). Latent normal values are sampled by inverse-CDF truncated normals. Coefficients are updated sequentially from Gaussian conditionals; the first annotation has a flat prior when `intercept_flat`. Other coefficients use step-specific variance. Each variance is `(sum alpha^2+b)/ChiSq(ncoef+a)`, lower bounded at `1e-12`. Updates occur every `annot_alpha_update_every` iterations when enabled; alpha/sigma retained means and final values are returned, while latent draws remain internal.

## 9. Genotype representation
The adapter validates SNP-major BED data and decodes selected rows/markers once with the blocked reader into fit-local packed marker-major storage. Marker maps preserve selected order, missing-value/scaling rules and allele frequencies. No MCMC-time file read occurs.

## 10. Ownership and lifetime
| object | owner/lifetime | mutability/copies |
|---|---|---|
| packed genotype/maps/order | fit adapter | immutable shared, one decode |
| phenotype and annotation | fit adapter | immutable shared; annotation no per-chain copy |
| effects/residual/assignments | chain helper | mutable chain-private |
| alpha/sigma/latent/workspaces | chain helper | mutable chain-private |
| RNG/distributions | chain helper/update site | chain/iteration bounded |
| task results | adapter task vector | one owned result per trait-chain |
| R objects/schema | adapter after OpenMP | binding-owned |

## 11. Marker traversal
Every iteration visits every marker once in `marker_order`, then optionally rebuilds residuals, updates `vb`, `ve`, and annotation coefficients, recalculates marker probabilities after alpha updates, and accumulates retained output. There is no adaptive scheduling, null-skip state, candidate state, due bucket or skipped marker.

## 12. RNG and distribution ownership
| stochastic object | owner/duration | risk |
|---|---|---|
| `mt19937` | one logical trait-chain | logical-chain-owned and safe |
| uniform + standard normal | one chain | safe, fit-bounded |
| marker inverse-CDF categorical | chain uniform, ordered cumulative loop | safe |
| truncated-normal uniform | local call | variable-parameter local and safe |
| coefficient normal | local coefficient update | variable-parameter local and safe |
| chi-square for alpha variance/vb/ve | local parameter update | variable-parameter local and safe |

No `static`, `thread_local`, or worker-index stochastic ownership exists.

## 13. Task dispatch and seed mapping
There are `nt*nchains` jobs, `trait=job%nt`, `chain=job/nt`, stored at the same job index under OpenMP `schedule(static)`. Seed is `seed + 1000003*(trait+1) + 9176*(chain+1)`, independent of worker assignment. Exceptions are captured in the chain result and translated after the parallel region.

## 14. Same-process reproducibility
Repeated A, A-B-A, normalized 1-2-2-1 cores, one/multiple/one chains, and intervening BED BayesR/BayesC calls are exact. Existing CSR SBayesRC coverage remains protected; the two models do not share fit-persistent RNG state.

## 15. Fresh-process reproducibility
The opt-in child-process configuration compares the nontrivial three-column, updated-alpha, two-chain/two-core fixture against its frozen raw and formatted reference with timing/core metadata normalization only.

## 16. Frozen references
Three permanent configurations are stored under `blr_phase14a_bed_bayesrc`: 1x1 intercept/fixed-alpha, 2x1 multiple annotations/updated alpha, and 2x2 multiple annotations/updated alpha. All use four components and CPO/NULL/schema coverage. The Phase 11A reference remains unchanged.

## 17. Probability and annotation identities
Alignment/reordering, unused rows, intercept addition, duplicate/missing/nonfinite/constant-column failures, finite nonnegative row-normalized stick probabilities, residual stick, `dm=1-P(null)`, valid states, `P x (K-1)` coefficients and final-prior column means are protected.

## 18. Reductions and nonreductions
Intercept-only fixed-alpha BayesRC is exactly the existing fixed-pi BayesR reduction because probit-stick intercepts reproduce the supplied probabilities. Zero coefficients imply `p=0.5`, yielding ordered probabilities `0.5,0.25,...`, not uniform. A constant non-intercept column is rejected. Fixed alpha leaves coefficients/probabilities fixed. Missing annotation is unsupported. CSR SBayesRC is conceptually/model related but execution-different (LD likelihood, selection/MAF/LD-swap policies), so general numerical equality is not asserted.

## 19. Aggregation
Chain results own marker/trace/component/alpha/sigma/final/CPO values. Inline native aggregation averages all chain summaries and final states, recomputes final marker priors from final alpha, derives component counts, and optionally computes `wy/r` from the aggregate effect. R list shaping and retained-chain construction remain entangled afterward. The future singular boundary begins immediately after task failure checks and consumes the typed chain-result vector.

## 20. Raw and formatted schemas
| category | principal fields/orientation | producer/future vocabulary |
|---|---|---|
| marker | `bm,dm,wy,r,b,state`: `m x nt` | aggregate result |
| trace | `vbs,vgs,ves,vle,vld,pis`: trace x trait | chain/native aggregate |
| variance | diagonal `covb,covg,cove,vb,vg,ve` | aggregate/converter |
| pi/component | `nt x K`, trait lists of `m x K`, names/scales/counts | aggregate plus binding labels |
| annotation | trait lists `P x (K-1)`, sigma `(K-1) x nt`, final marker priors | aggregate plus annotation names |
| diagnostics | retained counts, CPO, full-sweep flags, NULL LD swap | aggregate/presentation |
| chains | trait/chain compact lists or actual NULL | converter |

The common formatter preserves component/annotation aliases, `alpha`, `sigmaSqAlpha`, `comp_prob`, `dm_component_mean`, `ncomp`, stable marker summaries and actual NULL fields.

## 21. Typed audit contracts
`blr_bed_bayesrc_audit_types.h` defines binding-neutral component, annotation, coefficient-prior, execution and ownership specifications plus chain/aggregate result vocabularies. It validates component/null/stick dimensions, annotation/coefficient dimensions, positive finite priors, full-sweep controls and ownership. It is audit-only and invokes no sampler.

## 22. Shared infrastructure decision
| infrastructure | decision |
|---|---|
| packed BED reader/maps/view concept | reuse with BayesRC wrapper |
| task mapping/seed formula | reuse conceptually; preserve BayesRC constants |
| timing/failure/progress vocabulary | reuse conceptually |
| component specification | semantic wrapper; marker-specific probabilities require split |
| aggregation/converter conventions | reuse architecture only |
| BayesR scheduler controls | do not reuse |
| CSR SBayesRC annotation types | do not reuse directly |

## 23. Future extraction seam
Adapter preparation ends after `G`, maps/order, phenotype, annotation and priors are built. The extraction begins at `run_one_bayesrc_chain()` and ends at its returned `ChainResultBayesRC`; task dispatch remains outside. Aggregation begins after failure checks and R conversion begins with output-list construction. The recommended first unit is per-chain execution because genotype/annotation are fit-owned, state/RNG are chain-local, and dispatch/aggregation are separable.

## 24. Production behavior statement
Packed-BED BayesRC numerical source and shared annotation-prior helper remain byte-identical to `c446d1c`; no sampler, route, wrapper, schema or generated export changed.

## 25. Performance, memory, and I/O baseline
`blr_phase14a_bed_bayesrc_audit.R` provides four tiny five-repetition configurations. Medians were 0.00--0.02 seconds (ranges 0.02--0.41 seconds) and completed-fit RSS was 123.0--134.0 MiB, explicitly not peak RSS. Moderate/large and sampled-peak runs remain opt-in. BED reads occur once per fit; repeated timings benefit from page cache. No cross-model speed comparison is made.

## 26. Protected backends
Canonical BED BayesR/BayesC, CSR SBayesRC/BayesR, other CSR, block-eigen, multivariate, generated wrappers and NAMESPACE are unchanged and hash-protected.

## 27. Tests
Phase 14A adds route/contracts, annotation alignment, stick/probability identities, three exact references, same/fresh-process reproducibility, fixed-alpha reduction and protected-source checks. Focused validation passed 52 assertions with one opt-in skip; the enabled fresh-process run passed 54/54. The full suite passed 5,221 assertions with zero failures/errors and five intentional opt-in skips. The only warning was that `testthat` was built under R 4.4.3 while validation used R 4.4.1.

## 28. Confirmed defects or risks
No correctness or RNG defect was reproduced. Migration risk: the shared annotation helper calls `R::pnorm/qnorm`, so a future binding-neutral core must introduce an exactly equivalent native math boundary without changing tails or draw order. Aggregation and R conversion are inline and should remain so until after per-chain extraction. Tiny fixtures do not establish peak memory or realistic throughput.

## 29. Recommended Phase 14B
Mechanically extract the deterministic packed-BED BayesRC per-chain numerical execution and inseparable binding-neutral state into one guarded implementation header while preserving Phase 14A references, annotation alignment, probit-stick semantics, coefficient updates, RNG ownership, task dispatch, aggregation, and public schemas.

## 30. Readiness marker
PHASE 14A COMPLETE — PACKED-BED BAYESRC CONTRACT AND MIGRATION AUDIT READY
