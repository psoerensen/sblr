# Unified BLR Framework: Model, Architecture, and Migration Capability Matrix

**Status:** Phase 4 shared scalar infrastructure recorded
**Date:** 2026-07-13  
**Target location:** `docs/dev/blr_model_capability_matrix.md`

## 1. Purpose

This matrix tracks:

- current model capabilities;
- target kernel family;
- probability and scale policies;
- likelihood operators;
- migration treatment;
- current reference implementation;
- extension readiness;
- performance and memory requirements.

Status labels:

- **Current** — active working implementation.
- **Preserve** — retain kernel and refactor boundaries.
- **Extract** — move reusable infrastructure into shared core.
- **Correct** — localized fix required.
- **Rewrite** — new canonical implementation required.
- **Planned** — accepted future capability.
- **Deferred** — intentionally later.
- **Removed** — legacy implementation deleted after stabilization.

---

## 2. Kernel families

| Kernel family | Scope | Migration strategy |
|---|---|---|
| `scalar` | ST BayesC/R/RC, annotation, hierarchy | Preserve efficient kernels and extract shared infrastructure incrementally |
| `small_mt` | explicit multivariate patterns/components | New coherent implementation |
| `factor` | factor-analytic BLR and evidence models | New implementation with shared factor utilities |

---

## 3. Policy vocabulary

### State policies

| ID | Meaning |
|---|---|
| `binary` | Null versus active |
| `mixture` | BayesR/BayesRC component |
| `trait_pattern` | Explicit multivariate activity pattern |
| `pattern_mixture` | Trait pattern × mixture component |
| `sparse_factor` | Sparse row-factor membership |

### Probability policies

| ID | Meaning |
|---|---|
| `global_binary` | Global BayesC inclusion |
| `fixed_marker` | Fixed marker inclusion probabilities |
| `global_dirichlet` | Global mixture probabilities |
| `group` | Group-specific probabilities |
| `annot_logit` | Annotation-dependent binary probability |
| `annot_probit_stick` | Annotation-dependent mixture probability |
| `factor_membership` | Annotation-informed factor membership |

### Scale policies

| ID | Meaning |
|---|---|
| `unit` | No marker-specific scale |
| `fixed_marker` | Fixed marker multiplier |
| `maf_s` | `selection_s` MAF scaling |
| `component` | BayesR/BayesRC component multiplier |
| `group` | One categorical group layer |
| `hierarchy` | Multiple categorical layers |
| `log_additive` | Overlapping annotation log-scale model |
| `composite` | Product of compatible scales |

### Covariance policies

| ID | Meaning |
|---|---|
| `scalar` | Single-trait variance |
| `diagonal` | Independent trait variances |
| `full` | Full small-\(T\) covariance |
| `factor_analytic` | Low-rank plus diagonal |

---

## 4. Regression model matrix

| Model | Kernel | State | Probability | Scale | Covariance | Operator | Current implementation | Migration treatment | Extension readiness |
|---|---|---|---|---|---|---|---|---|---|
| ST BayesC CSR | scalar | binary | global_binary | unit | scalar | CSR | Canonical typed unscheduled CSR using proven shared scalar execution utilities | Migrated; canonical; legacy path removed | Active architecture reference; public schema unchanged; performance and memory validated |
| ST BayesC scheduled CSR | scalar | binary | global_binary | unit | scalar | CSR | Canonical corrected chain-owned RNG implementation with active scheduler contracts, typed borrowed context, callable core, typed result, one named converter, and one native aggregation path | Canonical; migration scaffolding removed; public route/schema unchanged; permanent Phase 10B corrected references and runtime/completed-fit-RSS baselines active | Fit-local/thread-assignment reproducibility exact; unsupported cases preserved; scheduled BayesR remains unsupported |
| ST BayesC BED | scalar | binary | global_binary | unit | scalar | BED | Current BED kernels | Preserve decoder/kernel, refactor boundary | After CSR architecture proven |
| ST BayesC scheduled packed BED | scalar | binary | global_binary | unit | scalar | SNP-major packed BED | Corrected logical-chain RNG active; public multichain 417-line execution body mechanically extracted; existing converter retained | Phase 11B references exact; typed callable boundary not active; migration in progress | Experimental single-chain route, BayesR and BayesRC unchanged |
| ST BayesC fixed marker prior | scalar | binary | fixed_marker/global | fixed_marker | scalar | CSR | Canonical typed implementation with borrowed `pi_marker`/`vb_multiplier`, callable core, typed result, one binding converter, and one wrapper aggregation path | Canonical; public route/schema unchanged; permanent exact references and runtime/completed-fit-RSS baselines active; trait restriction preserved; migration scaffolding removed | Stable canonical fixed-prior reference |
| ST BayesC + `selection_s` | scalar | binary | global_binary | maf_s | scalar | CSR | Current | Preserve and migrate | Good composite-scale test |
| ST group BayesC | scalar | binary | group | group | scalar | CSR | Canonical typed implementation with borrowed zero-based mapping, active group policy, callable core, typed result, one binding converter, and one wrapper aggregation path | Canonical; public route/schema and normalization semantics unchanged; permanent exact references and runtime/completed-fit-RSS baselines active; unsupported cases preserved; migration scaffolding removed | Stable canonical group reference; group normalization and mutable state remain specific |
| ST learned annotation BayesC | scalar | binary | annot_logit | learned annotation scale | scalar | CSR | Canonical typed implementation with borrowed annotation design, active learned policy, callable core, typed result, one binding converter, and one wrapper aggregation path | Canonical; public route/schema unchanged; permanent exact references and canonical runtime/completed-fit-RSS baselines active; unsupported cases preserved; migration scaffolding removed | Stable canonical learned-annotation reference; centered logistic/exponential MH policy, proposals, bounds, and update frequency remain specific |
| ST BayesR CSR | scalar | mixture | global_dirichlet | component | scalar | CSR | Canonical | Typed borrowed context, operator-aware templated core, typed result, one converter; permanent exact references and runtime/memory baseline | Stable canonical scalar mixture reference |
| ST BayesR BED | scalar | mixture | global_dirichlet | component | scalar | BED | Current | Preserve and migrate after CSR | Operator reuse test |
| ST BayesRC scheduled-chains BED | scalar | mixture | annotation probit stick | component | scalar | SNP-major packed BED | Unscheduled full-sweep annotation-component implementation; route/model, decoding and RNG audit complete | Production unchanged; migration not started | Adaptive scheduled controls do not apply; chain-local RNG reference active |
| ST SBayesRC CSR | scalar | mixture | annot_probit_stick | component/maf_s | scalar | CSR | Canonical typed borrowed CSR/annotation context with operator-aware templated core, typed result, and one ordinary-CSR converter | Ordered-probit and alpha behavior preserved; shared task/seed/status infrastructure active; public route/schema unchanged; permanent exact references and runtime/completed-fit-RSS baseline established | Canonical annotation-aware mixture reference |
| ST BayesRC BED | scalar | mixture | annot_probit_stick | component | scalar | BED | Current | Preserve and migrate | Cross-operator policy test |
| ST hierarchical BayesC | scalar | binary | global_binary | hierarchy | scalar | CSR first | Not implemented | New policy on migrated scalar core | Planned |
| ST hierarchical BayesR | scalar | mixture | global_dirichlet | component × hierarchy | scalar | CSR first | Not implemented | Compose after hierarchy | Planned |
| ST hierarchical BayesRC | scalar | mixture | annot_probit_stick | component × hierarchy | scalar | CSR first | Not implemented | Compose after BayesRC migration | Planned |
| MT BayesC legacy | small_mt | trait_pattern | global pattern | unclear/inconsistent | full/diagonal | dense/CSR variants | Legacy experimental | Retain only as reference; rewrite canonical model | Not extension base |
| Canonical MT BayesC | small_mt | trait_pattern | global pattern | unit/hierarchy later | fixed/diagonal/full | dense exact first | Not implemented | New implementation | Planned |
| MT BayesR | small_mt | pattern_mixture | global_dirichlet | component | diagonal/full | dense then CSR | Not implemented | New implementation | Planned later |
| MT BayesRC | small_mt | pattern_mixture | annotation-informed | component | diagonal/full | dense then CSR | Not implemented | New implementation | Deferred |
| Factor-analytic MT-BLR | factor | sparse_factor | factor policy | factor scale | factor_analytic | dense/summary | Not implemented | New implementation | Planned later |

---

## 5. Evidence-factor matrix

| Model | Row unit | Likelihood | Membership | Annotation role | Migration status |
|---|---|---|---|---|---|
| Gene Gaussian factor model | gene | Gaussian with uncertainty | dense or sparse | none initially | Planned |
| Annotation-informed gene factor | gene | Gaussian | sparse factor membership | gene-factor probability | Planned later |
| Multi-view gene factor | gene | view-specific Gaussian | shared factors | annotations and anchors | Deferred |
| SNP evidence factor | SNP | transformed/uncertainty-aware | sparse | SNP-factor probability | Deferred until gene model validated |
| LD-aware SNP factor | SNP | correlated residual | sparse | SNP-factor probability | Research extension |
| Pathway propagation | pathway | posterior mapping | inherited | pathway annotations | Post-processing first |

---

## 6. Operator capability matrix

| Operator | Current status | Strengths | Migration treatment | First supported kernels |
|---|---|---|---|---|
| CSR | Current and validated | scalable summary-statistics residual updates | Preserve and formalize ownership/interface | scalar BayesC/R/RC |
| Block eigen | Current internal/experimental | shared operator concept | Preserve implementation; migrate after CSR | scalar models |
| Packed BED | Current and efficient | direct genotype access and low memory | Preserve decoder and score/update logic; separate binding | scalar BayesC/R/RC |
| Dense individual | Partial/reference | exact shared-design likelihood | New clean operator | canonical MT reference |
| Dense summary | Legacy/partial | exact sufficient statistics when complete | New clean operator based on explicit contract | canonical MT reference |
| Evidence matrix | Not implemented | direct factor input | New frontend | factor models |

---

## 7. Migration treatment matrix

| Component | Treatment | Reason |
|---|---|---|
| Unscheduled CSR BayesC hot loop | Migrated and canonical | typed boundary active; one preserved production implementation; legacy path removed |
| CSR operator logic | Preserve/extract | strongest current operator abstraction |
| BayesR/BayesRC scalar kernels | Preserve/migrate | efficient and validated |
| BED decoder and residual updates | Preserve/extract | optimized and memory efficient |
| Block-eigen operator | Preserve/migrate | already separated conceptually |
| Scalar task/seed/status helpers | Extracted for canonical BayesC; matched in current BayesR | reusable and validated without rerouting BayesR |
| Recent chain-local RNG patterns | Extract | canonical reproducibility model |
| Persistent thread-local distributions | Correct | localized reproducibility risk |
| Raw `stblr_raw_v1` formatter | Preserve during migration | stable current public contract |
| Rcpp result construction in kernels | Move to binding boundary | prevents language-neutral core |
| Annotation stick-breaking math | Extract carefully | shared across operators |
| Group normalization | Preserve as legacy behavior | not identical to future hierarchy normalization |
| Legacy multivariate sampler | Rewrite | unresolved statistical contract |
| Positional MT output | Replace | unsuitable for extension |
| Obsolete duplicate MT covariance helpers | Remove after replacement | maintenance burden |

---

## 8. Shared infrastructure readiness

| Infrastructure | Current readiness | Phase |
|---|---|---|
| Typed specifications | Implemented and validated | Phase 1 |
| Typed result vocabulary | Implemented; CSR BayesC active | Phase 1–2 |
| Language-neutral error boundary | Active for canonical CSR BayesC | Phase 1–3 |
| Chain/RNG contract | Shared vocabulary active for CSR BayesC; migrated CSR BayesR preserves the proven task/seed semantics in its operator-aware core | Phase 2–5 |
| Posterior accumulation | Retained-iteration predicate shared; model-specific containers/finalization deferred | Phase 2–4 |
| Probability utility layer | BayesC binary and BayesR categorical semantics deliberately remain model-specific | Deferred until a second exact use exists |
| Scale-policy interface | Typed CSR BayesC controls active; broader extraction deferred | Phase 4–6 |
| Covariance-policy interface | Scalar CSR BayesC active, MT unresolved | Phase 3 and Phase 8 |
| CSR operator ownership | Shared borrowed immutable contract active | Phase 1–3 |
| BED operator ownership | Mixed core/binding | Phase 7 |
| Factor core | Not implemented | Phase 10–11 |

---

## 9. Performance and memory expectations

| Migration type | Runtime expectation | Memory expectation |
|---|---|---|
| Boundary-only refactor | Approximately unchanged | Unchanged |
| Shared utility extraction | No material regression | No material regression |
| RNG correction | Small acceptable variation only if justified | Unchanged |
| Operator migration | No material regression | No material regression |
| New MT implementation | Establish reference first; optimize before production | Explicit memory target |
| New factor implementation | Benchmark against scientific alternatives | Output-dependent and documented |

A working backend should not be replaced if the migrated implementation has an unexplained runtime or peak-memory regression.

---

## 10. Optional memory-heavy outputs

| Output | Default recommendation |
|---|---|
| Marker posterior mean | On |
| Marker PIP | On |
| Global parameter traces | On with thinning |
| Full marker samples | Off |
| Per-chain marker summaries | Off unless requested |
| Marker × component probability | Model-dependent; configurable |
| Final residual vector | Off unless diagnostic |
| Full CPO diagnostics | Off unless requested |
| Factor chain loadings | Retain only as needed for alignment/summary |

---

## 11. Trait-design contracts

| Design | Effect covariance | Residual covariance | Support |
|---|---:|---:|---|
| Same individuals and design | Full allowed | Full allowed | Planned canonical MT |
| Same design with structured missingness | Conditional | Conditional | Requires explicit sufficient statistics |
| Partially overlapping studies | Prior covariance possible | Only with overlap information | Future |
| Independent studies | Prior covariance possible | Usually diagonal | Planned restricted mode |
| Evidence matrix | Factor covariance | Observation-model residual | Separate factor frontend |

---

## 12. Replacement readiness checklist

A model/operator route may become canonical when:

- typed specification exists;
- core/binding boundary is clear;
- statistical reduction tests pass;
- reproducibility tests pass;
- runtime has no unexplained material regression;
- peak memory has no meaningful regression;
- output semantics are documented;
- legacy implementation has a removal plan;
- the extension path is clearer than before.
