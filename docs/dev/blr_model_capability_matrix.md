# Unified BLR Framework: Model, Architecture, and Migration Capability Matrix

**Status:** Phase 16A experimental packed-BED BayesC routes explicitly disposed
**Date:** 2026-07-18
**Target location:** `docs/dev/blr_model_capability_matrix.md`

## Current framework status

Canonical ordinary CSR: BayesC, BayesR, SBayesRC, fixed-prior BayesC, group
BayesC, and learned-annotation BayesC. Canonical scheduled: ordinary-CSR BayesC
and public scheduled packed-BED BayesC. Experimental/sparse BED BayesC and
block-eigen remain audited/protected and noncanonical. The authoritative public
default multivariate route is supported legacy and noncanonical; alternative
multivariate implementations retain separate experimental/research statuses.
Historical phase descriptions below retain their original point-in-time scope.

Packed-BED BayesR remains noncanonical and production-unchanged. Phase 13A
completed its route/model/component, scheduler/RNG, aggregation/schema,
deterministic-reference, and extraction-boundary audit.

Phase 13C activates the typed component specification, immutable genotype view,
per-chain execution context, callable numerical core, typed chain result, and
binding-neutral progress-event boundary. Phase 13A trajectories, task dispatch,
inline aggregation/conversion, public route, and schema remain unchanged.

Phase 13D activates `BedBayesRExecutionResult`, one native aggregation callable,
and one named binding converter. Fit-local/thread-independent reproducibility,
logical-chain RNG, adaptive scheduling, immutable genotype ownership, public
routing and schemas remain exact. The route is migrated, not yet canonical.

Phase 13E marks public packed-BED BayesR canonical with permanent Phase 13A
fixtures and the Phase 13D runtime/completed-fit-RSS/I/O baseline. Migration
scaffolding is absent; BayesRC and experimental/sparse BED BayesC remain
unchanged and noncanonical.

Phase 14A leaves packed-BED BayesRC production unchanged and noncanonical while
auditing its public route, full-sweep model, ordered components, annotation
alignment, probit-stick policy, coefficient updates, RNG ownership,
deterministic references, aggregation/schema and future extraction seam.

Phase 14B mechanically extracts the unchanged per-chain full-sweep execution to
one guarded implementation header. Probit-stick, latent/coefficient updates,
chain RNG and ownership remain exact; alignment, dispatch, aggregation and R
conversion remain inline. The typed callable boundary is deferred to Phase 14C.

Phase 14C activates binding-neutral typed component, annotation, coefficient-prior,
borrowed genotype/annotation view, per-chain context/result, and callable-core
contracts. A stateless normal-probability interface retains exact Rmath semantics;
alignment, dispatch, aggregation, and conversion remain adapter-owned.

Phase 14D adds the typed aggregate result, one native aggregation path, one final
marker-prior recomputation path, and one named binding converter. Full sweeps,
probit sticks, latent/alpha updates, logical-chain RNG, Phase 14A references,
the public route, and schemas remain exact. BayesRC is migrated, noncanonical,
and ready for canonicalization.

Phase 14E marks public packed-BED BayesRC canonical. The typed chain and
aggregate architecture, singular final-prior and converter paths, full-sweep
policy, Rmath-backed probability boundary, permanent Phase 14A fixtures, and
Phase 14D runtime/completed-fit-RSS/I/O baseline are permanent. Unsupported
annotation and scheduling policies remain unchanged.

Phase 15A classifies packed-BED family consolidation opportunities. Exact task
mapping and seed policy are candidates for narrow sharing; immutable genotype
ownership, dispatch/failure/timing, converter metadata, reference testing, and
benchmark reporting are convention-level opportunities. Canonical numerical
cores and all model-specific statistical contracts remain unchanged, and
production consolidation has not started.

Phase 15B makes the proven common layer canonical: exact task indexing, exact
logical-chain seeds, and the compatible BayesC/R immutable genotype view.
BayesRC representation-specific ownership and all model-specific dispatch,
failure/timing/progress, numerical, aggregation, conversion, scheduler, and
probability contracts remain separate. Shared test/reference and benchmark
conventions are active; all public model classifications remain canonical.

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
| ST BayesC scheduled packed BED | scalar | binary | global_binary | unit | scalar | SNP-major packed BED | Canonical corrected logical-chain RNG with borrowed typed per-chain context/result, callable core, typed aggregate result, one native aggregation path, and one named binding converter | Fit-owned immutable genotype borrowed; task dispatch retained in adapter; permanent Phase 11B references and fit-local/thread-independent reproducibility exact; public route/schema unchanged; runtime/completed-fit-RSS/I/O baselines established; migration scaffolding removed | Public multichain route canonical. Explicitly experimental lower-level scheduled single-chain and sparse routes are retained unchanged: the former is a deterministic historical oracle and the latter preserves a distinct `null_update_prob` scheduling experiment. Neither is canonical or supported for general use. Packed-BED BayesR and BayesRC remain canonical separate implementations. |
| ST BayesC fixed marker prior | scalar | binary | fixed_marker/global | fixed_marker | scalar | CSR | Canonical typed implementation with borrowed `pi_marker`/`vb_multiplier`, callable core, typed result, one binding converter, and one wrapper aggregation path | Canonical; public route/schema unchanged; permanent exact references and runtime/completed-fit-RSS baselines active; trait restriction preserved; migration scaffolding removed | Stable canonical fixed-prior reference |
| ST BayesC + `selection_s` | scalar | binary | global_binary | maf_s | scalar | CSR | Current | Preserve and migrate | Good composite-scale test |
| ST group BayesC | scalar | binary | group | group | scalar | CSR | Canonical typed implementation with borrowed zero-based mapping, active group policy, callable core, typed result, one binding converter, and one wrapper aggregation path | Canonical; public route/schema and normalization semantics unchanged; permanent exact references and runtime/completed-fit-RSS baselines active; unsupported cases preserved; migration scaffolding removed | Stable canonical group reference; group normalization and mutable state remain specific |
| ST learned annotation BayesC | scalar | binary | annot_logit | learned annotation scale | scalar | CSR | Canonical typed implementation with borrowed annotation design, active learned policy, callable core, typed result, one binding converter, and one wrapper aggregation path | Canonical; public route/schema unchanged; permanent exact references and canonical runtime/completed-fit-RSS baselines active; unsupported cases preserved; migration scaffolding removed | Stable canonical learned-annotation reference; centered logistic/exponential MH policy, proposals, bounds, and update frequency remain specific |
| ST BayesR CSR | scalar | mixture | global_dirichlet | component | scalar | CSR | Canonical | Typed borrowed context, operator-aware templated core, typed result, one converter; permanent exact references and runtime/memory baseline | Stable canonical scalar mixture reference |
| ST BayesR BED | scalar | mixture | global_dirichlet | component | scalar | SNP-major packed BED | Canonical typed component/genotype/context/core/chain-result/aggregate-result architecture with one native aggregation path and one named converter | Canonical; logical-chain RNG, adaptive scheduler, adapter-rendered progress, permanent Phase 13A fixtures, and Phase 13D runtime/completed-fit-RSS/I/O baseline active | Public route/schema unchanged; BayesRC is canonical separately; experimental/sparse BayesC remains noncanonical |
| ST BayesRC scheduled-chains BED | scalar | mixture | annotation probit stick | component | scalar | SNP-major packed BED | Canonical typed component/annotation/prior/views/context/core/chain-result/aggregate-result architecture; one aggregation, final-prior, and converter path; exact Rmath-backed probability interface | Canonical; ordered sticks, latent/alpha updates, logical-chain RNG, permanent Phase 14A references, and Phase 14D runtime/completed-fit-RSS/I/O baseline active | Public route/schema unchanged; alignment, dispatch, and optional genotype diagnostics remain adapter-owned; adaptive controls do not apply |
| ST SBayesRC CSR | scalar | mixture | annot_probit_stick | component/maf_s | scalar | CSR | Canonical typed borrowed CSR/annotation context with operator-aware templated core, typed result, and one ordinary-CSR converter | Ordered-probit and alpha behavior preserved; shared task/seed/status infrastructure active; public route/schema unchanged; permanent exact references and runtime/completed-fit-RSS baseline established | Canonical annotation-aware mixture reference |
| ST BayesRC BED | scalar | mixture | annot_probit_stick | component | scalar | BED | Alias inventory entry for the same canonical public full-sweep route above | Canonical through the scheduled-chains implementation row; not a second execution path | Permanent alias inventory only |
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

---

## 13. Phase 17A non-packed-BED status

| Family | Reachability | Current architecture | Reference strength | Disposition |
|---|---|---|---|---|
| ordinary and scheduled CSR BayesC | supported public | canonical | permanent exact | canonical |
| CSR BayesR and CSR SBayesRC | supported public | canonical | permanent exact | canonical |
| fixed-prior, group, and learned-annotation CSR BayesC | supported public | canonical | permanent exact | canonical |
| block-eigen BayesC/BayesR/SBayesRC | internal research | audited but noncanonical | smoke only | retain experimental |
| default legacy multivariate `mtblr` | authoritative supported public legacy | legacy and noncanonical | `updateB` corrected; explicit post-burn policy and accumulator-specific counts; historical and corrected deterministic references established | ready for guarded mechanical extraction; retain historical behavior as legacy evidence |
| CPG/Armadillo/eigen multivariate variants | public algorithm variants | legacy multivariate | smoke only | retain experimental |
| OpenMP CPG multivariate | public algorithm variant | legacy multivariate | worker-sensitive RNG risk | correct before migration |
| hybrid multivariate | native-only | legacy multivariate | none | retirement candidate |
| CSR OpenMP multivariate | native-only research | legacy multivariate | none | retain experimental pending evidence |

Phase 17A changes status metadata only. It does not alter numerical code,
routing, exports, schemas, or the canonical packed-BED classifications.

## 14. Phase 17B public multivariate audit status

| Route | Contract/reference status | Disposition |
|---|---|---|
| default `mtblr` | supported public legacy; authoritative; noncanonical; update control and retained-sample policy corrected; 3/3 corrected raw and 3/3 corrected formatted references within `1e-12`, exact structure; Phase 17B history retained | ready for mechanical extraction; migration not started |
| `mtblr_cpg` | historical redundant/defective progression variant | retired and removed in Phase 17G |
| `mtblr_cpg_arma` | historical redundant/defective Armadillo variant | retired and removed in Phase 17G |
| `mtblr_cpg_omp` | historical worker-sensitive branch | retired and removed in Phase 17G |
| `mtblr_eigen` | native-only legacy eigen evidence; unsupported | retain temporarily as internal research |
| `mtblr_cpg_omp_csr` | native-only local-CSR evidence; unsupported and noncanonical | retain temporarily as internal research |
| `mtblr_hybrid` | historical unused prototype | retired and removed in Phase 17G |

No `stblr_raw_v1` boundary, route, or production implementation changes in
Phase 17B.

## 15. Phase 17C corrected public multivariate status

All `B` updates honor `updateB`; post-burn begins at `nburn`; marker thinning is
post-burn-relative; and marker, covariance, and probability summaries use
distinct retained counts. Alternative multivariate classifications remain
unchanged, including the `mtblr_cpg_omp` worker-sensitive RNG risk.

## 16. Phase 17D lexical extraction status

| Route | Execution boundary | Result/conversion boundary | Status |
|---|---|---|---|
| default `mtblr` | one guarded lexical corrected single-fit implementation header; fit-local RNG and retained-count contract preserved | legacy 20-position finalization remains inline; R formatting unchanged | supported public legacy, noncanonical, migration in progress; typed context/result not active |

Phase 17C corrected references remain authoritative current behavior and Phase
17B fixtures remain immutable historical evidence. Alternative multivariate
and block-eigen dispositions are unchanged.

## 17. Phase 17E typed core-boundary status

| Route | Typed numerical boundary | Ownership | Remaining legacy boundary | Status |
|---|---|---|---|---|
| default `mtblr` | typed data/model/prior/execution/initial-state/result contracts and one callable binding-neutral core | dense summaries, models, sets, priors borrowed; mutable `b/B/E/pi` and core outputs owned | inline 20-position finalization and unchanged R formatting | authoritative supported public legacy; noncanonical; Phase 17C contract active; ready for typed finalization separation |

No native aggregation, named converter, chain/task framework, or generic data
operator is active. Alternative multivariate and block-eigen classifications
remain unchanged.

## 18. Phase 17F typed finalization and future-operator status

| Route | Core | Finalization | Compatibility adapter | Future operator status |
|---|---|---|---|---|
| default `mtblr` | one typed binding-neutral callable core | one typed binding-neutral posterior finalizer and owning finalized result | unchanged legacy 20-position positional adapter and R formatter | no operator active; canonical scalar CSR/block-eigen representation reuse and one operator per trait are required |

The route remains authoritative supported public legacy and noncanonical.
Multichain aggregation, a named R converter, a versioned raw schema, CSR/eigen
execution, and a generic operator abstraction remain inactive.
## Phase 17H storage note

All canonical scalar CSR families share `SparseLdCsrStorage`; ordinary CSR
BayesC additionally uses `SparseLdCsrView`. This is a storage/ownership
contract only. Dense corrected MT remains the sole public MT route, while the
local MT CSR and eigen implementations remain unsupported research.
## Phase 17I capability note

Trait-specific canonical MT CSR execution is validated internally against the
corrected dense oracle, including shared operators, different values/diagonals,
and independent sparsity. It is nonpublic and not yet a supported API.
| MT CSR BayesC | `mtblr_csr()` | supported public, canonical CSR | trait-specific canonical CSR views; serial Phase 17I core; named MT raw v1; no overlap likelihood |
Test portability: Phase 17J2 validates supported scientific/public capabilities
from both source and installed-package contexts without changing capability.
Phase 17K: scalar block-eigen BayesC, BayesR, and SBayesRC remain internal and
experimental. They share canonical block-filtered storage/view contracts; no
public selector, scheduled execution, LD-swap, or MT block-eigen route exists.
# Phase 17L internal capability

The internal MT BayesC block-eigen route supports shared and trait-specific operators, independent block boundaries, mixed filters, and trait-specific analysis sample sizes. It is serial, single-chain, internal-only, and does not model sample overlap. Public scalar and MT routes remain CSR/dense as previously documented.

# Phase 17M public capability

`mtblr_block_eigen()` is the supported public same-BED summary route for
serial multivariate BayesC. External GWAS/reference-panel combinations,
sample-overlap modeling, non-diagonal residual covariance, OpenMP, and
multichain execution remain unsupported.

# Phase 17N planned capability

An internal individual-level MT route now exists as `mtblr_bed_internal()`.
It is serial, one-chain BayesC with one shared packed genotype matrix, complete
finite centered pre-adjusted phenotypes, corrected joint patterns, canonical
full residual covariance, and a diagonal reduction mode. It has no public
adapter, missing-phenotype, covariate, CPO, LE/LD, OpenMP, or multichain support.

# Phase 17P public capability

`mtblr_bed()` is the supported public joint individual-level BayesC route. It
uses one BED-backed Glist and one shared individual set, standardized genotypes,
complete centered pre-adjusted phenotypes, full residual covariance by default,
and an explicit diagonal reduction mode. It remains serial and one-chain, with
no covariate fitting, missing phenotypes, CPO, LE/LD, predictions, OpenMP, or
multichain execution.

Phase 17Q does not change this capability matrix. It formalizes a future
BayesC-only chain-parallel route: each chain remains a complete joint MT fit;
no trait, marker, or set parallelism is specified.

Phase 17R makes chain-level multichain/OpenMP execution available through the
internal `mtblr_bed_chains_internal()` maintenance route. Public `mtblr_bed()`
remains serial and one-chain; no new statistical capability is introduced.

Phase 17S makes serial and static chain-level OpenMP multichain execution public
for individual-level MT BayesC. Full and diagonal residual covariance remain
supported; no convergence diagnostic or other model family is added.

Formal convergence diagnostics are contract-only in Phase 17T. Current
`bm_sd/dm_sd` remain chain-stability summaries rather than R-hat or ESS.
## Phase 17U convergence capability

| Route | Internal Tier 1 engine | Public diagnostics | Scope |
|---|---:|---:|---|
| MT packed BED multichain | yes | no | post-burn B/G/E diagonals |
| MT CSR/block eigen | no | no | not applicable |

Tier 2 covariance/probability and Tier 3 selected-marker traces remain
unsupported.

Phase 17V makes Tier 1 B/G/E-diagonal diagnostics public for MT packed BED.
Off-diagonal covariance, probability, and marker diagnostics remain unsupported.

Phase 17W formalizes, but does not implement, Tier 2 raw strict-lower B/G/E and
active/null/explicit-pattern diagnostics plus later Tier 3 selected-marker b/d.
