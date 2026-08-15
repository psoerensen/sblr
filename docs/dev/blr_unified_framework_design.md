# Unified BLR framework design

## 1. Status and scope

**Status:** `PHASE 2 IMPLEMENTED - READY FOR INDEPENDENT VERIFICATION`

This document defines a shared R/C++ architecture for three distinct analysis
modes:

1. one single-trait posterior;
2. several independent single-trait posteriors fitted in one invocation;
3. a genuinely joint multi-trait posterior.

This document defines the target architecture and records completed
checkpoints. Phase 1 implemented the shared production R specification,
validation, schema-v2 conversion, and formatting infrastructure described in
Section 31. Phase 2 adds the shared R global-marker, immutable operator-resource,
and likelihood-provider adapters and qualifies the existing native immutable
CSR, BED, and block-eigen views without changing a native call. The exact
implementation boundary is recorded in
[the Phase 2 checkpoint](blr_phase2_provider_operator_checkpoint.md). A
capability described only as proposed remains unimplemented; current support
is determined by executable public dispatch, source, schemas, tests, and
maintained current contracts.

The Phase 0 closeout was prepared from branch `master` at commit
`b930c723d2ebf7351b8d7509bf9f4752820d1a86`. Its central decision is that STBLR
and MTBLR must share infrastructure whenever the responsibility is the same,
while using separate statistical policies whenever the posterior target differs.
It must not create a second provider, operator, scheduler, retention, or result
stack for MTBLR.

The design covers:

- current-source inventory and migration boundaries;
- one resolved R specification;
- global marker, reusable operator-resource, and likelihood-provider
  contracts;
- dense, BED, CSR, and block-eigen likelihood operators;
- model, state, covariance, probability, and residual policies;
- deterministic logical-task scheduling and RNG ownership;
- raw schema version 2;
- formatted-fit responsibilities;
- reductions, tests, migration, and implementation phases.

Phases 1 and 2 did not change native interfaces, statistical kernels, posterior
targets, seed behavior, or retention behavior. Phase 2 multi-trait BED support
is representational only: it does not invoke the current MT sampler. Phase 3
scheduler/RNG work, heterogeneous multi-provider posterior likelihoods, and
corrected joint-MT capabilities remain proposed, and this document does not
approve a new statistical model merely by describing it.

## 2. Authority hierarchy

The authority rules in `docs/dev/README.md` and the repository root
`AGENTS.md` apply. For current executable behavior, evidence is ordered as:

1. public dispatch and executable R/C++ source;
2. generated registration and adapters;
3. current raw schemas and formatted-fit construction;
4. permanent tests;
5. maintained current developer contracts;
6. qualification records;
7. plans, research prototypes, and historical records.

For statistical correctness, the approved files under `docs/methods/` and an
independent derivation take precedence over existing implementation behavior.
The source audit in `docs/dev/methods_source_traceability_audit.md` identifies
known differences between those two kinds of authority. The validated MT
research design in `docs/dev/mtblr_covariance_design.md` supplies evidence for
future policies, but its standalone R samplers are not production
infrastructure.

When a maintained contract and current source disagree, the discrepancy must
be classified and resolved explicitly. A shared framework must not make an
incorrect current transition authoritative merely by generalizing it.

## 3. Terminology

### 3.1 Analysis mode

`analysis_mode` is statistical:

| Value | Meaning | Posterior structure |
|---|---|---|
| `single_trait` | One trait and one single-trait posterior | One scalar-effect model |
| `independent_traits` | Several single-trait posteriors resolved together | Product of traitwise posteriors |
| `joint_multitrait` | One multivariate posterior | Traits coupled through a declared multivariate prior or likelihood |

For independent traits,

$$
p(\boldsymbol{\beta}_1,\ldots,\boldsymbol{\beta}_T
\mid \mathcal D_1,\ldots,\mathcal D_T)
=
\prod_{t=1}^T
p(\boldsymbol{\beta}_t\mid\mathcal D_t).
$$

Joint MTBLR generally does not factorize this way. Coupling may be through
$V_b$, activity patterns, component states, covariance templates, shared
learned probability parameters, or an explicitly overlap-aware likelihood.

### 3.2 Execution mode

Execution policy is computational and is separate from `analysis_mode`. It has
two resolved fields:

| Field | Values | Meaning |
|---|---|---|
| `execution_mode` | `serial`, `parallel` | Whether more than one logical task may execute concurrently |
| `parallelization` | `none`, `chains`, `traits`, `trait_chains` | Which independent task axes may populate the parallel task pool |

Valid combinations are closed rather than inferred. `single_trait` permits
`none` or `chains`; `independent_traits` permits all four values;
`joint_multitrait` permits `none` or `chains`. `execution_mode = "serial"`
requires `parallelization = "none"`. Parallel traits do not create a joint
posterior. Within a joint chain, traits are statistically coupled and must not
be scheduled as independent trait tasks.

### 3.3 Operator resources, providers, blocks, and regions

- An **operator resource** owns one immutable genotype or LD representation,
  such as BED, dense cross-product, CSR, or block eigen. Several likelihood
  providers may reference the same resource without duplicating it.
- A **likelihood provider** owns one likelihood contribution, its scientific
  provenance, a nonempty set of traits, and a reference to an operator
  resource. Independent summary providers normally contain one trait;
  common-sample MT providers may contain several traits jointly.
- A **likelihood operator** evaluates the genotype cross-product or declared
  approximation owned by an operator resource.
- A **computational eigenblock** partitions an operator representation.
- A **statistical covariance region** assigns a prior covariance policy.
- A **marker set** is an explicitly declared grouping whose statistical or
  computational role must be named.
- An **annotation/source group** may change probability or scale policies.

These objects are not interchangeable. In particular, changing an eigenblock
boundary must not silently change the prior model.

## 4. Current architecture inventory

### 4.1 R entry points and resolution

| Responsibility | Current evidence | Classification | Unified-framework action |
|---|---|---|---|
| ST CSR public dispatch | `R/stblr-public.R`, `stblr_csr()` | Generalize behind a common interface | Retain as a compatibility wrapper that resolves a common specification |
| ST BED public dispatch | `R/stblr-public.R`, `stblr_bed()` | Generalize behind a common interface | Resolve to the same specification and provider contract |
| ST block eigen | `R/stblr-block-eigen.R`, `stblr_block_eigen()` | Generalize behind a common interface | Keep representation-specific preparation; remove model-resolution duplication |
| Annotation ST CSR | `R/stblr-csr-annot.R`, `stblr_csr_annot()` and provider-specific helpers | Model-specific policy plus duplicated wrapper infrastructure | Resolve probability/scale policies into common model namespaces |
| MT BED | `R/mtblr-bed.R`, `mtblr_bed()` | Generalize data preparation; replace statistical core | Reuse common sample/provider preparation and formatter envelope |
| MT CSR | `R/mtblr-csr.R`, `mtblr_csr()` | Generalize provider preparation; replace statistical core | Preserve validated alignment concepts, not the hybrid covariance transition |
| MT block eigen | `R/mtblr-block-eigen.R`, `mtblr_block_eigen()` | Generalize behind provider operators | Current reconstructed-dense route is an adapter, not the final operator abstraction |
| Common chain controls | `R/blr-unified.R`, `.blr_chain_controls()` | Shared already | Expand into `mcmc` and `compute` specification groups |
| ST model resolver | `R/blr-unified.R`, `.blr_resolve_st_model()` | Reusable concept, narrow implementation | Replace method-string logic with a validated model policy descriptor |
| MT model resolver | `.mtblr_resolve_public_method()` and helpers in MT wrappers | Duplicated | Merge semantic resolution without merging posterior policies |
| Convergence resolution | `R/blr-unified.R`, `.blr_convergence_controls()` and `R/mtblr-convergence.R` | Shared already in substance | Make it consume schema-v2 draw descriptors |
| ST logical seeds | `R/blr-unified.R`, `.blr_st_task_seeds()` | Reusable concept | Version the logical-task seed contract across all analysis modes |
| MT chain seeds | MT wrapper and native chain helpers | Generalize behind a common interface | Use the same task identity and seed resolver |
| ST raw validator/formatter | `R/sparse_ld_bed_helper.R`, `.validate_stblr_raw()`, `.as_stblr_fit()` | Replace with versioned common envelope | Keep strict failure on unsupported schemas |
| MT raw validator/formatter | `R/mtblr-csr.R`, `.validate_mtblr_raw()` and MT formatter | Replace with versioned common envelope | Remove separate ST/MT schema stacks |
| Final fit normalization | `R/blr-unified.R`, `.blr_finalize_fit()` | Reusable concept | Make it a schema-v2 formatter, not a semantic repair layer |
| Consistency checking | `R/check-stblr-backend-consistency.R`, `check_stblr_consistency()` | Reusable concept | Validate explicit schema-v2 namespaces and dimensions |

Current public wrappers use the canonical controls `nit`, `nburn`, `nthin`,
`nchains`, `ncores`, and `seed`, but validation is repeated in older helpers.
Some wrappers collect `...` and forward with `do.call()`. Phase 3A demonstrated
that exact pre-dispatch name validation is essential: unnamed, duplicate, and
partially matched arguments must never reach a scientific backend.

### 4.2 Data and likelihood representation

| Responsibility | Current evidence | Classification | Unified-framework action |
|---|---|---|---|
| BED decoding | `src/st_bed_decode.h` and BED data helpers | Reusable unchanged after a view adapter | One immutable packed-genotype provider view |
| ST CSR operator | `src/st_ld_operator.h`, `CsrOperator` | Shared already for scalar kernels | Promote its method contract, not its ST-specific name |
| Retained block eigen | `src/blr_block_eigen.h`, `src/st_ld_operator.h`, `BlockEigenDispatchOperator` | Reusable representation | Expose as a provider-local operator with retained residual state |
| Dense reconstructed eigen | block-eigen R/native preparation | Reusable adapter | Represent explicitly as a dense approximate cross-product operator |
| Current MT CSR data view | `src/blr_mt_default_core_impl.h` data-view templates | Generalize | Split immutable provider data from MT state policy |
| Current MT block eigen | reconstructed-dense MT dispatch | Replace as primary design | Use the same provider-specific retained/full-rank operator abstraction as ST |
| Marker alignment | ST/MT R preparation and `marker_policy` checks | Generalize behind a common interface | One global marker universe and explicit local maps |
| Multiple providers | research prototype under `tests/research/mtblr_covariance/` | Validated design, not current production | First-class provider collection |

`CsrOperator` and `BlockEigenDispatchOperator` already expose the beginnings of
a useful representation-neutral method set: diagonal access, corrected marker
right-hand side, residual update, fitted/operator quadratic, projected score
dot product, and residual reconstruction. Those methods should become the
operator concept. They do not imply that ST and MT marker transitions are the
same.

### 4.3 Samplers and statistical state

| Responsibility | Current evidence | Classification | Unified-framework action |
|---|---|---|---|
| ST CSR BayesC loop | `src/blr_csr_bayesc_core_impl.h` | Generalize around policies | Retain validated conditionals; separate loop mechanics from probability/scale policy |
| ST CSR BayesR loop | `src/blr_csr_bayesr_core_impl.h` | Generalize around policies | Reuse scheduling, retention, and operator access |
| Annotation probability/scale policies | annotation C++ policy types and provider helpers | Model-specific policy | Preserve as explicit `ProbabilityPolicy` and `MarkerScalePolicy` variants |
| MT activity-pattern state | `src/blr_mt_bayesr_types.h`, `src/blr_mt_bayesr_kernel_impl.h` and MT default core | Model-specific policy | Retain verified state geometry only after covariance transition replacement |
| Scalar marker variance | ST kernels | Model-specific covariance policy | `ScalarMarkerVariancePolicy` |
| Matrix marker covariance | current MT cores and `src/mtblr.cpp` | Replace | Implement the approved sampled-$V_b$ policy from the MT research design |
| Residual scalar variance | ST kernels/operator helpers | Reusable policy | `IndependentResidualVariancePolicy` |
| MT residual covariance | MT BED and summary routes | Model- and likelihood-specific policy | Separate full common-sample, diagonal summary, and future overlap-aware policies |
| Current MT set loop | `run_mt_bayesc_core_impl()` in `src/blr_mt_default_core_impl.h` | Replace | Computational traversal must not control global prior update frequency |
| `sampleBset()`/`sampleB()` heuristic | `src/mtblr.cpp` and MT cores | Retire | Never use as authoritative sampled $V_b$ |
| `sampleB_latent()` inside every set | MT default/BED cores | Replace | One scientifically valid covariance transition per completed iteration |
| `sampleBset_old()` and superseded paths | `src/mtblr.cpp` | Retire after reference search | Unreachable legacy code is not a capability |

The current MT core calls `sampleBset()` and then overwrites its result with
`sampleB_latent()` before marker updates, repeats an all-marker covariance draw
once per set, and later replaces that matrix with a heuristic covariance for
traces and output. `docs/dev/mtblr_covariance_design.md` establishes that this
is not a reusable covariance policy. The unified framework must introduce one
authoritative sampled covariance state and must not preserve this hybrid.

### 4.4 Scheduling, retention, results, and tests

| Responsibility | Current evidence | Classification | Unified-framework action |
|---|---|---|---|
| Static OpenMP task scheduling | scalar CSR cores and `src/blr_mt_bed_chains_execution_impl.h` | Reusable unchanged in principle | Schedule logical tasks, never statistical coordinates |
| Task-local `std::mt19937` | scalar cores and MT chain runners | Shared already in principle | Central versioned seed derivation and ownership contract |
| Worker diagnostics | current memory/worker and Phase 3A sampler-local diagnostics | Generalize | Distinguish requested, configured, actual team size, and worker IDs |
| Retention/thinning | repeated kernel-local allocation and loops | Generalize behind a common interface | One retention plan and shape descriptor |
| Native result structs | `src/blr_result.h`, `BlrResult` | Reusable starting point, currently narrow | Extend to schema-v2 scientific namespaces and explicit dimensions |
| Resolved native spec | `src/blr_spec.h`, `ResolvedSpec` | Reusable starting point, currently CSR BayesC-only | Replace restrictive enums with validated general descriptors |
| Generated Rcpp interfaces | `R/RcppExports.R`, `src/RcppExports.cpp` | Generated-interface consequence | Change only when vertical slices migrate |
| Raw schema tests | `test-stblr-raw-schema.R` and MT raw-contract tests | Reusable test ownership | Add version-2 structural and semantic tests |
| Public contract tests | `test-blr-unified-public-contract.R`, model/interface tests | Reusable test ownership | Add resolved-spec and migration assertions |
| Operator reductions | `test-blr-operator-reductions.R` and block-eigen tests | Reusable oracle | Extend across providers and analysis modes |
| RNG/reproducibility | `test-blr-unified-reproducibility.R` and MT chain tests | Reusable oracle | Key by logical task rather than backend traversal |
| Research MT covariance tests | `tests/research/mtblr_covariance/` | Scientific oracle only | Keep independent from production kernels |

## 5. Confirmed reuse opportunities

The following responsibilities should have one implementation family across
ST and MT:

1. exact R argument validation and resolution;
2. marker identifiers, alleles, and local-to-global maps;
3. immutable operator-resource preparation, provider preparation, and
   provenance;
4. BED decoding;
5. dense, CSR, and block-eigen operator representations;
6. logical-task enumeration, deterministic seeds, and OpenMP scheduling;
7. burn-in, thinning, retention, and memory planning;
8. convergence trace capture;
9. numerical and worker diagnostics;
10. raw-result dimensioning and construction;
11. formatted-fit construction and consistency checking.

Reuse means a shared contract and implementation, not a universal marker
kernel. The likelihood operator should supply representation-neutral linear
algebra; a statistical policy should decide which conditional distribution is
being evaluated.

## 6. Components requiring replacement or retirement

| Component | Action | Reason |
|---|---|---|
| Separate `stblr_raw` and `mtblr_raw` version-1 envelopes | Replace through a versioned migration | They encode similar quantities with different shapes and overloaded names |
| Positional or permissive `...` forwarding | Retire | Partial/positional matching can change scientific parameters |
| Current MT covariance hybrid | Replace | Marker transition, next-iteration state, trace, and returned covariance are not one coherent sampled parameter |
| Set-driven global covariance draws | Retire | Computational partition changes statistical update frequency and RNG path |
| Heuristic covariance replacing a draw | Retire | A descriptive estimate must not overwrite sampled $V_b$ |
| Separate MT block-eigen stack | Do not build | Provider-specific operators should be shared |
| Ambiguous `pis`, `pi`, and `pim` semantics | Replace with explicit names | They mix prior state mass, traces, and model-dependent meanings |
| Capability decisions from `.blr_model_capability_matrix()` | Retire | The matrix is known to be stale |
| Unreachable legacy MT covariance functions | Retire after active-reference audit | Source existence is not support |

The current scalar ST kernels and current MT state calculations should be
treated differently: validated scalar conditionals are migration candidates;
the MT state geometry is a reference for a new vertical slice, while its
covariance scheduling and output flow are not.

## 7. Proposed shared conceptual architecture

```text
Public R interfaces
        |
        v
Resolved BLR specification
        |
        v
Global marker map + operator resources + provider collection
        |
        v
Shared likelihood-operator layer
        |
        v
Shared chain engine and scheduler
        |
        +-- single/independent-trait statistical policy
        +-- joint-multitrait statistical policy
        |
        v
Versioned raw-result builder
        |
        v
Formatted fit
```

### 7.1 Layer responsibilities

| Layer | Owns | Must not own |
|---|---|---|
| Public R interface | User vocabulary, deprecation messages, exact argument capture | Native shapes or posterior shortcuts |
| R resolver | Defaults, aliases, cross-field validation, capability resolution | Sampling or result reinterpretation |
| Global marker map | Canonical markers, alleles, indices, presence masks | Effect estimates |
| Operator resources | Immutable BED/dense/CSR/eigen storage and representation provenance | Trait likelihoods or priors |
| Provider collection | Likelihood groups, trait sets, local maps, resource references, and scientific provenance | Cross-trait prior state |
| Likelihood operator | Linear algebra for one operator resource | BayesC/BayesR state probabilities |
| Statistical policy | Joint model, state space, conditionals, parameter updates | File formats or scheduling |
| Covariance/residual policy | $V_b$/$V_e$ prior, update, state, and summaries | Operator traversal order |
| MCMC engine | Update order, iteration lifecycle, policy calls | R objects and presentation |
| Scheduler/RNG | Logical tasks, workers, seeds, task-local RNG | Posterior target |
| Retention | Burn-in, thinning, draw shape, chain identity | Posterior summaries that replace draws |
| Raw builder | Typed scientific arrays and metadata | Dimension dropping or compatibility aliases |
| Formatter | User conveniences, tables, explicit simplification | Overwriting sampled parameters |
| Diagnostics/provenance | Runtime, safeguards, approximation and build lineage | RNG consumption or target changes |

### 7.2 Composition examples

| Model | State/probability policy | Effect/covariance policy | Residual policy |
|---|---|---|---|
| ST BayesC | binary state + scalar inclusion probability | scalar slab variance | scalar residual variance |
| ST BayesR | null-plus-component categorical state | scalar base variance + $\gamma_kq_j$ scale | scalar residual variance |
| Annotation ST | selected ST state policy + annotation `P` and/or `Q` policy | scalar slab variance with declared marker multiplier | scalar residual variance |
| Independent traits | one ST policy instance per trait | independent scalar states | one residual state per trait |
| Cheng MT-BayesC$\Pi$ | joint activity-pattern categorical state | completed latent vector + one sampled $V_b$ | full common-sample or diagonal summary policy |
| Shared-component MT-BayesR | joint activity pattern and shared component | completed vector with $\gamma_kq_jV_b$ | declared MT residual policy |
| Future regional covariance | existing state policy plus region assignment | persistent $V_{b,r}$ policy updated once per iteration | unchanged likelihood-appropriate policy |
| Future covariance templates | template/component state policy | fixed or learned template library | unchanged likelihood-appropriate policy |

Trait-specific-component MT-BayesR is not covered by the shared-component
row. The architecture reserves a state-policy extension point without claiming
that its posterior has been derived.

For common-sample Cheng MT-BayesC$\Pi$, one multi-trait likelihood provider
references one shared BED operator resource and owns the aligned phenotype
matrix and joint residual policy. For independent traits using the same BED,
several singleton-trait providers may reference that same resource. Thus shared
infrastructure does not require duplicated genotype decoding, and a full
common-sample $V_e$ is not incorrectly factorized into marginal providers.

## 8. Proposed R architecture

### 8.1 Public strategy

Existing public wrappers should initially remain as exact-name adapters. They
should call one internal resolver and then one execution entry point. A future
structured public function may be introduced after the internal contract is
stable, but it does not block Phase 1.

Every wrapper must:

- capture supplied arguments before forwarding;
- reject unnamed, empty, `NA`, duplicate, ambiguous, and unknown names;
- use exact matching only;
- resolve aliases once and record the migration;
- produce one immutable resolved specification;
- avoid backend-specific `do.call()` argument binding.

No compatibility alias may silently change a posterior parameter. Deprecated
aliases should either be transformed deterministically with a warning or fail
with a migration message.

### 8.2 R-side objects

R should own:

- `BlrResolvedSpec`: a uniquely named, versioned list;
- `GlobalMarkerMapSpec`: marker/allele identity and provider maps;
- `OperatorResourceSpec`: immutable BED/dense/CSR/eigen representation and
  approximation descriptors;
- `ProviderSpec`: trait set, likelihood regime, scientific provenance, local
  map, and operator-resource reference;
- `ModelPolicySpec`: state, probability, scale, covariance, and residual policy;
- `McmcSpec`, `ComputeSpec`, and `OutputSpec`;
- validation error messages and capability checks;
- conversion from native typed buffers to raw schema version 2;
- user-facing formatting and migration aliases.

R should not construct or repair posterior draws after native execution.

## 9. Proposed C++ architecture

### 9.1 Minimal concepts

| Concept | Recommended form | Responsibility |
|---|---|---|
| `GlobalMarkerMap` | validated immutable data structure/view | Global IDs, alleles, provider-local indices, presence masks |
| `OperatorResource` | immutable tagged resource/view | Dense/BED/CSR/block-eigen storage and provenance, reusable by providers |
| `ProviderCollection` | immutable vector of likelihood-provider views | Trait sets, resource references, local maps, and overlap groups |
| `LikelihoodOperator` | tagged variant plus hot-path templated adapter | Operations on an `OperatorResource` |
| `EffectState` | plain task-local structure | Latent, realised, component, and pattern states |
| `StatePolicy` | compile-time policy where hot; small runtime descriptor at dispatch | State enumeration, prior mass, conditional update |
| `ProbabilityPolicy` | model-specific policy | Fixed/learned state probabilities and their parameters |
| `MarkerScalePolicy` | model-specific policy | Unit, MAF-S, annotation, external $q_j$, component scale |
| `MarkerCovariancePolicy` | model-specific policy | Scalar variance, sampled $V_b$, future regional/template covariance |
| `ResidualPolicy` | model/likelihood-specific policy | Scalar, diagonal, full common-sample, future overlap-aware residual state |
| `ChainState` | plain task-local structure | Scientific parameter state and operator residual state |
| `ChainScheduler` | shared nonstatistical component | Logical tasks, static worker assignment, lifecycle |
| `RngPolicy` | shared task-local component | Versioned seed derivation and engine ownership |
| `RetentionPolicy` | shared component | Burn-in, thinning, typed draw buffers |
| `DiagnosticsCollector` | task-local plain buffers | Runtime, counters, safeguards, actual worker use |
| `RawResultBuilder` | serial post-worker builder | Validate shapes and construct Rcpp objects after parallel regions |

### 9.2 Runtime versus compile-time boundaries

Use runtime descriptors in R and at the native boundary because provider
representation and model family are user-resolved. Dispatch once per provider
or logical task. Inside marker loops, use templates or `std::variant` visitation
resolved outside the innermost loop so virtual calls are not paid per marker.

Use compile-time policies for validated hot conditional kernels where this
materially simplifies and speeds code. Use simple structures for dimensions,
maps, state, controls, and diagnostics. Do not construct a deep inheritance
hierarchy.

An operator must not know whether a marker is null. A state policy must not know
whether its cross-product is CSR or block eigen. A covariance policy must not
be called once per eigenblock unless the statistical model explicitly defines
one covariance per block.

### 9.3 Operator concept

The common operator surface should cover, with representation-specific state
where necessary:

- marker diagonal on the declared scale;
- marker score/right-hand-side calculation;
- apply effect difference to maintained residual state;
- rebuild or validate residual state;
- operator application and quadratic form;
- residual and fitted quadratics with explicit exact/approximate semantics;
- block traversal that is computational only;
- immutable provenance and approximation metadata.

BED may expose streaming genotype operations rather than materializing
$X^\top X$, while satisfying the same statistical requests. Dense and CSR
operators maintain marker-scale residual scores. Retained block eigen maintains
reduced-coordinate residual state. These different representations share a
contract, not necessarily one storage layout.

## 10. Provider and operator-resource contract

An operator resource $o$ owns an immutable genotype or LD representation

$$
\mathcal O_o
=
\left(
\mathcal M_o,
\mathcal G_o,
\mathcal A_o,
\mathcal P_o
\right),
$$

where $\mathcal M_o$ is the local marker set, $\mathcal G_o$ is BED storage,
a cross-product operator, or its declared approximation, $\mathcal A_o$
contains allele/coding information, and $\mathcal P_o$ records provenance.

A likelihood provider $d$ references one operator resource and supplies a
nonempty ordered trait set $A_d$. Write

$$
\mathcal D_{d,A_d}
=
\left(
\mathcal Y_{d,A_d},
\mathcal O_{o(d)},
N_d,
\mathcal L_d,
\mathcal P_d
\right),
$$

where $\mathcal Y_{d,A_d}$ contains the phenotype or summary sufficient
statistics, $\mathcal L_d$ declares the likelihood regime and residual/error
contract, and $\mathcal P_d$ records provider provenance.

For an independent marginal-summary provider, $A_d=\{t\}$ and
$\mathcal Y_{d,A_d}$ includes

$$
\mathbf s_{dt}=X_{dt}^\top\mathbf y_{dt},
\qquad
C_{dt}=X_{dt}^\top X_{dt},
$$

or a declared approximation on the cross-product scale. For a common-sample
joint MT provider, $A_d$ may contain several traits, $\mathcal Y_{d,A_d}$ may
contain the aligned phenotype matrix, and its likelihood is evaluated jointly
under the declared $V_e$. It must not be decomposed into traitwise providers
when a full residual covariance couples those traits.

Separating operator resources from likelihood providers prevents duplicated
BED decoding or LD storage when several traits or providers use the same
genotype representation.

### 10.1 Required provider metadata

Each provider owns:

- `provider_id`, a nonempty ordered `trait_ids` set, and an
  `operator_resource_id`;
- likelihood regime and declared use of the referenced operator resource;
- aligned phenotypes or score vectors and any phenotype cross-products
  required by its residual policy;
- $N_d$ and explicit cross-product/covariance scale;
- provider-local marker IDs and local-to-global map;
- effect and non-effect alleles and allele orientation;
- genotype coding, centering, and standardization;
- phenotype/effect scale;
- source and target population;
- operator/LD provenance;
- approximation status and retained-rank metadata;
- overlap or summary-error group declaration.

The global marker universe is the union required by the resolved analysis. A
missing provider marker supplies no likelihood term. It is not a zero effect or
a negative observation.

### 10.2 Supported operator representations

One operator resource may use:

- individual-level BED access;
- dense cross-product;
- CSR sparse LD/cross-product;
- full-rank block eigen;
- retained-rank block eigen.

Providers and resources may differ in $N_d$, marker coverage, LD population,
block boundaries, eigenvectors, eigenvalues, and retained rank. A full-rank
block-eigen representation is exact for its declared block-diagonal operator.
It equals an original dense cross-product only when omitted cross-block entries
are zero or are explicitly excluded by the declared likelihood. A retained-
rank representation equals its reconstructed approximate operator, not
automatically the unfiltered dense operator.

### 10.3 Multiple providers and overlap

For independent singleton-trait providers estimating the same effect
$\mathbf b_t$,

$$
p(\{\mathcal D_{dt}\}_d\mid\mathbf b_t)
=
\prod_d p(\mathcal D_{dt}\mid\mathbf b_t).
$$

Operator contributions may be accumulated without materializing
$\sum_d C_{dt}$. Combination requires aligned markers and alleles, compatible
genotype and effect scales, compatible phenotype scale, a common effect
interpretation, and independent summary-estimation errors.

Three likelihood regimes must be explicit:

1. common-sample individual-level likelihood;
2. heterogeneous independent-provider summary likelihood;
3. overlap-aware summary likelihood with declared cross-provider error
   information.

Marginal provider-specific block-eigen decompositions do not identify
cross-provider summary-error covariance. Overlap-aware likelihoods are a later
gated model, not an operator flag.

For common-sample MT data, one multi-trait provider may reference one shared
BED resource and evaluate the joint residual likelihood directly. For future
overlap-aware summaries, a likelihood group may couple several providers; that
group is a statistical object above the marginal operator resources and must
declare its cross-provider error information.

## 11. Scheduling and RNG contract

### 11.1 Logical task graph

| Analysis mode | Logical task | Permitted outer parallelism |
|---|---|---|
| `single_trait` | `chain` | `chains` |
| `independent_traits` | `trait × chain` | `traits`, `chains`, or `trait_chains` |
| `joint_multitrait` | `chain` | `chains` |

All traits inside a joint chain share state and are updated by the joint
statistical policy. They are not independent OpenMP tasks. Provider or
eigenblock traversal may be parallelized later only after proving that the
reduction order, RNG contract, and posterior transition remain unchanged.

When `ncores < number_of_logical_tasks`, tasks are assigned using deterministic
static scheduling. A worker may run several logical tasks sequentially, but
each task retains its own state and RNG.

### 11.2 Seed derivation

Seed-contract version 1 is backend-independent and uses unsigned modular
arithmetic. The public user seed is an integer in $[0,2^{32}-1]$; the chain
index is zero-based and in the same range. Analysis-mode codes are
`single_trait = 1`, `independent_traits = 2`, and `joint_multitrait = 3`.
`single_trait` uses the UTF-8 sentinel `sblr:single_trait`,
`joint_multitrait` uses `sblr:joint_multitrait`, and `independent_traits` uses
the stable trait ID itself.

The 64-bit FNV-1a hash of UTF-8 bytes is exact:

```text
h = 0xcbf29ce484222325
for byte in utf8(identity):
    h = (h xor uint64(byte)) * 0x00000100000001b3 mod 2^64
```

The exact SplitMix64 output transform is:

```text
splitmix64(x):
    z = x + 0x9e3779b97f4a7c15 mod 2^64
    z = (z xor (z >> 30)) * 0xbf58476d1ce4e5b9 mod 2^64
    z = (z xor (z >> 27)) * 0x94d049bb133111eb mod 2^64
    return z xor (z >> 31)
```

Seed-contract version 1 derives one native seed as follows:

```text
x = uint64(user_seed)
x = splitmix64(x xor uint64(1))
x = splitmix64(x xor uint64(analysis_mode_code))
x = splitmix64(x xor fnv1a64_utf8(identity))
x = splitmix64(x xor uint64(zero_based_chain_index))
native_seed = uint32((x xor (x >> 32)) & 0xffffffff)
```

All multiplication and addition are modulo $2^{64}$. The final seed is passed
unchanged to `std::mt19937::result_type`; zero is valid and is never replaced.
The fixed vectors in
`tests/research/blr_framework_contract/blr_contract_fixtures.R` are part of the
contract.

An explicit resolved `task_seeds` table bypasses derivation and supplies the
final native seeds. Its dimensions are `chain` for `single_trait` and
`joint_multitrait`, and `trait × chain` for `independent_traits`. A convenience
`chain_seeds` vector may be accepted by public independent-trait wrappers, but
it remains combined with each stable trait ID through the versioned derivation;
it is not itself a final task-seed table. The resolved task-seed table is
recorded in provenance. Changing any derivation rule requires a new seed-
contract version.

Current wrappers do not yet use this derivation. Phase 1 preserves each
wrapper's source-traced legacy seed rule while it resolves the new
specification. Migration to seed-contract version 1 is a later explicit
checkpoint with trajectory tests and release notes. A changed deterministic
trajectory is not, by itself, a changed posterior target.

### 11.3 RNG ownership and reproducibility

- Each logical task owns one C++ RNG engine.
- R RNG is never called inside worker code.
- No mutable RNG is shared across workers.
- Diagnostics and provenance consume zero random draws.
- Disabled optional policies are zero-RNG no-ops.
- Serial and parallel execution of the same logical tasks must produce
  identical task-level results under the declared contract.
- Changing execution mode must not change retained-draw semantics or the
  posterior target.

If a policy intentionally changes draw order, the result is a versioned
scientific/RNG migration, not an implementation detail.

### 11.4 Worker diagnostics

Record separately:

- requested cores;
- configured task workers;
- OpenMP availability and maximum available threads;
- actual team size inside the relevant sampler region;
- worker identifier assigned to each logical task;
- per-task and dispatch runtime.

A worker identifier is not a team size. These diagnostics are collected in
plain task-local memory and assembled into R objects after parallel regions.

## 12. Resolved input design

### 12.1 Envelope

The internal R specification is a uniquely named, versioned list:

```text
blr_resolved_spec
├── schema
├── data
├── model
├── prior
├── mcmc
├── compute
└── output
```

The resolver validates all cross-field invariants before native dispatch. The
C++ boundary receives only resolved values, not aliases or unvalidated `...`.

| Group | Required | Optional/defaulted | Validation owner |
|---|---|---|---|
| `data` | Analysis mode, traits, markers, operator resources, providers, likelihood regime | Statistical regions absent; no implicit overlap | R provider/marker resolver |
| `model` | Family, state, probability, scale, covariance, residual, and storage policies | Model-registered defaults only | R model registry plus cross-policy validator |
| `prior` | Every prior needed by selected policies | Fixed-state values only when their update is disabled | R statistical validator; native defensive checks |
| `mcmc` | Burn-in, sampling iterations, thinning, chains, seed | Task seeds absent; update flags use model-contract defaults | Shared R MCMC resolver |
| `compute` | Execution mode and cores | Memory limit and numerical operator controls | Shared R compute resolver |
| `output` | Retention and summary policy | Default compact summaries and core convergence quantities | Shared R output/memory resolver |

The resolved specification contains no missing scientific default: optional
means that the resolver supplies and records an approved value or that the
selected policy declares the field inapplicable.

### 12.2 Field groups

#### `schema`

| Field | Requirement |
|---|---|
| `name` | Exactly `blr_resolved_spec` |
| `version` | Positive supported integer |
| `compatibility_id` | Resolver/API migration identifier |

#### `data`

| Field | Requirement |
|---|---|
| `analysis_mode` | One of the three declared modes |
| `trait_ids` | Unique, nonempty, ordered character vector |
| `global_markers` | Unique IDs plus allele/coding metadata |
| `operator_resources` | Nonempty immutable BED/dense/CSR/eigen descriptors with unique IDs |
| `providers` | Nonempty validated provider list |
| `provider_maps` | Integer local-to-global maps with explicit missingness |
| `likelihood_regime` | Common-sample, independent-summary, or overlap-aware |
| `statistical_regions` | Optional declared marker-to-region map, separate from blocks |

Provider-specific sample sizes belong inside `providers`; a single top-level
sample size is permitted only as a validated common-sample reduction.

#### `model`

| Field | Requirement |
|---|---|
| `family` | Explicit BayesC, BayesR, BayesRC, or future registered family |
| `state_space` | Named states/patterns with exact ordering |
| `effect_storage_convention` | `base_latent`, `scaled_latent`, or `realised`, with realised output defined |
| `probability_policy` | Fixed/global/trait-specific/annotation/pattern policy descriptor |
| `marker_scale_policy` | Unit, component, MAF-S, annotation, or external $q_j$ descriptor |
| `marker_covariance_policy` | Scalar, global matrix, regional, or template policy |
| `residual_policy` | Scalar, diagonal, full common-sample, or overlap-aware |
| `update_order_version` | Registered versioned model transition identifier |

#### `prior`

The prior namespace contains mathematically named objects, including:

- scalar marker-variance prior parameters;
- $V_b$ inverse-Wishart degrees of freedom and scale matrix;
- residual variance or $V_e$ prior parameters;
- inclusion, component, or activity-pattern prior parameters;
- annotation-coefficient priors;
- marker multipliers $q_j$ and component multipliers $\gamma_k$;
- fixed versus sampled flags.

The public inverse-Wishart contract must name the parameterization. For

$$
V\sim\operatorname{IW}_T(\nu,\Psi),
$$

with density proportional to

$$
|V|^{-(\nu+T+1)/2}
\exp\left\{-\frac{1}{2}\operatorname{tr}(\Psi V^{-1})\right\},
$$

store `degrees_of_freedom = nu` and `scale = Psi`. Do not accept an ambiguous
"prior covariance" whose conversion depends on the backend.

#### `mcmc`

Use one canonical meaning for each control:

| Canonical field | Definition | Current alias |
|---|---|---|
| `burn_in_iterations` | Number of initial transitions discarded | `nburn` |
| `sampling_iterations` | Number of post-burn transitions executed | `nit` |
| `thin_interval` | Retain post-burn transitions divisible by this interval | `nthin` |
| `retained_draws` | Derived, never independently conflicting | current `nsamples`-like outputs |
| `chains` | Number of logical chains | `nchains` |
| `seed` | User seed before task derivation | `seed` |
| `task_seeds` | Optional resolved overrides | `chain_seeds` |
| `update_flags` | Named policy-specific update switches | `updateB`, `updateE`, `updatePi`, others |

Let post-burn transition indices be $u=1,\ldots,n_{\mathrm{sampling}}$. Retain
transition $u$ exactly when

$$
u\bmod n_{\mathrm{thin}}=0.
$$

Therefore,

$$
n_{\mathrm{retained}}
=
\left\lfloor
\frac{n_{\mathrm{sampling}}}{n_{\mathrm{thin}}}
\right\rfloor.
$$

The resolver requires $n_{\mathrm{burn}}\geq0$,
$n_{\mathrm{sampling}}\geq1$, and $n_{\mathrm{thin}}\geq1$, all finite scalar
integers. When retained draws are requested, it also requires
$n_{\mathrm{retained}}\geq1$. It records the exact retained transition
indices. `retained_draws` is derived and verified against allocated output.
Neither `niters` nor `nsamples` remains an alternative input meaning.
Migration adapters must source-trace each current wrapper's meaning of `nit`;
a wrapper whose historical meaning differs must convert explicitly or fail
rather than silently reinterpret it.

#### `compute`

| Field | Requirement |
|---|---|
| `execution_mode` | `serial` or `parallel` |
| `parallelization` | `none`, `chains`, `traits`, or `trait_chains`, validated against analysis mode |
| `cores` | Positive requested worker count |
| `scheduler_version` | Explicit deterministic scheduling-contract integer; seed version belongs in `schema` |
| `memory_limit_bytes` | `NULL`, a finite nonnegative scalar, or positive `Inf` |
| `operator_numerical_controls` | Uniquely named representation-specific numerical controls only |

Scientific controls must not be hidden in `compute`. For example, retained rank
changes the approximate likelihood and therefore also appears in provider
provenance and the resolved model record.

#### `output`

| Field | Requirement |
|---|---|
| `posterior_summaries` | Logical scalar selecting the contracted posterior-summary set |
| `retained_parameters` | Named scientific parameters to retain |
| `effect_draw_policy` | `none`, selected markers, compact, or full |
| `state_draw_policy` | `none`, selected markers, compact, or full |
| `convergence_policy` | Registered mode and selected-quantity policy |
| `derived_quantities` | Explicit requested derived statistics |
| `preserve_chains` | Whether chain identity is retained; raw v2 always preserves chain axes |

Defaults belong to the resolver and are recorded after resolution. Native code
does not invent omitted defaults.

## 13. Model-specific input extensions

Shared infrastructure does not imply one generic `pi` argument. Model policies
must declare their own inputs:

| Concept | Explicit internal name | Example scope |
|---|---|---|
| Scalar inclusion probability | `inclusion_probability` | ST BayesC |
| Traitwise inclusion probabilities | `trait_inclusion_probabilities` | Independent ST traits |
| Activity-pattern definition | `activity_patterns` | Cheng MT-BayesC$\Pi$ |
| Activity-pattern mass | `activity_pattern_probabilities` | Cheng MT models |
| Component multipliers | `component_multipliers` | BayesR |
| Component prior mass | `component_prior_mass` | ST or declared MT component model |
| Marker multipliers | `marker_variance_multipliers` | $q_j$, MAF-S, annotation `Q` |
| Probability architecture | `marker_probability_policy` | Annotation `P` |
| Marker covariance prior | `marker_covariance_prior` | $V_b$ policy |
| Residual covariance prior | `residual_covariance_prior` | $V_e$ policy |
| Region assignment | `statistical_region_id` | Future regional covariance |
| Covariance templates | `covariance_templates` and `template_weights` | Future mash-like policy |

PIPs are posterior summaries and are never an input probability. Prior state,
component, and activity-pattern masses retain their separate names.

## 14. Raw schema version 2

### 14.1 Envelope

```text
blr_raw
├── schema
├── model
├── input
├── posterior
├── draws
├── final
├── derived
├── diagnostics
└── provenance
```

The stable R raw object should use one class, `blr_raw`, for all analysis modes.
Analysis mode and model family are fields, not separate schema classes.

### 14.2 Namespace contract

#### `schema`

- `name = "blr_raw"`;
- `version = 2L`;
- `compatibility_id` identifying formatter/migration support;
- `dimension_contract_version`.

#### `model`

- `analysis_mode`;
- `family` and data-level model name;
- ordered state-space definition and null-state index;
- probability, scale, marker-covariance, and residual policies;
- effect-storage convention;
- update-order/transition version;
- approximation-sensitive scientific choices.

#### `input`

- resolved trait and global marker IDs;
- operator-resource descriptors, provider summaries, trait sets, resource
  references, and local maps;
- validated priors;
- resolved MCMC and compute controls;
- requested output policy;
- initial scientific states and whether supplied or resolved;
- migration/deprecation actions applied by wrappers.

Large immutable matrices need not be duplicated in the result. Store stable
identifiers, dimensions, hashes where appropriate, and sufficient provenance
to reconstruct the contract.

#### `posterior`

- `realised_effect_mean` for realised marker effects;
- `latent_effect_mean` only when scientifically requested and interpretable;
- `pips`;
- explicitly named traitwise or joint marker-state probabilities;
- `activity_pattern_probabilities` where applicable;
- explicitly separate traitwise or joint component-assignment probabilities;
- covariance/variance posterior summaries derived from draws;
- uncertainty summaries with estimator and chain aggregation metadata.

Posterior means are not final states, and posterior component probabilities are
not prior component masses.

#### `draws`

- sampled marker-variance or marker-covariance parameters;
- sampled residual variance/covariance parameters;
- sampled probability parameters in policy-specific arrays with declared axes;
- optional latent, realised-effect, and state draws under the output policy;
- chain-preserving convergence draws;
- draw iteration indices.

An actual sampled $V_b$ draw must remain an actual sampled draw. A descriptive
effect covariance belongs in `derived` with a different name.

#### `final`

For every logical chain:

- final realised effects;
- final latent/base effects if part of the chain state;
- final marker states, activity patterns, and components;
- final probability parameters;
- final marker and residual variance/covariance states;
- final operator residual state only if continuation is intentionally
  supported.

RNG continuation state is absent by default. If later supported, its engine and
serialization version must be explicit.

#### `derived`

- predictions and prediction provenance;
- identified genetic variances and genomic covariances;
- heritability and polygenicity with definitions;
- component-specific genetic contributions;
- operator-relative residual/genomic quadratics for indefinite or approximate
  operators;
- descriptive cross-operator bilinear forms with non-covariance names;
- CPO or other model diagnostics where supported.

Only a genuine genotype cross-product or compatible positive-semidefinite
covariance operator supports unqualified variance/covariance terminology.

#### `diagnostics`

- convergence results and retained convergence traces;
- acceptance/proposal diagnostics;
- ESS/MCSE/R-hat where defined;
- runtime and memory;
- requested/configured/actual workers and worker IDs;
- numerical safeguards, repairs, failures, and counts;
- operator approximation and indefiniteness warnings;
- stable present-but-`NULL` fields for inapplicable contracted diagnostics.

#### `provenance`

- package version;
- source Git SHA when available;
- dirty-development-build indicator when available;
- compiler, architecture, OpenMP, BLAS, and build identifiers;
- operator and provider provenance;
- marker/allele alignment decisions;
- block/eigen retained-rank details;
- approximation status;
- seed resolution and task identity table;
- optional analysis timestamp.

### 14.3 Git provenance

Do not invoke Git during every fit. At source-build or development-load time,
resolve a build provenance record containing package version, SHA, and dirty
indicator, then expose it through an installed package constant or compiled
metadata. A release source archive without Git metadata records the package
version and archive/build identifier with `git_sha = NULL`; it must not invent a
SHA. A development build may record `dirty = TRUE` without enumerating user
files in fit objects.

## 15. Dimension conventions

### 15.1 R-facing order

Raw schema version 2 preserves all scientific axes:

| Quantity | Canonical dimensions |
|---|---|
| Retained realised effects | `draw × chain × marker × trait` |
| Retained latent effects | `draw × chain × marker × trait` |
| Independent-trait marker states | `draw × chain × marker × trait` |
| Joint marker-state indices | `draw × chain × marker` |
| Traitwise binary activity | `draw × chain × marker × trait` |
| PIPs | `marker × trait` |
| Traitwise posterior state probabilities | `marker × trait × state` |
| Joint posterior state probabilities | `marker × joint_state` |
| Activity-pattern posterior probabilities | `marker × activity_pattern` |
| Traitwise state-parameter draws | `draw × chain × trait × state` |
| Joint state-parameter draws | `draw × chain × joint_state` |
| Activity-pattern parameter draws | `draw × chain × activity_pattern` |
| Traitwise component-assignment probabilities | `marker × trait × component` |
| Joint component-assignment probabilities | `marker × component` |
| $V_b$ and $V_e$ draws | `draw × chain × trait_row × trait_col` |
| Scalar trait variances | `draw × chain × trait` |
| Final effects | `chain × marker × trait` |
| Final independent-trait states | `chain × marker × trait` |
| Final joint states | `chain × marker` |
| Final covariance states | `chain × trait_row × trait_col` |
| Predictions | `draw × chain × observation × trait`, or a documented summary without draw axis |
| Regional covariance states | `draw × chain × region × trait_row × trait_col` |
| Provider diagnostics | `provider` plus named block/task axes as applicable |

Every axis has stable `dimnames` or an adjacent ID vector. Every retained array
also carries an exact named-axis descriptor validated against the registered
field contract. `T=1`, one chain, one state, and one draw retain their axes in
raw output. A formatter may simplify explicitly and record the convenience
representation. Policy-specific arrays use distinct names and registered fixed
axes; one field must not switch rank or axis meaning between models.

### 15.2 Native memory layout

Native code may use marker-major, trait-major, or task-major contiguous buffers
for performance. `RawResultBuilder` owns the one explicit conversion into the
R-facing order. Native layout is recorded in code contracts but is not exposed
as raw-schema semantics.

Avoid matrices whose dimensions change meaning by model. If an axis is not
applicable, use `NULL` under the named field or an explicit zero-length axis as
the field contract requires; do not repurpose another axis.

## 16. Naming semantics

| Name | Reserved meaning |
|---|---|
| `state_probabilities` | Retired as an unqualified variable-rank field; use an explicit traitwise, joint-state, activity-pattern, or component name in the appropriate namespace |
| `traitwise_state_probabilities` | Markerwise posterior mass with explicit `marker × trait × state` structure |
| `joint_state_probabilities` | Mass over one declared joint marker-state space |
| `pattern_probabilities` | Retired as an unqualified field; use `activity_pattern_probabilities` with a prior/posterior/draw namespace |
| `prior_state_mass` | Resolved prior mass over the declared complete marker state space |
| `posterior_state_probabilities` | Conceptual term only; raw fields use the explicit traitwise or joint-state name and fixed axes |
| `activity_pattern_probabilities` | Probability parameters over trait-activity patterns |
| `traitwise_component_assignment_probabilities` | Marker-by-trait posterior mass over BayesR-type scale components |
| `joint_component_assignment_probabilities` | Markerwise posterior mass over one shared/joint component state |
| `pips` | Marker-trait posterior inclusion probabilities |
| `joint_non_null_probabilities` | Posterior probability that a marker is non-null in any declared joint sense, with definition metadata |
| `pleiotropic_probabilities` | Posterior probability of a declared multi-trait active pattern or set of patterns |
| `marker_variance` | Scalar latent/base marker-effect variance |
| `marker_covariance` | Latent/base marker-effect covariance $V_b$ |
| `residual_variance` | Scalar observation-level residual variance |
| `residual_covariance` | Observation-level residual covariance $V_e$ under an identified likelihood |
| `genetic_variance` | Variance of genetic values under a valid PSD genotype operator |
| `genomic_covariance` | Covariance of genetic values under a valid common/target-population PSD operator |
| `operator_relative_genomic_quadratic` | Algebraic quadratic under an operator not established as PSD |
| `operator_relative_residual_quadratic` | Algebraic residual quadratic under such an operator |
| `descriptive_cross_operator_bilinear` | Cross-trait descriptive summary not identified as covariance |
| `prior_probabilities` versus `posterior_probabilities` | Prior parameter/state mass and posterior marker/state summaries must occupy different namespaces and never share one field |

`pis` is retired from raw schema version 2. Current `fit$pis`, `pi_trace`,
`pi_final`, `pi_mean`, `pi`, and `pim` have model-dependent meanings. Migration
must map each source field using model metadata; no single schema-v2 field may
mean prior state mass in one model and PIP in another.

## 17. Raw versus formatted objects

### 17.1 Native result

The native result consists of validated typed buffers and dimension metadata.
It contains scientific states, retained draws, sufficient posterior
accumulators, and task-local diagnostics. It contains no R presentation aliases
and constructs no R objects in worker regions.

### 17.2 Stable R raw object

`blr_raw` version 2 is the canonical language-level scientific object. It:

- preserves all axes and chains;
- distinguishes posterior summaries, draws, and final states;
- uses explicit model-specific probability names;
- stores sampled parameters without replacement;
- contains enough provenance for interpretation;
- fails strict schema validation before formatting.

### 17.3 Formatted fit

The formatted fit owns:

- optional dimension simplification;
- user-facing tables and summaries;
- plotting inputs;
- documented compatibility aliases;
- convenience correlations derived from covariance draws/summaries;
- model-specific convenience fields;
- `check_stblr_consistency()`-style validation.

It must not overwrite a sampled covariance with an empirical covariance,
manufacture unavailable fields, or collapse final/draw/posterior meanings.
Present-but-`NULL` stable fields are inserted with list semantics so the names
remain present.

## 18. Migration from current interfaces and outputs

### 18.1 Common current R arguments

| Current argument/concept | Resolved specification | Migration action |
|---|---|---|
| `nit` | `mcmc$sampling_iterations` | Retain wrapper alias initially; document exact meaning |
| `nburn` | `mcmc$burn_in_iterations` | Retain wrapper alias initially |
| `nthin` | `mcmc$thin_interval` | Retain wrapper alias initially |
| `nchains` | `mcmc$chains` | Retain wrapper alias initially |
| `ncores` | `compute$cores` | Retain wrapper alias initially |
| `seed`, `chain_seeds` | `mcmc$seed`, resolved `mcmc$task_seeds` | Version task derivation |
| `method` | `model$family` plus data-level/operator semantics | Reject ambiguous aliases |
| `stats`, `Glist`, `ld_prefix` | `data$providers` | Resolve local maps and provenance once |
| `y`, BED `Glist` | common-sample provider collection | Separate phenotype and genotype provenance |
| `pi` | policy-specific probability field | Remove generic internal meaning |
| `models`, `pimodels` | `model$state_space`, activity-pattern policy | Normalize ordering and null index |
| `mixture_var` | `prior$component_multipliers` | Use $\gamma_k$ semantics explicitly |
| `joint_pi` | declared activity/component state mass | Split according to model geometry |
| `h2`, `vg`, `vb`, `ve` | named initialization/prior calibration fields | Distinguish initial values from priors |
| `ssb_prior`, `sse_prior`, `nub`, `nue` | explicit prior distributions | Migrate with parameterization checks |
| `sets` | computational partition or `statistical_region_id` | Require one declared role |
| `representation`, eigen controls | provider operator descriptor | Record exact/approximate status |
| annotation arguments | `model$probability_policy`/`marker_scale_policy` | Resolve `P` and `Q` separately |
| `keep_chains`, convergence controls | `output` namespace | One retention/convergence contract |

Current flat wrappers should resolve into the structured form before a new
structured public API is considered. Once migration is complete, a public
`blr_fit(specification)`-style entry point may be added without removing
specialized convenience wrappers.

#### ST wrapper mapping

| Current ST input | Resolved destination | Decision |
|---|---|---|
| `stblr_csr(stats, Glist, ld_prefix, ...)` | Summary providers with CSR operators and local maps | Preserve wrapper; resolve all providers before dispatch |
| `stblr_bed(y, Glist, ...)` | Individual-level provider and trait collection | Preserve wrapper; use shared BED view |
| `stblr_block_eigen(stats, Glist, block_start, ...)` | Summary provider with full/retained eigen operator | Move block/eigen controls into provider descriptor |
| `stblr_csr_annot(..., annotation_model, annotations)` | Base model plus explicit probability/scale policy | Preserve provider-specific scientific validation |
| ST `method` | Model family plus likelihood/data-level name | Stop deriving statistical meaning from an `s` prefix internally |
| ST `effect_maf`, `maf_effect_s` | Marker scale policy and provenance | Preserve scale and population alignment |
| ST probability/variance annotation controls | Named `P` or `Q` policy inputs | Do not merge into generic annotations |

#### MT wrapper mapping

| Current MT input | Resolved destination | Decision |
|---|---|---|
| `mtblr_bed(y, Glist, residual_covariance, ...)` | Common-sample individual provider plus residual policy | First new MT vertical slice |
| `mtblr_csr(stats, Glist/ld_prefix, trait_metadata, ...)` | One or more summary providers per trait | Generalize beyond shared LD while preserving explicit no-overlap regime |
| `mtblr_block_eigen(..., operator_sharing, block_start)` | Provider-specific block-eigen descriptors | Replace shared/reconstructed assumptions with local owners |
| `sample_overlap` | `data$likelihood_regime` and overlap contract | Current `not_modeled` maps to independent-summary regime |
| `models`, `pimodels` | Ordered activity-pattern state space | Preserve explicit null state and names |
| `joint_pi`, `joint_pi_prior` | Activity-pattern or joint pattern-component policy | Split pattern and component semantics explicitly |
| `beta`, `b`, `state`, `component` | Named initial latent, realised, pattern, and component states | Reject storage ambiguity |
| `B`, `vb`, `ssb_prior`, `nub` | Initial marker covariance and explicit inverse-Wishart prior | Remove current heuristic/IW scale collision |
| `E`, `ve`, `sse_prior`, `nue` | Initial residual state and likelihood-compatible prior | Separate full common-sample from diagonal summary policy |
| `sets` | Declared computational partition or statistical region | Never infer one role from the other |

### 18.2 Current native/raw fields

| Current field | Schema-v2 destination | Status |
|---|---|---|
| `bm` | `posterior$realised_effect_mean` | Retain meaning, rename explicitly |
| `dm` | `posterior$pips` when verified | Retain only with model-specific definition |
| `b`, `beta_final`, `d`, `component_final` | `final` named states | Preserve chain axis and storage convention |
| `vbs` | `draws$marker_variance` or `draws$marker_covariance` | Split scalar/matrix semantics |
| `ves` | `draws$residual_variance` or covariance diagonal summary | Split by residual policy |
| `cov_b_*` | sampled `draws/final$marker_covariance` or posterior summary | Correct current MT hybrid before migration |
| `cov_g_*` | `derived$genomic_covariance` only when identified | Otherwise rename descriptive bilinear |
| `vle`, `vld` | `derived` with exact/operator-relative qualification | Preserve algebraic decomposition metadata |
| `pis`, `pi_*`, `pi`, `pim` | policy-specific probability fields | Deprecate ambiguous names |
| `component_probabilities` | traitwise or joint component-assignment field after model-specific verification | Split rather than infer rank |
| `chains` | chain axis in `draws` and `final` | Replace nested shape variability |
| `convergence_traces` | `draws$convergence` plus metadata | Preserve chain and iteration axes |
| `diagnostics` | `diagnostics` | Stabilize names and present-but-`NULL` fields |
| `input`, `data` | resolved `input` and `provenance` | Separate scientific inputs from lineage |

### 18.3 Current formatted fields

| Current formatted field | New formatted source | Migration status |
|---|---|---|
| `bm` | `posterior$realised_effect_mean` | Retain temporary alias |
| `dm` | `posterior$pips` | Retain temporary alias after semantic validation |
| `b`, `d`, `component_final` | `final` chain states | Retain aliases only with explicit primary-chain policy |
| `vbs`, `vgs`, `ves`, `vle`, `vld` | Named `draws` or `derived` quantities | Deprecate compact ambiguous names |
| `cov_b_mean/final` | Posterior/final sampled marker covariance | Retain only after MT covariance correction |
| `cov_g_mean/final` | Identified derived genomic covariance | Rename when only a descriptive cross-operator bilinear is available |
| `cov_e_mean/final` | Posterior/final residual covariance under the declared likelihood | Use `NULL` off-diagonals when not identified |
| `pi_trace/final/mean`, `pi`, `pim` | Explicit policy-specific probability fields | Deprecate generic aliases |
| `component_probabilities` | Posterior component probabilities | Retain clearer alias if unambiguous |
| `chains` | Chain-preserving `draws` and `final` | Replace nested model-varying objects |
| `diagnostics`, `convergence`, `convergence_traces` | Same schema-v2 namespaces | Retain names; normalize dimensions and metadata |
| `input`, `data`, `model_parameters` | Resolved `input`, `model`, and `provenance` | Replace duplication with explicit namespaces |

The first formatter may provide deprecated aliases such as `bm`, `dm`, `vbs`,
and `pi_trace`, but each alias must be a read-only view of one unambiguous
schema-v2 quantity. `fit$pis` should not be manufactured as a universal field.
If retained temporarily, it must have one documented meaning or be `NULL` for
models where that meaning does not exist.

Serialized version-1 objects require an explicit reader/converter or a clear
unsupported-version error. Silent positional fallback remains prohibited.

### 18.4 Generated interfaces

Future native entry points will change `R/RcppExports.R` and
`src/RcppExports.cpp` as generated consequences. The current
`mtblr_csr_chains_raw_internal()`, `mtblr_block_eigen_chains_raw_internal()`,
`mtblr_bed_chains_internal()`, and corresponding scalar scheduled entry points
are expected to become adapters to a small number of versioned
resolved-spec/result boundaries. Old wrappers remain until their vertical slice
is qualified. Generated files are not changed at the design checkpoint.

## 19. Reductions and validation strategy

### 19.1 Independent-trait execution

For identical logical task identities and seed-contract version,

$$
\text{serial independent traits}
=
\text{parallel independent traits}
$$

at the task-level scientific output and retained-draw level. Tests must include
several traits and chains so that more than one worker can actually be used.

### 19.2 Backend representation

Where preprocessing and likelihoods are exactly equivalent,

$$
\text{BED}
=
\text{dense cross-product}
=
\text{CSR}
=
\text{full-rank block eigen}
$$

within predeclared floating-point tolerances and, where RNG draw order is
contracted, exact logical-chain equality. Retained-rank block eigen must equal
its reconstructed retained operator, not the original dense cross-product.
Full-rank block eigen equals the dense cross-product only when the declared
dense operator is block diagonal under the same partition; otherwise it equals
the block-diagonal projection and must be labelled accordingly.

### 19.3 Ordering invariance

Require invariance to:

- provider ordering;
- consistent provider-local marker permutation and map permutation;
- eigenblock traversal order;
- computational set order when sets have no statistical role.

Statistical region order may alter storage order but not the posterior when IDs,
priors, seeds, and update schedule are mapped consistently.

### 19.4 Joint-MT reduction to independent traits

Diagonal $V_b$ alone is insufficient. A joint model reduces to independent
trait posteriors only when all coupling is removed, including:

1. the likelihood factorizes by trait/provider, including diagonal or
   independent residual/error structure;
2. marker-effect covariance and all other continuous priors factorize;
3. activity-state probabilities factorize into traitwise states;
4. there are no shared learned probability, scale, annotation, or covariance
   parameters;
5. priors on all traitwise parameters factorize;
6. update schedules implement compatible traitwise conditionals;
7. logical seeds are mapped to the same trait-chain identities;
8. requested derived quantities do not feed back into the chain.

Under those conditions, compare joint-policy output with independent ST fits.
Otherwise, do not advertise an independence reduction.

### 19.5 Scientific and schema oracles

Mandatory implementation tests include:

- analytical BayesC/BayesR marker conditionals;
- the independent Cheng MT research references for $T=2$;
- null collapse and conditional completion;
- one authoritative sampled $V_b$ without heuristic replacement;
- provider dense/CSR/eigen and ordering reductions;
- exact-zero and no-op RNG reductions for optional policies;
- serial/parallel logical-task equality;
- stable schema names and dimensions;
- no dropped `T=1`, one-chain, or one-draw axes in raw objects;
- strict final/draw/posterior separation;
- prior/posterior probability naming;
- build and operator provenance;
- rejection of unnamed, duplicate, partial, and unknown arguments;
- raw-to-fit consistency without scientific reinterpretation.

Permanent expected values must come from analytical identities, independent
reference code, or validated external parameterization crosswalks—not from the
same production implementation on both sides.

## 20. Phased implementation roadmap

### Phase 0: approve unified contracts

| Item | Detail |
|---|---|
| Likely files | This design, `docs/dev/blr_architecture.md`, output/model/RNG contracts |
| Tests | Static contract crosswalks; proposed schema fixtures independent of native code |
| Scientific acceptance | Analysis/execution modes, resource/provider contract, dimensions, naming, $V_b$ and staged $V_e$ policies approved |
| Migration risk | Ambiguous current probability and iteration names |
| Checkpoint | Contract-only commit before executable changes |

Approve the resolved specification, provider/operator contract, logical-task
RNG contract, raw schema version 2, and naming/dimension conventions before
changing a kernel.

### Phase 1: shared R specification and result envelope

| Item | Detail |
|---|---|
| Likely files | New R resolver/validator/schema modules; existing public wrappers; `R/blr-unified.R`; formatter and consistency helpers; schema tests |
| Tests | Exact-name validation, alias migration, schema dimensions, v1/v2 failure/conversion, formatted aliases |
| Scientific acceptance | Existing kernels produce unchanged scientific values through the new envelope |
| Migration risk | `...` behavior, overloaded `pi`, dropped dimensions |
| Rollback boundary | Wrappers can retain the old execution adapter behind one resolver switch |

No statistical kernel changes in this phase.

### Phase 2: shared provider and operator adapters

| Item | Detail |
|---|---|
| Likely files | R provider/resource preparation; `src/st_ld_operator.h`; BED/block-eigen preparation; new immutable resource/provider/map headers; operator tests |
| Tests | Shared-resource reuse; multi-trait BED provider; BED/dense/CSR/full-rank eigen reductions; retained-rank reconstruction; maps, alleles, missing markers, provider order |
| Scientific acceptance | Existing ST likelihoods and trajectories preserved under identical task seeds |
| Migration risk | Cross-product scale, allele orientation, residual-state differences |
| Rollback boundary | One operator representation at a time behind the common contract |

### Phase 3: shared scheduler, RNG, retention, and diagnostics

| Item | Detail |
|---|---|
| Likely files | `R/blr-unified.R`; scheduled execution headers; convergence modules; native result builder; reproducibility tests |
| Tests | Trait×chain task seeds; serial/parallel equality; zero-RNG diagnostics; burn/thin/draw indices; actual worker use |
| Scientific acceptance | Execution mode leaves target and retained draws unchanged |
| Migration risk | RNG ordering and chain aggregation |
| Rollback boundary | Versioned scheduler/seed contract selectable internally during qualification |

### Phase 4: first production MTBLR vertical slice

| Item | Detail |
|---|---|
| Likely files | New MT state/covariance policy headers and core; MT BED adapter; MT public resolver; schema-v2 formatter; focused tests |
| Tests | Common-sample $T=2$ Cheng reference, joint categorical states, null collapse, conditional completion, fixed-full-$V_e$ reference, then sampled-full-$V_e$ inverse-Wishart and diagonal reductions |
| Scientific acceptance | Matches independent research reference; one sampled $V_b$ is used, retained, summarized, and returned |
| Migration risk | Existing MT results change; current hybrid outputs are not backward compatible scientifically |
| Rollback boundary | New method/version identifier; old MT route remains clearly deprecated until qualification completes |

The first slice is common-sample BED MT-BayesC$\Pi$. Its first scientific
sub-checkpoint fixes a supplied symmetric positive-definite full $V_e$, which
isolates and validates the Cheng marker-state and $V_b$ transition. Before the
vertical slice is promoted as the maintained production replacement, a second
sub-checkpoint adds a sampled full inverse-Wishart $V_e$ under the explicitly
declared parameterization; a diagonal policy is retained as a separate
reduction and for likelihood regimes that do not identify off-diagonals. It
does not reuse the obsolete MT covariance hybrid.

### Phase 5: heterogeneous summary-statistic MTBLR

| Item | Detail |
|---|---|
| Likely files | Summary provider resolvers; CSR provider adapter; MT summary likelihood policy; no-overlap contracts/tests |
| Tests | Different $C_{dt}$, $N_{dt}$, maps, populations, multiple providers; provider-order invariance; dense/CSR equality |
| Scientific acceptance | Trait/provider likelihood contributions combine without a common materialized operator; residual policy is identifiable |
| Migration risk | Existing descriptive off-diagonal `cov_g` fields and assumed shared LD |
| Rollback boundary | Independent-summary regime gated separately from common-sample MT |

### Phase 6: provider-specific block eigen

| Item | Detail |
|---|---|
| Likely files | Block-eigen provider preparation and views; retained residual policy; provenance/result modules; reduction tests |
| Tests | Provider-specific blocks/eigensystems/ranks; full-rank reduction; retained-operator equality; traversal invariance |
| Scientific acceptance | Each provider evaluates its declared operator without forcing common blocks or eigenvectors |
| Migration risk | Confusing retained approximation with dense truth; block residual policies |
| Rollback boundary | Full-rank promotion before retained-rank promotion |

### Later gated phases

| Phase | Gate |
|---|---|
| Shared-component MT-BayesR | BayesC MT slice accepted; $\gamma_kq_j$ storage/statistic contract fixed |
| MT-BayesRC | Probability architecture and sharing policy independently validated |
| Regional covariance | Persistent $V_{b,r}$ prior/regularization selected; update once per iteration; low-information tests pass |
| Covariance templates | Template library and singular active-subspace semantics approved |
| Overlap-aware summaries | Cross-provider error information and likelihood identified |
| Larger-$T$ state models | Restricted patterns/factors have independent target derivations |

Each phase ends with a manual checkpoint containing its contracts, production
changes, independent tests, and migration note. Rollback means returning to the
prior checkpoint, not silently mixing old and new schema or posterior logic.

## 21. Decision register

| Decision | Options | Evidence | Recommendation | Consequence | Blocks |
|---|---|---|---|---|---|
| Public structured API timing | Immediate new public API; internal spec first | Existing wrappers are numerous and stable; current resolution is duplicated | Build the internal structured spec first; keep exact-name wrappers; consider public `blr_fit()` after Phase 2 | Lowers migration risk without preserving unsafe forwarding | No Phase 1; public API later |
| Analysis/execution modes | One overloaded mode; separate fields | Current ST schedules trait×chain while MT schedules chains | Separate enumerations exactly as in Section 3 | Prevents parallel traits being mislabeled as MT | Yes, Phase 1 |
| Execution-policy encoding | Open set of mode strings; closed mode/axis fields | Trait parallelism is invalid for a joint chain; combined trait-chain pools are useful for independent traits | `execution_mode = serial/parallel` plus closed `parallelization = none/chains/traits/trait_chains`, validated by analysis mode | Removes ambiguous compatible sets | Yes, Phase 1 |
| Seed derivation | Backend-specific arithmetic; unspecified hash; versioned stable mixing | Reproducibility must survive worker and backend changes | Seed-contract v1 uses FNV-1a for stable trait IDs and ordered SplitMix64 mixing; explicit task-seed arrays have contracted shapes | Requires frozen reference vectors | Yes, Phase 3 contract; specification fixed in Phase 0 |
| Retention indexing | Retain first post-burn draw; retain every divisible post-burn index; backend-specific | Current kernels have repeated local rules | Index post-burn transitions from one and retain $u\bmod n_{\mathrm{thin}}=0$; source-trace wrapper aliases during migration | Fixes draw counts and iteration indices | Yes, Phase 1 |
| R-facing dimension order | Task-major; variable-major; draw-first | R convergence use and current arrays favor explicit draw/chain axes | Draw, chain, scientific object axes; marker before trait | Native builder performs one conversion | Yes, Phase 1 |
| Raw schema versioning | Extend v1 in place; v2 envelope | ST/MT v1 schemas differ and names are overloaded | New `blr_raw` v2 with strict compatibility ID | Requires explicit v1 handling | Yes, Phase 1 |
| Git provenance | Runtime Git; build-time record; omit | Runtime Git is unreliable for installed packages | Once-per-session load/development cache; unavailable `git_sha` and `dirty_build` remain present with `NULL`; future build-time injection is a promotion task | Reproducible without fit-time subprocesses | Phase 1 cache implemented; build-time injection remains |
| Initial MT $V_e$ policy | Fixed full; sampled full; sampled diagonal; staged combination | Common-sample identifies full; fixing $V_e$ isolates the marker/$V_b$ oracle; independent summaries do not identify off-diagonals | Phase 4a uses supplied fixed SPD full $V_e$; Phase 4b adds sampled full inverse-Wishart $V_e$ before promotion; diagonal is a declared reduction and the default class for no-overlap summaries | Gives an exact first oracle without leaving residual uncertainty out of the maintained model | Yes, MT slice |
| Inverse-Wishart public prior | Mean/df shortcut; explicit $(\nu,\Psi)$; both | Current scale collision is documented | Canonical explicit df and scale; optional R calibration helper with recorded conversion | Removes backend ambiguity | Yes, MT slice |
| Effect draw retention | Always full; never; policy-based | Full marker draws are costly but selected convergence is useful | Default summaries plus selected/compact draws; full draws opt-in with memory guard | Raw dimensions remain stable when present | No Phase 1 |
| State-probability naming | Generic `pi`; explicit model names | Current `pis/pi` meanings conflict | Explicit prior/posterior state, pattern, and component names | Requires migration aliases | Yes, Phase 1 |
| `fit$pis` | Universal field; model-dependent alias; retire | Audit confirms no consistent universal meaning | Retire as raw field; formatted alias only if one documented meaning applies | Prevents PIP/prior confusion | Yes, schema contract |
| Compatibility aliases | Preserve all; none; targeted | Package is developmental, but specialized wrappers are useful | Targeted, versioned read-only aliases with warnings; no positional fallback | Intentional migration without semantic ambiguity | No Phase 1 |
| C++ policy boundary | All virtual; all templates; hybrid | Hot scalar templates exist; representations are runtime choices | Runtime resolved variants, compile-time hot policies, plain state structs | Reuse without class hierarchy overhead | Yes, operator/core design |
| Initial MT covariance transition | Current hybrid; full latent; null-collapsed completion | Research design validates full latent and completed sampler; hybrid is incoherent | Cheng full latent as oracle; joint categorical + null collapse + conditional completion in production | Existing MT covariance results change | Yes, MT slice |
| Provider/operator ownership | One operator embedded per provider; separate reusable resources | Common-sample MT has a joint multi-trait likelihood, and several providers may share BED/LD storage | Immutable operator resources referenced by one or more likelihood providers; providers own nonempty trait sets | Supports full $V_e$ without factorization and avoids duplicated decoding/storage | Yes, Phase 2 |
| Regional covariance | Independent unrestricted matrices now; shared/template later; defer | Low-information regularization unresolved | Defer production; retain a gated policy extension point | Avoids encoding unstable regional default | No core phases |

## 22. Acceptance criteria for production implementation

Production work may begin only when Phase 0 approves:

1. analysis and execution modes as separate concepts;
2. the resolved specification and exact argument-validation rules;
3. global marker, reusable operator-resource, and multi-trait-capable provider
   contracts;
4. provider-specific dense, BED, CSR, and block-eigen operator semantics,
   including the exact block-diagonal boundary;
5. the versioned logical-task RNG/scheduler contract;
6. schema-v2 namespaces, fixed named axes, dimensions, and probability names;
7. raw versus formatted responsibilities;
8. explicit inverse-Wishart and staged fixed-full/sampled-full $V_e$ policies
   for the first MT slice;
9. migration behavior for current wrappers and schema-v1 objects;
10. independent acceptance tests and checkpoint boundaries.

The first implementation must preserve current supported ST scientific results
while moving shared infrastructure behind the new contracts. The first MT
vertical slice must match the independent Cheng reference and must not reuse the
current hybrid covariance transition. Heterogeneous $X_t^\top X_t$ and
provider-specific block eigen are first-class requirements of the architecture,
even though they arrive after the common-sample MT slice.

No production phase is accepted if it:

- builds a separate MT provider/operator/result stack;
- conflates independent traits with joint MT;
- lets execution mode change the posterior;
- lets blocks or traversal order change global parameter update frequency;
- calls all probability quantities `pi`;
- drops raw chain, trait, marker, state, or covariance dimensions;
- replaces sampled covariance with a descriptive estimator;
- treats missing provider markers as zero effects;
- reports retained-rank operators as exact dense operators;
- silently accepts positional or partially matched scientific arguments.

## 23. Phase 0 current-behavior trace

This section freezes the migration evidence. It does not promote current
schema version 1 or the current MT covariance transition into the target.

### 23.1 Current iteration and retention semantics

Current public wrappers use `nit`, `nburn`, and `nthin`. No maintained wrapper
uses `niters`; `nsamples` is a derived native/output count, not an independent
MCMC input. Native transition index `it` is zero-based and runs from zero
through `nit + nburn - 1`. Final state is the state after the last transition.

| Route class | Executable evidence | Current marker-summary retention | Current traces/convergence | Migration |
|---|---|---|---|---|
| ST CSR BayesC/BayesR and scheduled variants | `src/blr_scalar_execution.h`, `scalar_iteration_is_retained()`; scalar core implementations | `it >= nburn && (it - nburn) % nthin == 0`; post-burn indices $1,1+n_{\mathrm{thin}},\ldots$ and count $\lceil n_{\mathrm{sampling}}/n_{\mathrm{thin}}\rceil$ | Scientific variance/probability traces have `nit + nburn` entries including burn-in; formal convergence capture has `nit` unthinned post-burn entries | Legacy retention-contract `st_scalar_v1`; adapter records exact legacy indices until deliberate migration |
| ST BED BayesC/BayesR/BayesRC | `src/blr_bed_scheduled_bayesc_core_impl.h`, `src/blr_bed_bayesr_core_impl.h`, `src/blr_bed_bayesrc_core_impl.h` | Same first-post-burn rule | Same full scientific traces and unthinned post-burn convergence capture | Legacy retention-contract `st_bed_v1` |
| ST annotation providers | learned/group/prior/log-variance/SBayesRC cores under `src/` | Same repeated first-post-burn rule in maintained routes | Provider-specific parameter traces may be full or retained; formal convergence remains unthinned post-burn | Preserve route-specific v1 descriptor; do not infer from array length |
| MT CSR/block-eigen default cores | `src/blr_mt_default_core_impl.h` | Marker and component summaries use the first-post-burn rule | Probability/covariance means use every post-burn transition; scientific traces include burn-in | Legacy retention-contract `mt_default_v1`; mixed aggregation is recorded, not normalized silently |
| MT BED chains | `src/blr_mt_bed_core_impl.h` and chain execution wrappers | Marker and component summaries use the first-post-burn rule | Covariance/probability means and convergence use every post-burn transition; displayed traces are chain means | Legacy retention-contract `mt_bed_v1` |

The route-level boundary trace is:

| Route | R entry | Native boundary/core owner | Current output shape | Cross-route consistency | Migration action |
|---|---|---|---|---|---|
| ST CSR | `stblr_csr()` | `stblr_cpg_omp_csr()`, `stblr_cpg_omp_csr_bayesr()`, scalar CSR cores | marker × trait summaries; scientific traces `nit + nburn` × trait; convergence `nit × chain × quantity` | Same legacy retention helper in maintained scalar cores | Resolve into spec while retaining `st_scalar_v1` |
| ST block eigen | `stblr_block_eigen()` | block-eigen adapters ending in scalar CSR/BayesR/SBayesRC cores | same scalar shapes plus block diagnostics | Retention follows the dispatched scalar kernel | Preserve kernel-specific legacy descriptor |
| ST BED | `stblr_bed()` | scheduled chain BayesC/BayesR/BayesRC boundaries and BED cores | marker × trait summaries; full scientific traces; unthinned post-burn convergence | Same first-post-burn summary rule, separate BED seed implementation | Preserve `st_bed_v1` and compare backend reductions |
| ST annotation | `stblr_csr_annot()` plus registered annotation routes | `stblr_cpg_omp_csr_annot()`, group/prior/log-variance/SBayesRC boundaries | common marker shapes plus provider-specific low-dimensional traces | Common convergence contract, heterogeneous provider trace retention | Record each provider trace descriptor explicitly |
| MT CSR | `mtblr_csr()` | `mtblr_csr_chains_raw_internal()`, `run_mt_default_core()` | marker × trait summaries; full `nit + nburn` trait traces; covariance/probability summaries unthinned post-burn | Mixed thinned marker and unthinned covariance aggregation | Preserve `mt_default_v1`; do not normalize silently |
| MT block eigen | `mtblr_block_eigen()` | `mtblr_block_eigen_chains_raw_internal()`, MT default core | MT shapes plus per-owner block metadata | Same MT default retention; operator representation differs | Preserve `mt_default_v1` and provenance |
| MT BED | `mtblr_bed()` | `mtblr_bed_chains_internal()`, `run_mt_bed_bayesc_core()` | marker × trait summaries; full traces; chain-pooled means and primary-chain final state | Matches MT summary aggregation policy but uses distinct chain seed rule | Preserve `mt_bed_v1`; later replace covariance transition |

Current raw marker matrices are $M\times T$. Current full scientific traces are
stored as iteration-by-trait objects, often including burn-in; formal
convergence bundles are unthinned post-burn `iteration × chain × quantity`.
These shapes are evidence for converters, not target raw-v2 axis rules.

Under target retention-contract version 1, post-burn transitions are
$u=1,\ldots,n_{\mathrm{sampling}}$, retained exactly when
$u\bmod n_{\mathrm{thin}}=0$. Thus the indices are
$n_{\mathrm{thin}},2n_{\mathrm{thin}},\ldots$ and the count is the floor in
Section 12. The legacy and target rules coincide when `thin_interval = 1` but
not generally. Phase 1 wrappers preserve legacy indices and label them; the
switch to target retention-contract version 1 requires an explicit migration
checkpoint.

### 23.2 Current RNG semantics

| Route class | Current rule | Classification | Migration |
|---|---|---|---|
| Maintained scalar trait × chain kernels | Without explicit seeds: `seed + 1000003 * (trait_index + 1) + 9176 * (chain_index + 1)` modulo $2^{32}$; the one-chain helper includes the chain-zero `9176` term. With explicit chain bases: `chain_seed[chain] + 1000003 * (trait_index + 1)` | Reproducible but route-specific and trait-order dependent | Preserve as `legacy_st_arithmetic_v1` in Phase 1; later migrate to unified v1 using stable trait IDs |
| Older sequential annotation chain wrappers | Per-chain arithmetic is resolved outside or around single-chain provider kernels; maintained paths remain `std::mt19937` based | Reproducible but route-specific | Record the resolved final seed rather than claim one universal legacy formula |
| Current single-chain MT CSR/default | `std::mt19937(seed)` | Reproducible route-specific | Preserve as `legacy_mt_single_v1` |
| Current MT BED multichain | Chain zero uses `uint32(seed)`; later chains use `uint32(seed + 9176 * chain_index)`; explicit signed values are reinterpreted as final uint32 seeds | Reproducible across the current static chain schedule | Preserve as `legacy_mt_bed_chains_v1` |

Current logical tasks own `std::mt19937` engines, scheduled workers do not call
R RNG, and current static task scheduling supports serial/parallel equality on
qualified routes. Traversal order can still affect draws inside a logical task;
the obsolete MT set loop is additionally posterior-state dependent. Unified
seed version 1 removes worker and trait-position identity from task derivation,
but it does not promise invariance to a scientifically meaningful update-order
change.

## 24. Frozen analysis, execution, and task contract

| Analysis mode | Valid `parallelization` | Logical task ID | Canonical task order |
|---|---|---|---|
| `single_trait` | `none`, `chains` | `(chain_id)` | increasing zero-based chain index |
| `independent_traits` | `none`, `chains`, `traits`, `trait_chains` | `(trait_id, chain_id)` | declared trait-ID order, then increasing chain index |
| `joint_multitrait` | `none`, `chains` | `(chain_id)` | increasing zero-based chain index |

`execution_mode = "serial"` requires `parallelization = "none"`.
`execution_mode = "parallel"` requires a non-`none` value valid for the
analysis mode. No other pair is accepted. Providers, sets, blocks, regions,
and eigenblocks are traversal objects, not execution modes. Traits inside a
joint chain are never separate logical tasks.

Resolved final task seeds have named shape `chain` for `single_trait` and
`joint_multitrait`, and `trait × chain` for `independent_traits`. An explicitly
resolved table bypasses derivation only after exact shape, names, uint32 range,
and uniqueness of task identity have been validated. A public independent-
trait `chain_seeds` convenience vector supplies per-chain base seeds to the
same versioned trait-ID derivation; it is not copied across traits as a final
table.

## 25. Frozen `blr_resolved_spec` version 1

The top-level names, in order, are exactly `schema`, `data`, `model`, `prior`,
`mcmc`, `compute`, and `output`. In the tables, **R** means required, **O** means
optional with a named owner-supplied default, and **N** means required-present
with value `NULL` when inapplicable. `Copy` identifies whether the resolved
value is copied into `blr_raw$input`. All listed fields are copied. `Sci`
identifies fields that can change a posterior, likelihood, transition, or
retained scientific object.

### 25.1 Schema and data fields

| Path | Type/shape | Presence | Valid values/default owner | Validation owner | Sci |
|---|---|---|---|---|---|
| `schema$name` | character scalar | R | exactly `blr_resolved_spec` | schema resolver | no |
| `schema$version` | integer scalar | R | exactly `1` | schema resolver | no |
| `schema$compatibility_id` | character scalar | R | registered Phase 0 compatibility ID | schema resolver | no |
| `schema$seed_contract_version` | integer scalar | R | legacy ID or unified `1` | seed resolver | yes |
| `schema$retention_contract_version` | integer scalar | R | legacy ID or target `1` | MCMC resolver | yes |
| `schema$dimension_contract_version` | integer scalar | R | exactly `1` for raw v2 | schema resolver | no |
| `data$analysis_mode` | character scalar | R | three frozen modes | cross-field resolver | yes |
| `data$trait_ids` | ordered character vector, length $T$ | R | unique, nonempty, non-`NA` | trait resolver | yes |
| `data$global_markers` | ordered character vector, length $M$ | R | unique, nonempty, non-`NA` | marker resolver | yes |
| `data$global_alleles` | table, $M$ rows | R | marker ID plus aligned effect/other allele and coding | alignment resolver | yes |
| `data$operator_resources` | uniquely named list | R | one or more validated immutable resources | resource resolver | yes |
| `data$providers` | uniquely named list | R | one or more providers with nonempty `trait_ids` | provider resolver | yes |
| `data$provider_maps` | provider-named integer vectors | R | resource-local order to global indices; no implicit zeros | marker resolver | yes |
| `data$likelihood_regime` | character scalar | R | `common_sample`, `independent_summary`, `overlap_aware` | likelihood resolver | yes |
| `data$statistical_regions` | named marker-to-region vector or `NULL` | N | `NULL` unless a regional policy is selected | model/data cross-validator | yes |

Operator resources require `resource_id`, `operator_type`, local `marker_ids`,
alleles/coding, centering, standardization, `operator_scale`, immutable storage
or view descriptor, block/eigen metadata or `NULL`, approximation status, and
provenance. Providers require `provider_id`, nonempty ordered `trait_ids`,
`operator_resource_id`, `local_to_global`, phenotype or summary sufficient
statistics, trait-named sample size, likelihood regime, residual/error
contract, population, effect scale, overlap/error group or `NULL`, and
provenance. Resource resolution validates storage and operator algebra;
provider resolution validates statistics, maps, trait ownership, scales, and
likelihood regime. C++ receives immutable views. Several providers may
reference one resource without duplicating BED or LD storage.

A common-sample full-$V_e$ MT likelihood uses one provider owning all coupled
traits. It must not be factorized into singleton providers. Independent
summary providers normally own one trait and may have different
$X_t^\top X_t$, sample sizes, marker sets, populations, blocks, eigenvectors,
eigenvalues, and ranks. Missing markers contribute no likelihood term.

### 25.2 Model and prior fields

| Path | Type/shape | Presence | Contract/default owner | Validation owner | Sci |
|---|---|---|---|---|---|
| `model$family` | character scalar | R | registered family | model registry | yes |
| `model$state_space` | ordered named table | R | exact state/pattern/component definition | model registry | yes |
| `model$null_state_index` | integer scalar | R | one-based index into `state_space` | model registry | yes |
| `model$effect_storage_convention` | character scalar | R | `realised`, `base_latent`, or `scaled_latent`, with realised output always defined | model registry | yes |
| `model$probability_policy` | named descriptor | R | explicit inclusion/component/activity policy | probability resolver | yes |
| `model$marker_scale_policy` | named descriptor | R | unit, component, MAF-S, annotation, or external $q_j$ | scale resolver | yes |
| `model$marker_covariance_policy` | named descriptor | R | scalar, global matrix, regional, or template | covariance resolver | yes |
| `model$residual_policy` | named descriptor | R | scalar, diagonal, fixed full, sampled full, or overlap-aware | residual resolver | yes |
| `model$update_order_version` | integer scalar | R | registered transition version | model registry | yes |
| `prior$probability` | named list or `NULL` | N | policy-specific prior | probability resolver | yes |
| `prior$component_multipliers` | named numeric vector or `NULL` | N | finite nonnegative $\gamma_k$ with declared null | scale resolver | yes |
| `prior$marker_multipliers` | marker-named numeric vector or `NULL` | N | finite positive $q_j$ when applicable | scale resolver | yes |
| `prior$scalar_variance` | named list or `NULL` | N | explicit shape/scale convention | variance resolver | yes |
| `prior$marker_covariance` | named list or `NULL` | N | `degrees_of_freedom`, `scale`, and sampled/fixed flag | covariance resolver | yes |
| `prior$residual_covariance` | named list or `NULL` | N | explicit fixed value or inverse-Wishart parameters | residual resolver | yes |
| `prior$annotation` | named list or `NULL` | N | registered annotation-policy prior | annotation resolver | yes |

### 25.3 MCMC, compute, and output fields

| Path | Type/shape | Presence | Contract/default owner | Validation owner | Sci |
|---|---|---|---|---|---|
| `mcmc$burn_in_iterations` | integer scalar | R | $\geq0$; wrapper adapter supplies | shared MCMC resolver | yes |
| `mcmc$sampling_iterations` | integer scalar | R | $\geq1$ | shared MCMC resolver | yes |
| `mcmc$thin_interval` | integer scalar | R | $\geq1$ | shared MCMC resolver | yes |
| `mcmc$retained_draws` | integer scalar | R | derived floor under target v1 | shared MCMC resolver | no |
| `mcmc$retained_transition_indices` | integer vector | R | exact indices under declared retention version | shared MCMC resolver | no |
| `mcmc$chains` | integer scalar | R | $\geq1$ | shared MCMC resolver | yes |
| `mcmc$seed` | numeric integer scalar | R | uint32 user seed for unified v1 | seed resolver | yes |
| `mcmc$task_seeds` | named vector or matrix | R | final native uint32 seeds with contracted axes | seed resolver | yes |
| `mcmc$update_flags` | uniquely named logical list | R | model-registry defaults, all recorded | model/MCMC cross-validator | yes |
| `compute$execution_mode` | character scalar | R | `serial` or `parallel` | compute resolver | no |
| `compute$parallelization` | character scalar | R | closed combination table in Section 24 | compute resolver | no |
| `compute$cores` | integer scalar | R | $\geq1$ | compute resolver | no |
| `compute$scheduler_version` | integer scalar | R | registered deterministic scheduler | compute resolver | no |
| `compute$memory_limit_bytes` | finite nonnegative scalar, `Inf`, or `NULL` | N | output resolver default | memory resolver | no |
| `compute$operator_numerical_controls` | uniquely named list | R | resource-specific registered defaults | resource resolver | potentially |
| `output$posterior_summaries` | logical scalar | R | `TRUE` | output resolver | no |
| `output$retained_parameters` | unique character vector | R | model-applicable names | output resolver | no |
| `output$effect_draw_policy` | character scalar | R | `none`, `selected`, `compact`, or `full` | output resolver | no |
| `output$state_draw_policy` | character scalar | R | `none`, `compact`, or `full` | output resolver | no |
| `output$convergence_policy` | named descriptor | R | current convergence resolver | output resolver | no |
| `output$derived_quantities` | unique character vector | R | empty unless requested | output resolver | no |
| `output$preserve_chains` | logical scalar | R | `TRUE` for raw v2 | output resolver | no |
| `output$memory_estimate_bytes` | finite nonnegative scalar | R | derived | memory resolver | no |

Every list at every scientific forwarding boundary has exact unique nonempty,
non-`NA` names. Unnamed, empty, duplicate, partial, positional, and unknown
scientific arguments fail before native dispatch. Public wrappers may keep
their flat formals initially, but they must resolve into this exact internal
form rather than forward `...` positionally.

## 26. Frozen `blr_raw` schema version 2

The top-level order is exactly `schema`, `model`, `input`, `posterior`,
`draws`, `final`, `derived`, `diagnostics`, and `provenance`. `schema`, `model`,
`input`, `posterior`, `draws`, `final`, `derived`, `diagnostics`, and
`provenance` are required. Model-inapplicable scientific fields listed below
are required-present-`NULL`; assignment must preserve their names.

| Namespace | Required fields | Required-present-`NULL` when inapplicable |
|---|---|---|
| `schema` | `name = "blr_raw"`, `version = 2`, `compatibility_id`, `dimension_contract_version` | none |
| `model` | copied resolved family, modes, state space, null index, covariance/residual/probability policies, effect storage, update-order version | none |
| `input` | complete validated `blr_resolved_spec` v1 | none |
| `posterior` | `realised_effect_mean`, `pips` | `latent_effect_mean`, `scaled_effect_mean`, `traitwise_state_probabilities`, `joint_state_probabilities`, `activity_pattern_probabilities`, `traitwise_component_assignment_probabilities`, `joint_component_assignment_probabilities`, `marker_variance_mean`, `marker_covariance_mean`, `residual_variance_mean`, `residual_covariance_mean`, `uncertainty` |
| `draws` | `realised_effects`, `traitwise_activity` when retained by output policy | `latent_effects`, `scaled_effects`, `independent_trait_states`, `joint_states`, `traitwise_probability_parameters`, `joint_probability_parameters`, `activity_pattern_parameters`, `traitwise_component_probability_parameters`, `joint_component_probability_parameters`, `marker_variance`, `marker_covariance`, `residual_variance`, `residual_covariance`, `regional_marker_covariance`, `convergence` |
| `final` | `realised_effects` plus the applicable state and variance/covariance state | `latent_effects`, `scaled_effects`, `independent_trait_states`, `joint_states`, `traitwise_probability_parameters`, `joint_probability_parameters`, `activity_pattern_parameters`, `traitwise_component_probability_parameters`, `joint_component_probability_parameters`, `marker_variance`, `marker_covariance`, `residual_variance`, `residual_covariance`, `rng_continuation` |
| `derived` | namespace present | `predictions`, `genetic_variance`, `genomic_covariance`, `operator_relative_quadratics`, `descriptive_bilinear_forms` |
| `diagnostics` | namespace present | `convergence`, `acceptance`, `runtime`, `memory`, `workers`, `numerical_safeguards`, `approximation_warnings` |
| `provenance` | `package_version`, `operator_resources`, `marker_alignment`, `seed_contract_version`, `task_seeds` | `git_sha`, `dirty_build`, `compiler`, `timestamp` when unavailable |

Required array axes are fixed:

| Field class | Axis order |
|---|---|
| retained realised/latent/scaled effects | `draw × chain × marker × trait` |
| independent-trait states | `draw × chain × marker × trait` |
| joint states | `draw × chain × marker` |
| traitwise activity | `draw × chain × marker × trait` |
| PIPs | `marker × trait` |
| traitwise state probabilities | `marker × trait × state` |
| joint-state probabilities | `marker × joint_state` |
| activity-pattern probabilities | `marker × activity_pattern` |
| traitwise component-assignment probabilities | `marker × trait × component` |
| joint component-assignment probabilities | `marker × component` |
| traitwise probability-parameter draws | `draw × chain × trait × state` |
| joint probability-parameter draws | `draw × chain × joint_state` |
| activity-pattern parameter draws | `draw × chain × activity_pattern` |
| traitwise component-probability parameter draws | `draw × chain × trait × component` |
| joint component-probability parameter draws | `draw × chain × component` |
| marker/residual covariance draws | `draw × chain × trait_row × trait_col` |
| scalar trait variance draws | `draw × chain × trait` |
| final effects and independent states | `chain × marker × trait` |
| final joint states | `chain × marker` |
| final traitwise probability parameters | `chain × trait × state` |
| final joint probability parameters | `chain × joint_state` |
| final activity-pattern parameters | `chain × activity_pattern` |
| final traitwise component-probability parameters | `chain × trait × component` |
| final joint component-probability parameters | `chain × component` |
| final covariance states | `chain × trait_row × trait_col` |
| regional marker-covariance draws | `draw × chain × region × trait_row × trait_col` |
| predictions | `draw × chain × observation × trait` when retained; summaries use `observation × trait` |
| provider diagnostics | `provider` first, followed by a diagnostic-specific named axis; no anonymous list order |

Every axis has its canonical axis name plus stable dimnames. Size-one axes are
never dropped in raw v2. Native memory may remain task-major; `RawResultBuilder`
owns the checked conversion. The formatted layer may explicitly simplify
dimensions while recording that transformation.

Raw v2 retires unqualified `pi`, `pis`, `pim`, `state_probabilities`, and
`pattern_probabilities`. The explicit concepts are:

| Name | Fixed meaning |
|---|---|
| `prior_state_mass` | resolved prior mass before marker data |
| `traitwise_probability_parameters` | sampled trait-specific state-mass parameters |
| `joint_probability_parameters` | sampled joint-state simplex parameters |
| `activity_pattern_parameters` | sampled Cheng activity-pattern simplex parameters |
| `traitwise_state_probabilities` | markerwise posterior probabilities over traitwise states |
| `joint_state_probabilities` | markerwise posterior probabilities over declared joint states |
| `activity_pattern_probabilities` | markerwise posterior probabilities over activity patterns |
| `traitwise_component_assignment_probabilities` | marker-by-trait posterior component allocation probabilities |
| `joint_component_assignment_probabilities` | markerwise posterior shared/joint component allocation probabilities |
| `pips` | marker-by-trait posterior non-null probabilities |
| `joint_non_null_probabilities` | markerwise probability that at least one trait is active |
| `pleiotropic_probabilities` | markerwise probability of the explicitly declared pleiotropic pattern set |

A temporary formatted alias is a read-only view of exactly one documented v2
field. There is no universal `fit$pis`.

Git provenance is resolved at build/load time, not by invoking Git for every
fit. Development builds record the SHA and dirty indicator when the build
system supplies them. Installed release packages use package version and build
ID and store `git_sha = NULL` and `dirty_build = NULL` if unavailable.

## 27. Covariance and initial residual contracts

For $V_b\sim\operatorname{IW}_T(\nu_b,\Psi_b)$, the density is proportional to

$$
|V_b|^{-(\nu_b+T+1)/2}
\exp\left\{-\frac12\operatorname{tr}(\Psi_bV_b^{-1})\right\}.
$$

The only canonical internal/public parameters are
`degrees_of_freedom = nu_b` and `scale = Psi_b`. A distribution is proper for
$\nu_b>T-1$. Its mean is finite only for $\nu_b>T+1$, in which case a helper
may accept a desired mean $M_b$ and record the conversion
$\Psi_b=(\nu_b-T-1)M_b$ in the resolved spec and provenance. No backend may
reinterpret the scale as a mean.

Phase 4a uses a supplied fixed symmetric positive-definite full $V_e$ with a
common-sample multi-trait BED provider. Raw v2 records the fixed policy and
final fixed state but does not manufacture posterior $V_e$ draws. Phase 4b
adds sampled full inverse-Wishart $V_e$ and is required before the new MT route
replaces maintained production MT. A diagonal $V_e$ is an explicit common-
sample reduction. Independent no-overlap summaries use diagonal/fixed marginal
residual policies because they do not identify off-diagonal residual
covariance. Overlap-aware error covariance remains separately gated.

The current MT covariance hybrid is not a sampled schema-v2 $V_b$. Its
full-latent/set updates and later heuristic replacement cannot be converted
into an authoritative draw; affected analyses require rerunning under the
corrected model.

## 28. Wrapper and raw-v1 migration

### 28.1 Current argument mapping

| Current argument | Target path | Phase 1 action |
|---|---|---|
| `nit` | `mcmc$sampling_iterations` | Preserve current transition count and legacy retention version |
| `nburn` | `mcmc$burn_in_iterations` | Exact alias |
| `nthin` | `mcmc$thin_interval` | Preserve current first-post-burn rule until explicit retention migration |
| `nchains` | `mcmc$chains` | Exact validated alias |
| `seed` | `mcmc$seed` | Preserve route-specific legacy seed version initially |
| `chain_seeds` | wrapper convenience, then final `mcmc$task_seeds` | Preserve current route behavior initially; never treat an independent-trait chain vector as a final trait × chain table under unified v1 |
| `ncores` | `compute$cores` plus execution-policy resolution | Derive valid mode/policy explicitly |
| `keep_chains` | `output$preserve_chains` and compact draw policy | Preserve current request, but raw v2 always preserves chain axes |
| `updateB`, `updateE`, `updatePi`, provider flags | `mcmc$update_flags` with policy-specific names | Exact-name model resolver; no generic partial matching |
| `method` | `model$family` plus data-level/provider semantics | Resolve aliases explicitly; reject ambiguous values |
| `stats`, `Glist`, BED inputs, LD/eigen inputs | `data$operator_resources`, `data$providers`, and maps | Separate reusable representation from likelihood statistics |
| `pi`, `joint_pi`, component and annotation controls | `model$probability_policy` and `prior$probability` | Split by scientific meaning; no generic internal `pi` |
| `vb`, `ve`, `ssb_prior`, `sse_prior`, `nub`, `nue` | explicit initial states and covariance-prior objects | Convert only where parameterization is source-known; otherwise fail |

`niters` and `nsamples` are not accepted target aliases. Structured public
arguments are deferred; current wrappers are compatibility adapters into the
structured internal spec during Phase 1.

### 28.2 Raw and formatted migration

| Current v1 field | V2 treatment |
|---|---|
| `marker$bm`, formatted `bm` | rename to `posterior$realised_effect_mean`; expand trait axis if needed |
| `marker$dm`, formatted `dm` | rename to `posterior$pips` only where it is verified binary non-null probability |
| `marker$b`, `marker$beta`, final fields | split into realised, latent, or scaled `final` fields under declared storage convention |
| `marker$state`, `component_final` | split into independent-trait or joint-state fields; never infer from rank alone |
| `component_probabilities` | map only when verified as markerwise posterior allocation probabilities |
| `pi_*`, `pis`, `pim`, pattern fields | split into explicit sampled-parameter and markerwise posterior fields; ambiguous objects are unsupported |
| `vbs`, `ves` | map to scalar variance draws only with known iteration/chain ownership; otherwise require rerun |
| `cov_b_*`, `cov_e_*` | expand to explicit draw/final axes only when they are actual sampled states with known semantics |
| `cov_g_*`, `vgs`, `vle`, `vld` | place under derived quantities using identified or operator-relative names; do not relabel descriptive bilinear forms as covariance |
| convergence bundle v1 | retain through a versioned converter when iteration and chain identities are explicit |
| diagnostic/provider metadata | normalize into diagnostics/provenance without changing scientific arrays |

A compatibility wrapper resolves new calls. A schema converter accepts a
serialized v1 object only when field meaning, axes, retention, and chain
ownership are known from its schema and compatibility ID. A formatted alias is
a read-only view. A scientific migration requires rerunning. There is no
positional fallback, no dimension-only guessing, and unsupported versions fail
with the exact missing semantic requirement.

## 29. Phase 0 fixture evidence

The base-R files under `tests/research/blr_framework_contract/` are independent
contract fixtures, not production validators. They cover one-chain/one-axis
single-trait raw arrays, two-trait/two-chain independent and joint modes, one
BED resource shared by several providers, one joint multi-trait provider,
heterogeneous singleton summary providers, required-present-`NULL` fields,
fixed and sampled full $V_e$, exact retention indices, execution-policy
failures, and fixed task-seed vectors.

Seed-contract version-1 reference vectors are:

| User seed | Analysis mode | Trait ID | Zero-based chain | Native uint32 seed |
|---:|---|---|---:|---:|
| 0 | `single_trait` | sentinel | 0 | 830191578 |
| 0 | `single_trait` | sentinel | 1 | 160141543 |
| 0 | `independent_traits` | `traitA` | 0 | 226943096 |
| 0 | `independent_traits` | `traitB` | 1 | 286956759 |
| 0 | `joint_multitrait` | sentinel | 0 | 3100589946 |
| 17 | `single_trait` | sentinel | 0 | 3397578794 |
| 17 | `independent_traits` | `traitB` | 1 | 1132619387 |
| 17 | `joint_multitrait` | sentinel | 1 | 3700933392 |

Stable trait-ID reordering leaves each independent trait's row unchanged.
Serial versus parallel execution consumes the same resolved table.

## 30. Phase 0 completion table

| Contract | Decision | Evidence | Remaining implementation phase |
|---|---|---|---|
| Analysis/execution | Separate closed enumerations and task identities | Source task topology plus executable fixtures | Phase 1 resolver; Phase 3 scheduler |
| Retention | Target version 1 uses divisible post-burn indices; legacy rules remain versioned | Native source trace and exact fixture calculations | Phase 1 adapters; deliberate migration checkpoint |
| RNG | Exact FNV-1a/SplitMix64 seed-contract version 1; zero valid | Constants, pseudocode, and eight fixed vectors | Phase 3 implementation; wrapper trajectory migration later |
| Resolved input | `blr_resolved_spec` version 1 with seven exact namespaces | Field tables and fixtures | Phase 1 |
| Resource/provider | Immutable reusable resources; providers own nonempty trait sets and statistics | Shared BED and joint-provider fixtures | Phase 2 |
| Raw result | `blr_raw` version 2 with fixed names and axes | Single/independent/joint raw fixtures | Phase 1 builder and converters |
| Dimensions | Draw × chain first; no size-one dropping | Exact dim/dimname assertions | Phase 1 builder/validator |
| Probability names | Explicit parameter versus posterior concepts; no universal `pis` | Schema fixture rejects ambiguous names | Phase 1 schema/formatter migration |
| Inverse-Wishart | Explicit degrees of freedom and scale; mean conversion recorded | Density and propriety/mean conditions | Phase 4 covariance policies |
| Initial MT residual | Fixed full SPD in Phase 4a, sampled full IW in Phase 4b before promotion | Research oracle design and identifiability boundary | Phase 4a/4b |
| V1 migration | Convert only source-known semantics; current hybrid requires rerun | Current raw/source audit | Phase 1 converters and release notes |

All Phase 0 decisions have recommendations. This Phase 0 table does not by
itself mark a production phase implemented.

## 31. Phase 1 R-infrastructure checkpoint

Phase 1 implements the R-side contracts without changing native interfaces or
scientific kernels:

| Contract | Implemented evidence | Boundary retained |
|---|---|---|
| Resolved specification | `R/blr-resolved-spec.R`; exact namespaces, names, resources/providers, execution policy, task-seed shapes, and legacy retention indices | Current wrappers and native arguments remain unchanged |
| Raw schema v2 | `R/blr-raw-v2.R`; fixed axes, required-present-`NULL`, probability semantics, covariance checks, and cached provenance | Native backends still return schema v1 |
| ST conversion | Explicit one-chain BayesC/BayesR/BayesRC-family converter used by maintained CSR, BED, block-eigen, and supported annotation wrappers | Multi-chain final effects unavailable in v1; those fits remain legacy |
| Formatting | Validated v2 precedes Phase 1 alias construction; each alias records one schema-v2 source | No universal ambiguous `pis` field |
| Legacy trajectories | Seed and retention contract version 0 plus exact final task seeds and retained indices | Unified seed and retention version 1 remain Phase 3 |
| Provenance | Once-per-session cached package/load metadata; absent Git metadata is present with `NULL` | Build-time Git injection remains a promotion task |
| MT restriction | Current MT raw is rejected by the v2 converter | Corrected Cheng MT implementation remains Phase 4 |

Phase 1 does not implement shared native operator adapters, unified native
scheduling, the unified seed derivation, or a corrected joint-MT covariance
transition. Those roadmap boundaries remain unchanged.
