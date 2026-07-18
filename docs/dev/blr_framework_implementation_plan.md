# Unified BLR Framework: Architecture and Incremental Refactoring Plan

**Status:** Phase 13B packed-BED BayesR per-chain block extracted; migration in progress
**Date:** 2026-07-18
**Target location:** `docs/dev/blr_framework_implementation_plan.md`

## Current framework status

Canonical ordinary CSR routes are BayesC, BayesR, SBayesRC, fixed-prior
BayesC, group BayesC, and learned-annotation BayesC. Canonical scheduled routes
are ordinary-CSR BayesC and the public scheduled packed-BED BayesC route.
Experimental/sparse BED BayesC, packed-BED BayesR, packed-BED BayesRC,
block-eigen routes where applicable, and legacy multivariate implementations
remain audited or protected but noncanonical. Later phase sections are the
current record; earlier status statements are historical snapshots.

Phase 13A leaves packed-BED BayesR production unchanged while establishing its
route/statistical contract, component ordering, scheduler and RNG ownership,
aggregation/schema inventory, deterministic references, and future per-chain
extraction seam.

Phase 13B mechanically relocates the deterministic per-chain helper and its
binding-neutral result/helper vocabulary into one guarded implementation header.
Logical-chain RNG, component and scheduler semantics, genotype ownership, inline
task aggregation and inline conversion are unchanged. A typed callable boundary
is not yet active.

## 1. Purpose

This document defines the implementation strategy for a maintainable, extensible, and language-neutral Bayesian linear regression framework in `sblr`.

The current package already contains fast and memory-efficient native implementations for several important single-trait models and data representations. Those implementations represent substantial validated work and should not be rewritten merely to obtain a cleaner architecture.

The selected strategy is therefore:

> Preserve efficient and statistically validated computational kernels, refactor their boundaries and shared infrastructure incrementally, and rewrite only components that are statistically incoherent, structurally unsuitable, or impossible to extend safely.

The current implementations serve as:

- production methods during migration;
- statistical reference implementations;
- performance and memory baselines;
- sources of validated algorithms and fixtures;
- migration targets for shared infrastructure.

The long-term result should contain one canonical implementation for each supported model/operator combination, but the path to that result is an incremental architectural refactor rather than a package-wide greenfield rewrite.

---

## 2. Priority order

The framework should optimize for the following priorities, in order:

1. statistical correctness;
2. ease of maintenance;
3. ease of extension;
4. preservation of current speed;
5. preservation of current memory efficiency;
6. reproducibility;
7. language-neutrality for future R and Python interfaces;
8. performance improvements where profiling identifies a real bottleneck.

Speed and memory are regression constraints during migration. A refactor is not required to be faster than the current native code, but it must not introduce an unexplained material regression.

---

## 3. Core principles

### 3.1 Preserve validated hot loops

Working marker-update kernels should remain unchanged initially when they are:

- statistically coherent;
- reproducible;
- well tested;
- already fast and memory efficient.

Refactoring should first occur around them:

- typed model specifications;
- data/operator ownership;
- chain and RNG infrastructure;
- result accumulation;
- result conversion;
- probability and scale utilities;
- common validation.

Only after exact or statistical equivalence is established should a hot loop be changed.

### 3.2 Rewrite only where necessary

A clean replacement is appropriate when the existing implementation is:

- statistically inconsistent;
- based on conflicting updates;
- impossible to validate through reductions;
- tightly coupled to obsolete assumptions;
- duplicated in ways that prevent reliable extension;
- non-reproducible in a way that cannot be fixed locally.

The legacy multivariate implementation falls into this category and should not become the canonical new multivariate engine.

### 3.3 Share infrastructure, not every instruction

The framework should share:

- specifications;
- operators;
- RNG;
- chain scheduling;
- state-probability utilities;
- scale policies;
- covariance utilities;
- posterior accumulators;
- diagnostics;
- typed results;
- binding conversion.

It should not force statistically different marker kernels into one universal runtime-polymorphic implementation.

### 3.4 One canonical implementation after migration

Temporary old/new duplication is acceptable during validation.

After a model has:

- migrated to the shared architecture;
- passed statistical and reproducibility tests;
- passed performance and memory regression gates;
- stabilized in production use;

the superseded legacy implementation should be removed.

### 3.5 Thin language bindings

The C++ statistical core should be independent of R and Python.

Rcpp should remain a thin current binding. A future Python interface should use the same core through a thin binding such as `pybind11`.

No Python package or dependency is required during the current refactor.

---

## 4. Phase 0 findings

The repository audit in `docs/dev/blr_framework_phase0_audit.md` established the following architectural constraints.

### 4.1 Strong reusable components

The strongest reusable current components include:

- unscheduled single-trait CSR BayesC;
- CSR BayesR and SBayesRC chain infrastructure;
- the CSR/block-eigen operator abstraction;
- packed-BED decoding and score/update helpers;
- current efficient residual updates;
- chain-task seed helpers;
- recent chain-local RNG implementations;
- online posterior accumulation;
- named raw output and common formatting;
- BayesRC stick-breaking utilities;
- typed sufficient-statistic preparation examples.

These components should be migrated or extracted, not automatically replaced.

### 4.2 Legacy multivariate code is not canonical

The current multivariate routes contain unresolved statistical conflicts, including:

- multiple incompatible covariance updates in one iteration;
- inconsistent use of update flags;
- incomplete use of full residual precision;
- unclear active-subvector versus masked-latent-vector priors;
- positional output semantics;
- no automated multivariate reduction suite.

The canonical multivariate model must therefore be newly specified and implemented.

### 4.3 Scheduled BayesC RNG requires a scoped correction

Some scheduled BayesC paths contain persistent distribution state.

This should be fixed as a separate reproducibility task. It should not be mixed into the first architecture migration.

### 4.4 Full multivariate residual covariance needs an explicit data contract

The framework must distinguish:

- same individuals and same design;
- trait-specific designs;
- partially overlapping studies;
- independent studies;
- summary-statistics inputs;
- evidence matrices.

Full residual covariance may be enabled only when the supplied data identify cross-trait residual products.

---

## 5. Target architecture

```text
R interface now / Python interface later
                 |
                 v
       user-facing validation and alignment
                 |
                 v
       typed language-neutral specification
                 |
                 v
          static model/operator dispatch
                 |
                 v
   shared native infrastructure and operators
                 |
       +---------+----------+
       |                    |
validated migrated      new kernels where
existing kernels        redesign is required
       |                    |
       +---------+----------+
                 |
                 v
         typed native result
                 |
                 v
        thin language conversion
```

The architecture has five major layers:

1. binding and public interface;
2. resolved specification;
3. shared execution infrastructure;
4. optimized model kernels;
5. typed result and formatting.

---

## 6. Kernel families

The framework should use a small number of optimized kernel families.

### 6.1 Scalar regression kernels

Used for:

- single-trait BayesC;
- BayesR;
- BayesRC;
- fixed marker priors;
- `selection_s`;
- group and hierarchy models;
- learned annotation models.

Existing efficient scalar kernels should be migrated incrementally.

They may remain specialized if merging them would reduce clarity or performance.

### 6.2 Small-\(T\) multivariate kernel

Used for:

- explicit multivariate BayesC patterns;
- multivariate mixtures;
- hierarchical multivariate models;
- full or diagonal covariance policies.

This kernel should be newly implemented from a coherent mathematical model.

### 6.3 Factor kernels

Used for:

- factor-analytic multivariate BLR;
- sparse marker-factor effects;
- annotation-informed gene/SNP evidence factors;
- multi-view factor models.

Regression and evidence-matrix likelihoods may share factor utilities but remain distinct frontends.

---

## 7. Resolved model specification

A fit should be represented by typed specifications before entering the sampler.

Core specification concepts:

```text
DataSpec
ModelSpec
McmcControl
OutputSpec
ExecutionSpec
```

### 7.1 `DataSpec`

Defines:

- data representation;
- marker and trait counts;
- marker and trait order;
- sample/design contract;
- scaling convention;
- sample-size information;
- operator ownership;
- optional cross-trait sufficient statistics.

### 7.2 `ModelSpec`

Defines:

- kernel family;
- state family;
- state-probability policy;
- marker-scale policy;
- trait covariance policy;
- residual covariance policy;
- hyperparameters;
- fixed/update flags.

### 7.3 `McmcControl`

Defines:

- iterations;
- burn-in;
- thinning;
- number of chains;
- seeds and optional explicit chain seeds;
- number of cores;
- residual rebuild schedule;
- validation controls.

### 7.4 `OutputSpec`

Defines requested output.

Illustrative structure:

```cpp
struct OutputSpec {
  bool marker_mean = true;
  bool marker_pip = true;
  bool component_probabilities = false;
  bool parameter_traces = true;
  bool final_state = true;
  bool final_residual = false;
  bool keep_chain_summaries = false;
  bool keep_marker_chain_summaries = false;
  bool cpo = false;
};
```

Existing output behavior should be preserved during migration. More selective output can be introduced only with explicit tests and documentation.

---

## 8. Language-neutral C++ boundary

The native statistical core should use:

- standard C++ types;
- Armadillo value types where useful;
- standard exceptions;
- chain-owned RNG.

It should not use:

- `Rcpp::List`;
- `Rcpp::NumericMatrix`;
- `Rcpp::stop`;
- `Rcpp::Rcout`;
- R RNG;
- R normal CDF/quantile calls;
- Python or NumPy types.

The binding flow should be:

```text
R arguments
  -> R validation and alignment
  -> Rcpp adapter
  -> typed C++ specification
  -> native execution
  -> typed C++ result
  -> R result conversion
```

A future Python interface should use the same middle layers.

A stable public C++ ABI is not required during refactoring.

---

## 9. Static dispatch and performance preservation

Model and operator choices should be resolved before the MCMC loop.

A small outer dispatch may select a concrete kernel/policy combination.

Avoid within marker loops:

- repeated enum switches;
- virtual calls;
- `std::function`;
- binding callbacks;
- input validation;
- avoidable allocation.

However, existing fast code should not be templated or reorganized solely for stylistic uniformity. The benefit of extraction must be demonstrated through maintenance or extension value.

---

## 10. Data and operator layer

The operator represents immutable likelihood data and optimized residual operations.

A conceptual scalar operator supports:

```text
diag(marker)
score(marker, chain_state)
apply_effect_change(marker, delta, chain_state)
rebuild(chain_state)
```

Optional capabilities may include:

- LD-neighbor enumeration;
- direct individual residual access;
- cross-trait residual products;
- block diagnostics;
- memory-mapped access.

Planned operators:

- CSR;
- packed BED;
- block eigen;
- dense individual;
- dense summary;
- factor evidence frontend.

The existing CSR/block-eigen abstraction should be retained and improved incrementally.

The existing BED decoder and score/update helpers should be separated from Rcpp conversion without rewriting the underlying optimized decoding logic.

---

## 11. Chain state, workspace, and ownership

Each chain owns mutable state:

- effects;
- residual or score residual;
- states/components;
- sampled parameters;
- posterior accumulators;
- scratch workspace;
- RNG and distributions.

Immutable data should be shared safely across chains:

- CSR arrays;
- BED bytes;
- marker diagonals;
- allele frequencies;
- annotations;
- hierarchy indices;
- component scales;
- trait-pattern definitions.

No architecture refactor should introduce copying of full operator data per chain.

### 11.1 Workspace rule

New or substantially rewritten kernels should preallocate workspaces.

Existing kernels need not be rewritten merely to satisfy a stylistic no-allocation rule, but any material allocation inside their hot loop should be measured before extension.

### 11.2 Posterior accumulation

Default marker summaries should be accumulated online.

Full marker-by-iteration storage should not be introduced.

---

## 12. RNG and parallelism

Each chain must own:

- its engine;
- all stateful distributions;
- proposal distributions;
- marker-order state.

No stochastic state may be:

- global;
- static;
- persistent thread-local;
- shared between chains;
- supplied by R or Python;
- silently supplied by Armadillo global RNG.

Parallelism should initially use independent trait-chain or chain tasks.

A chain's trajectory must not depend on the worker thread executing it.

Existing kernels with validated deterministic behavior should preserve their RNG ordering during migration.

If a kernel's RNG is corrected, that correction must be isolated, characterized, and tested separately from architectural refactoring.

---

## 13. Typed result model

The long-term core result should be typed.

Illustrative structure:

```cpp
struct BlrResult {
  MarkerSummary marker;
  StateSummary state;
  TraceSummary trace;
  VarianceSummary variance;
  Diagnostics diagnostics;

  OptionalComponentSummary component;
  OptionalPatternSummary pattern;
  OptionalAnnotationSummary annotation;
  OptionalHierarchySummary hierarchy;
  OptionalFactorSummary factor;

  std::vector<ChainSummary> chains;
};
```

Canonical dimensions:

```text
marker effects:            markers × traits
marker PIP:                markers × traits
traces:                    retained samples × parameter dimension
trait covariance:          traits × traits
component probability:     optional markers × components × traits
pattern probability:       optional markers × patterns
factor row loadings:       rows × factors
factor trait loadings:     traits × factors
```

During migration, the typed result may be converted into the existing `stblr_raw_v1` representation.

A new public result schema should be introduced only after enough model families have migrated to determine stable common semantics.

---

## 14. Model composition

A general marker prior can be represented as:

\[
P(z_j\mid A_j,\Theta),
\]

\[
\mathbf b_j\mid z_j=(r,k)
\sim
N_T\left(0,\;q_jD_r\Sigma_kD_right).
\]

This separates:

- discrete state;
- state probabilities;
- marker scale;
- covariance;
- likelihood.

This is a conceptual decomposition. It does not require one universal marker-update function.

---

## 15. Marker-scale policies

Initial scale policies:

- unit;
- fixed marker scale;
- MAF-dependent `selection_s`;
- group scale;
- categorical hierarchy;
- component scale;
- compatible composite scale.

For a categorical hierarchy:

\[
q_j=\prod_{h=1}^{H}\theta_{h,g_h(j)}.
\]

A scale-preserving weighted geometric normalization is:

\[
c_h =
\exp\left[
\frac{\sum_g n_{hg}\log\theta_{hg}}
{\sum_g n_{hg}}
ight],
\]

\[
\theta_{hg}\leftarrow\theta_{hg}/c_h,
\qquad
v_b\leftarrow v_bc_h.
\]

This preserves \(v_bq_j\).

The current group model's normalization should remain protected as legacy behavior until the hierarchy model is implemented and its intended reduction is explicitly decided.

---

## 16. Multivariate strategy

The canonical multivariate implementation should be newly designed.

Initial model requirements:

- small \(T\);
- exact shared-design likelihood;
- explicit restricted trait patterns;
- fixed, diagonal, and one approved full covariance policy;
- complete use of the selected residual precision;
- explicit active-subvector or masked-latent prior;
- no automatic large \(2^T\) state enumeration;
- exact \(T=1\) reduction to scalar BLR;
- diagonal reduction to independent scalar fits.

The initial data operator should use exact dense shared-design sufficient statistics or exact individual-level data.

CSR and other approximations should follow only after the reference model is validated.

---

## 17. Factor strategy

The factor framework should share:

- membership policies;
- annotation effects;
- factor scales;
- normalization;
- sign conventions;
- ordering;
- cross-chain alignment;
- posterior accumulation.

Separate frontends:

1. factor-analytic regression;
2. evidence-matrix factor model.

The evidence model should begin with gene-by-trait input. SNP-to-gene and gene-to-pathway propagation should initially be posterior post-processing.

---

## 18. Migration classes

Every existing component should be classified as one of:

### Preserve and wrap

Use when code is efficient, coherent, and well tested.

Examples:

- unscheduled CSR BayesC hot loop;
- current efficient residual updates;
- selected BayesR/BayesRC kernels;
- BED decoding;
- operator storage.

### Extract shared infrastructure

Use when the same concept is implemented repeatedly.

Examples:

- chain/RNG ownership;
- stable log-weight normalization;
- categorical sampling;
- variance utilities;
- output accumulation;
- state counts;
- typed result assembly;
- scale log-prior calculations.

### Correct in place

Use for localized defects.

Examples:

- persistent stateful distributions;
- inconsistent `NULL` conversion;
- binding-specific math in a reusable utility.

### Rewrite

Use when existing code cannot be made canonical safely.

Examples:

- legacy multivariate BLR;
- conflicting covariance machinery;
- obsolete positional MT output;
- experimental implementations replaced by validated shared infrastructure.

### Remove after stabilization

Use when the migrated implementation has passed all gates.

---

## 19. Migration acceptance gates

### 19.1 Statistical gate

- unchanged model preserves its posterior target;
- reduction tests pass;
- output identities pass;
- no unexplained discrepancy remains.

### 19.2 Reproducibility gate

- repeated fixed-seed equality;
- one-core/multi-core equality;
- reversed call order;
- unrelated intervening fit;
- explicit chain seed behavior.

### 19.3 Performance gate

For a migration of an already efficient backend:

- no unexplained material runtime regression;
- representative workloads benchmarked;
- initialization and binding overhead measured separately where relevant.

A default practical threshold is approximately 5–10% runtime variation, interpreted with benchmark uncertainty. A larger regression requires explicit approval and justification.

### 19.4 Memory gate

- no meaningful peak-memory regression;
- no operator data copied per chain;
- no new full marker-sample storage;
- optional outputs measured separately.

### 19.5 Maintenance gate

- responsibilities are clearer than before;
- extension points are documented;
- duplicate infrastructure is reduced;
- legacy code is removed after stabilization;
- tests identify the reduction and reference model.

---

## 20. Benchmark strategy

Benchmarks are regression protection, not the primary purpose of the refactor.

Representative benchmarks should record:

- wall time;
- peak resident memory;
- marker updates per second where available;
- one-chain and multi-chain behavior;
- one-core and multiple-core behavior;
- minimal and rich output;
- initialization/conversion overhead.

Benchmarks should compare:

- pre-migration backend;
- migrated backend;
- future new implementation when one is required.

No claim of performance improvement should be made without measurement.

---

## 21. Revised implementation phases

## Phase 0 — repository audit

Completed in:

```text
docs/dev/blr_framework_phase0_audit.md
```

The audit remains a factual record and should not be rewritten.

## Phase 1 — contracts, boundaries, and regression baselines

**Status: complete.** Typed specifications, typed results, CSR ownership, and
the binding round trip are implemented and validated.

Implement:

- typed `DataSpec`;
- typed `ModelSpec`;
- typed `McmcControl`;
- typed `OutputSpec`;
- typed result vocabulary;
- language-neutral core rules;
- R binding conversion skeleton;
- CSR operator ownership contract;
- deterministic legacy CSR BayesC fixtures;
- runtime and memory regression baselines.

Do not rewrite the sampler.

Readiness gate:

- specifications compile without Rcpp;
- conversion is tested;
- baseline results are reproducible;
- benchmark variability is documented;
- no public behavior changes.

## Phase 2 — migrate unscheduled CSR BayesC

**Status: complete.** Ordinary unscheduled CSR BayesC executes behind the
typed binding-neutral boundary with exact reference, runtime, and memory
validation.

Refactor the existing efficient kernel into the new boundaries while preserving:

- mathematics;
- marker traversal;
- RNG ordering;
- chain aggregation;
- residual operations;
- performance;
- memory use.

Use typed input/result boundaries and existing or extracted shared infrastructure.

Readiness gate:

- exact old/new comparisons where possible;
- no material runtime or memory regression;
- no Rcpp in migrated reusable core;
- full tests pass.

## Phase 3 — canonicalize and stabilize unscheduled CSR BayesC

**Status: complete.** The migrated implementation is the sole ordinary CSR
BayesC production implementation. The public route always reaches the typed
core, raw-result conversion is centralized, the implementation header is
single-translation-unit guarded, deterministic references remain permanent,
and the legacy ordinary-CSR path is absent. The public schema is unchanged.

Readiness gate:

- one canonical ordinary CSR marker loop;
- no old/new runtime selector or fallback;
- exact reference and reproducibility tests pass;
- shared CSR ownership remains borrowed and immutable;
- runtime and memory validation is recorded.

## Phase 4 — extract proven shared scalar infrastructure

**Status: complete.** Trait-major scalar chain tasks, exact seed resolution,
chain execution status/RNG ownership vocabulary, and retained-iteration timing
are binding-neutral and active in canonical CSR BayesC. The matching current
CSR BayesR uses and the future migration seam are documented, but BayesR
production execution remains unchanged and is not migrated.

Extract only utilities proven to be shared by at least two migrated models:

- chain execution;
- RNG ownership;
- categorical/log-weight utilities;
- scalar variance updates;
- online accumulation;
- result assembly;
- probability/scale interfaces.

Do not generalize speculative future behavior.

Use this infrastructure to prepare migration of unscheduled CSR BayesR without
altering the canonical CSR BayesC hot loop.

## Phase 5 — migrate unscheduled CSR BayesR (complete; canonicalized in Phase 6)

Move the existing efficient unscheduled CSR BayesR implementation behind typed
execution and result boundaries, adopting only the shared scalar
infrastructure validated in Phase 4.

Preserve the specialized marker kernel, component mathematics, RNG ordering,
public schema, speed, and memory use. BayesRC/SBayesRC remains a later task.

Ordinary CSR BayesR now uses a borrowed typed input, one operator-aware
templated execution core, a typed execution result, and one binding-layer raw
converter. The public route and schema are unchanged; six raw and six formatted
references are exact. Block-eigen BayesR remains externally unchanged and is
not itself marked migrated. Phase 6 may canonicalize and remove migration-only
aliases after performance stabilization.

Phase 6 canonicalized this route: one guarded numerical implementation header,
one marker loop, one binding converter, permanent exact fixtures, and a
moderate post-migration timing/RSS baseline. Ordinary CSR BayesR is canonical;
the block-eigen instantiation remains preserved but is not independently
migrated.

## Phase 6 — migrate annotation and scale models

Migrate:

- fixed marker priors;
- `selection_s`;
- group model;
- learned annotation model;
- categorical hierarchy;
- hierarchy × mixture/annotation composition.

## Phase 7 — migrate operators

Migrate or adapt:

- packed BED;
- block eigen;
- dense individual;
- dense summary.

Preserve efficient operator-specific implementations.

## Phase 8 — new multivariate BayesC

Implement a coherent new small-\(T\) model.

Do not port the legacy MT sampler.

## Phase 9 — multivariate mixtures and hierarchy

Add restricted pattern × mixture models and hierarchy.

## Phase 10 — factor-analytic BLR

Implement low-rank covariance and alignment.

## Phase 11 — evidence-factor model

Implement gene-level evidence factors and later biological propagation.

---

## 22. Source organization

A logical target is:

```text
src/
  core/
    blr_spec.h
    blr_result.h
    blr_rng.h
    blr_chain.h
    blr_accumulator.h
    blr_diagnostics.h

  operators/
    csr_operator.h
    bed_operator.h
    block_eigen_operator.h
    dense_operator.h

  scalar/
    bayesc/
    bayesr/
    bayesrc/
    scale/
    annotation/

  multivariate/
    mt_state.h
    mt_covariance.h
    mt_bayesc.h

  factor/
    factor_state.h
    factor_alignment.h

  bindings/
    r/
      spec_from_r.cpp
      result_to_r.cpp
      rcpp_entries.cpp
```

The physical layout may remain flatter if required by the R build system. Logical separation matters more than directory depth.

---

## 23. Work split

### Design and review

Responsible for:

- model definitions;
- covariance contracts;
- normalization;
- extension points;
- reduction tests;
- migration gates;
- review of reports and benchmarks.

### Codex

Responsible for:

- repository-grounded refactoring;
- compilation;
- focused tests;
- package tests;
- regression benchmarks;
- documentation reports;
- scoped deletion when explicitly requested.

### Maintainer

Responsible for:

- approving architecture decisions;
- choosing public switch points;
- committing clean checkpoints;
- approving performance trade-offs;
- deciding when legacy code is removed.

---

## 24. Codex task rules

Each task should:

- have one principal objective;
- name the preserved behavior;
- name untouched model families;
- define statistical and regression gates;
- avoid unrelated cleanup;
- report exact commands and results;
- avoid commit/push unless requested.

Do not combine in one task:

- architecture redesign;
- sampler rewrite;
- new statistical model;
- operator migration;
- performance optimization;
- legacy deletion.

---

## 25. Adopted decisions

- Incremental architecture-first refactor.
- Preserve fast and memory-efficient validated kernels.
- Rewrite only statistically or structurally unsuitable code.
- Legacy implementations are temporary references during migration.
- One canonical implementation after stabilization.
- Language-neutral C++ core.
- Thin R binding and future Python binding.
- Shared infrastructure with specialized hot loops.
- Performance and memory are regression constraints.
- Multivariate BLR will be newly implemented.
- Public schema redesign is deferred until common semantics are stable.

---

## 26. Current migration boundary

Ordinary-CSR SBayesRC is canonical. Its typed borrowed CSR/annotation context,
operator-aware templated core, typed result, and sole ordinary-CSR Rcpp
converter are active. Ordered-probit stick-breaking and alpha behavior are
unchanged; all six permanent raw and formatted references are exact, and a
canonical timing and completed-fit RSS baseline is recorded. Shared scalar
task, seed, and status infrastructure is active. The public route and schema
are unchanged, and migration-only active-path wording has been removed. CSR
BayesC, BayesR, and fixed-prior CSR BayesC are canonical. Fixed-prior BayesC
uses borrowed typed priors, one callable core, one converter, and one preserved
wrapper aggregation path. The existing block-eigen SBayesRC
compile-time instantiation is preserved but is not an independently migrated
typed public adapter. Group CSR BayesC is canonical behind typed execution
boundaries. Learned-annotation CSR BayesC has completed migration behind a
typed borrowed annotation context, active learned policy, callable core, typed
result, one binding converter, and one preserved wrapper aggregation path.
Centered-logistic probabilities, exponential multipliers, proposal, bound, and
update-frequency semantics, public routing, and schema remain unchanged; exact
references and post-migration runtime/completed-fit-RSS checks are active. It
is ready for canonicalization. BED BayesRC, scheduled, and multivariate
backends remain unchanged.

### Phase 9A comparative boundary

Fixed-prior, group, and learned-annotation CSR BayesC production routes were
frozen here before migration. They share binding-neutral borrowed-CSR ownership,
scalar controls, LD-swap controls, result vocabulary, and policy tags. Fixed
marker vectors are immutable inputs; group state has indexed Beta/variance
updates and optional normalization; learned annotations use centered logistic
and exponential links with random-walk Metropolis updates. Nine permanent
raw/formatted reference pairs and per-backend benchmarks protect these
boundaries. Fixed-prior BayesC subsequently completed Phase 9B migration;
group subsequently completed Phase 9D migration; learned annotation completed
Phase 9F migration without changing its public route or policy semantics, and
Phase 9G subsequently established it as canonical.

### Phase 9B3 fixed-prior migration closure

Fixed-prior CSR BayesC now uses the explicit borrowed
`CsrPriorBayesCExecutionContext`, callable `run_csr_prior_bayesc()` numerical
core, typed `CsrPriorBayesCExecutionResult`, and one named binding-layer result
converter. Borrowed immutable `pi_marker` and `vb_multiplier` semantics are
active. Wrapper-level multichain aggregation, the public route, native public
signature, `stblr_raw_v1`, and formatted fit remain unchanged. All frozen
references are exact, runtime and completed-fit RSS were checked, and the
existing unsupported trait-dimension construction remains rejected. The route
completed Phase 9C canonicalization with permanent exact fixtures and canonical
runtime/completed-fit-RSS baselines. Group CSR BayesC is canonical, and
learned-annotation CSR BayesC is migrated behind its typed execution boundary
with exact references retained.

### Phase 9C fixed-prior canonical implementation

Fixed-prior CSR BayesC is canonical. Borrowed immutable `pi_marker` and
`vb_multiplier`, the callable numerical core, typed execution result, sole
binding converter, and sole wrapper-level aggregation path are active. The
public route and schema are unchanged; permanent exact references and runtime
and completed-fit RSS baselines are established; the current trait-dimension
restriction remains protected. Migration scaffolding and active-path wording
have been removed.

### Phase 9D3 group BayesC migration closure

Group CSR BayesC is migrated behind typed execution boundaries. The typed
borrowed zero-based group mapping and `GroupBayesCPolicyView`, callable
`run_csr_group_bayesc()` numerical core, typed execution result, one named
binding-layer result converter, and one wrapper-level multichain aggregation
path are active. Group order, update timing, both normalization modes, RNG,
scheduling, public routing, and public schemas are unchanged. All three raw
and three formatted frozen references are exact; runtime and completed-fit RSS
were checked, and existing unsupported cases remain protected. Phase 9E
subsequently canonicalized this route. Learned-annotation CSR BayesC remains
production unchanged with contracts and references ready, but is not migrated.

### Phase 9E canonical group BayesC implementation

Group CSR BayesC is canonical. Its typed borrowed zero-based group mapping,
active `GroupBayesCPolicyView`, callable numerical core, typed execution result,
sole binding converter, and sole wrapper-level multichain aggregation path are
permanent. Group probability, multiplier, ordering, normalization, global
parameter, RNG and scheduling semantics are preserved. The public route,
native signature, `stblr_raw_v1` and formatted schema are unchanged. Three raw
and three formatted fixtures are permanent and exact; canonical timing and
completed-fit RSS baselines are established; existing unsupported cases remain
protected. Migration scaffolding is removed. Learned-annotation CSR BayesC is
also canonical behind its distinct typed learned-policy boundary.

### Phase 9G canonical learned-annotation BayesC implementation

Learned-annotation CSR BayesC is canonical. Its typed borrowed annotation
design, active `LearnedAnnotationBayesCPolicyView`, callable numerical core,
typed execution result, sole binding converter, and sole wrapper-level
multichain aggregation path are permanent. Centered-logistic probability,
exponential multiplier, coefficient prior/proposal, bound, clipping,
update-frequency, diagnostic, RNG, and scheduling semantics are preserved. The
public route, native signature, `stblr_raw_v1`, and formatted schema are
unchanged. Three raw and three formatted fixtures are permanent and exact;
canonical runtime and completed-fit RSS baselines are established; existing
unsupported cases remain protected; migration scaffolding is removed.

### Phase 10A scheduled ordinary-CSR audit

Scheduled ordinary CSR currently means scheduled BayesC only; scheduled CSR
BayesR is rejected and no scheduled CSR BayesRC/SBayesRC route exists. Phase
10A leaves production execution unchanged while establishing binding-neutral
sweep, skip, candidate, neighbor-wakeup, execution, mutable-state, diagnostic,
and result vocabularies. Three fresh-process raw and formatted scheduled BayesC
references protect dense and nontrivial scheduling configurations.

The audit confirms that production's `static thread_local`
`std::normal_distribution<double>` is worker-thread-owned rather than fit- or
chain-owned. Its cached second variate survives reconstruction/reseeding of the
chain engine, producing same-process call-order dependence and potential
thread-assignment dependence. Scheduled RNG ownership correction is therefore
required before execution migration. All canonical ordinary-CSR models remain
unchanged and scheduled CSR is neither canonical nor migrated.

### Phase 10B scheduled ordinary-CSR RNG correction

Scheduled BayesC now constructs one `ScheduledChainRng` after each final
trait-chain seed is resolved. Its `std::mt19937`, normal distribution, and
uniform distribution live for exactly one chain execution and are never owned
by or shared through an OpenMP worker. Variable-parameter chi-square and gamma
distributions remain locally constructed at their unchanged logical draw sites.

The scheduler, seed formulas, task order, OpenMP static scheduling, marker
traversal, skipped-marker policy, posterior formulas, public route, and schema
are unchanged. Repeated, intervening-fit, fresh/reused-process, different-chain,
and 1/2-core sequences are exact with post-correction references. Scheduled
BayesR remains unsupported and scheduled execution migration has not started.

### Phase 10C scheduled ordinary-CSR execution migration

Phase 10C1 mechanically extracted the corrected deterministic scheduled BayesC
execution body into one guarded implementation header. Phase 10C2 activates the
Phase 10A scheduler contracts through a binding-neutral
`CsrScheduledBayesCExecutionContext`, callable
`run_csr_scheduled_bayesc()`, and typed
`CsrScheduledBayesCExecutionResult`. CSR/statistic/prior/order inputs remain
borrowed immutable; scheduler, sampler, accumulator, and
`ScheduledChainRng` state remain logical-chain-owned. Phase 10C3 closes the
migration with one named binding-layer converter and one native aggregation
path. Corrected references and fit-local, fresh-process, explicit-seed, and
worker-assignment reproducibility remain exact. The public route and schema are
unchanged. Scheduled BayesC migration is complete and ready for
canonicalization; scheduled BayesR remains unsupported.

### Phase 10D scheduled ordinary-CSR BayesC canonicalization

Scheduled ordinary-CSR BayesC is canonical. The corrected fit-bounded,
chain-owned RNG design, Phase 10A scheduler contracts, borrowed typed context,
callable numerical core, typed result, one native aggregation path, and one
named binding converter are the sole production architecture. Phase 10B
corrected fixtures are permanent canonical references; Phase 10A defective
fixtures remain historical audit evidence only. Fit-local and thread-assignment
reproducibility, public routing, `stblr_raw_v1`, and formatted schemas remain
unchanged. Runtime and completed-fit-RSS baselines are established, unsupported
cases are preserved, and migration scaffolding is absent. Scheduled ordinary-
CSR BayesR remains unsupported.

### Phase 11A individual-level and packed-BED backend audit

The remaining public individual-level route is `stblr_bed()`, backed entirely
by SNP-major packed-BED marker execution rather than an in-memory genotype
matrix. BayesC exposes historical sparse and scheduled single-chain native
entries plus the active scheduled multichain route; BayesR uses adaptive
scheduled multichain execution; BayesRC uses unscheduled full sweeps with an
annotation-component prior. Production execution is unchanged.

Packed-BED BayesC scheduling is semantically close to packed-BED BayesR and to
canonical scheduled CSR for sweeps, adaptive null skips, candidates, due
buckets and skipped-marker RNG avoidance, but lacks CSR neighbor wake-up and
has backend-specific decoding and seed ownership. BayesRC does not share the
adaptive scheduler. Shared scheduled infrastructure is therefore only
partially appropriate: control vocabulary may be reused as a subset for BED
BayesC/BayesR after correction, while genotype access and BayesRC execution
remain backend-specific. Active packed-BED BayesC single/multichain samplers
retain worker-owned `static thread_local` normal/uniform distributions and
require fit-local RNG correction before migration. BayesR and BayesRC construct
distributions per logical chain. Deterministic references are established
where valid; migration has not started.

### Phase 11B scheduled packed-BED BayesC RNG correction

Scheduled single-chain and multichain packed-BED BayesC now construct one
`BedScheduledBayesCChainRng` per logical chain. Its `std::mt19937`, normal and
uniform distributions have one-chain lifetimes; no stateful distribution is
static, thread-local, worker-owned or fit-persistent. Seed formulas, logical
draw sites, scheduler transitions, genotype decoding, I/O, public routing and
schemas are unchanged. Post-correction deterministic references and fit-local,
worker-independent reproducibility are active. Numerical-core migration has
not started; packed-BED BayesR and BayesRC remain unchanged.

### Phase 11C1 packed-BED BayesC mechanical extraction

The corrected 417-line logical-chain numerical body behind the public
multichain packed-BED BayesC route is mechanically located in
`blr_bed_scheduled_bayesc_core_impl.h`. The guarded implementation header holds
the sole active public-path MCMC/marker loop and chain RNG construction. BED
loading/decoding, task orchestration, native aggregation and existing R result
construction remain in the selected binding source. The experimental
single-chain route is unchanged. Typed callable execution boundaries are not
yet active; migration remains in progress.

### Phase 11C2 packed-BED BayesC typed per-chain boundary

The corrected public multichain route now constructs a binding-neutral
`BedScheduledBayesCChainExecutionContext` for each logical trait-chain and calls
`run_bed_scheduled_bayesc_chain()`, which returns a typed
`BedScheduledBayesCChainExecutionResult`. Fit-owned decoded packed genotype
storage, marker maps, phenotypes and priors are borrowed immutable. The Phase
10A sweep, null-skip and candidate contracts are reused with their identical
BED semantics; neighbor wake-up is not activated. Logical-chain RNG ownership
and Phase 11B corrected references remain exact. OpenMP task dispatch,
cross-chain aggregation and existing result conversion remain in the adapter;
migration remains in progress. The experimental route, BayesR and BayesRC are
unchanged.

### Phase 11C3 packed-BED BayesC migration closure

The public multichain packed-BED BayesC route now feeds its typed per-chain
results through one binding-neutral `aggregate_bed_scheduled_bayesc_results()`
path and one typed `BedScheduledBayesCExecutionResult`. One named binding
converter, `stblr_bed_scheduled_bayesc_result_to_raw()`, constructs the
unchanged `stblr_raw_v1`. Fit-local decoding and OpenMP task dispatch remain in
the adapter; the per-chain numerical core remains binding neutral and borrows
immutable packed genotype storage. Logical-chain RNG ownership, scheduler
semantics, Phase 11B corrected references, and fit-local/thread-independent
reproducibility remain exact. The public route is migrated and ready for
canonicalization; the experimental route, BayesR, and BayesRC are unchanged.

### Phase 11D public scheduled packed-BED BayesC canonicalization

The public multichain packed-BED BayesC route is canonical. Its fit-owned
decoded packed genotype is borrowed immutable through a typed per-chain
context; one callable chain core, one typed chain result, one typed aggregate
result, one native aggregation path, and one named binding converter remain.
Logical-chain-owned RNG state, scheduler semantics, public routing and
`stblr_raw_v1` are unchanged. Phase 11B corrected fixtures are permanent
canonical references, and Phase 11C3/11D establish runtime, completed-fit RSS,
BED-size and page-cache-qualified I/O baselines. Migration scaffolding is
absent. The experimental scheduled single-chain and sparse BayesC routes,
packed-BED BayesR, and packed-BED BayesRC remain unchanged and noncanonical.

### Phase 12A cross-cutting hardening

Current-status metadata is synchronized; maintained Markdown rejects control
bytes; structural source counts distinguish zero, one, and duplicates; and
`full_sweep_every=0` is consistently accepted as full sweep every iteration
for scheduled routes already implementing that behavior. Posterior aggregation
and converter responsibilities are explicit, unconditional scheduled-CSR core
console writes are removed, canonical commented snapshots are removed, fast
and extended CI are established, and peak-RSS tooling distinguishes sampled
peak from completed-fit RSS. Canonical numerical references remain unchanged.
