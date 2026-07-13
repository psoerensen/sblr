# Unified BLR Framework: Architecture and Implementation Plan

**Status:** Draft architectural plan  
**Date:** 2026-07-13  
**Target location:** `docs/dev/blr_framework_implementation_plan.md`

## 1. Purpose

This document defines the implementation strategy for a maintainable and extensible Bayesian linear regression framework in `sblr`.

The framework must support, without duplicating complete samplers for every combination:

- single-trait Bayesian linear regression;
- multivariate Bayesian linear regression;
- BayesC, BayesR, and BayesRC state models;
- fixed, learned, grouped, and hierarchical annotation priors;
- MAF-dependent effect-size scaling through `selection_s`;
- dense, sparse-CSR, block-eigen, and packed-BED likelihood representations;
- future factor-analytic multivariate models;
- a separate annotation-informed evidence-factor model for genes, pathways, and target discovery.

The central architectural decision is:

> A model is a composition of independently testable policies operating inside a shared chain runner and likelihood operator, rather than a separate end-to-end backend for every model/data combination.

This plan should be treated as the authoritative design reference before major multivariate refactoring or implementation of hierarchical and factor models.

---

## 2. Design goals

### 2.1 Statistical correctness

Every model component must have:

- an explicit mathematical definition;
- an explicit parameterization and scale convention;
- an independently testable implementation;
- reduction tests to a simpler validated model;
- deterministic fixed-seed behavior;
- a documented interpretation of returned quantities.

No heuristic covariance update, shrinkage rule, normalization, or approximation should be hidden behind a generic function name.

### 2.2 Maintainability

The framework should minimize duplicated code for:

- residual updates;
- marker scores;
- marker-state sampling;
- MCMC scheduling;
- chain seeds and random-number generation;
- burn-in and thinning;
- chain aggregation;
- raw result construction;
- output formatting;
- annotation preprocessing;
- marker alignment.

### 2.3 Extensibility

Adding a new prior or data representation should usually require implementing one component rather than a new complete sampler.

Examples:

- hierarchical BayesC should primarily add a marker-scale policy;
- BayesRC should primarily add an annotation-dependent state-probability policy;
- multivariate BayesC should add a trait-pattern state and trait-covariance policy;
- a new LD representation should add a likelihood operator;
- factor-analytic multivariate BLR should add a trait-covariance representation and factor state;
- the evidence-factor model should reuse the factor core but have its own likelihood frontend.

### 2.4 Reproducibility

All stochastic state must belong to a chain.

The framework must guarantee:

- identical results for repeated calls with identical inputs and seeds;
- identical results under one and multiple OpenMP threads when chain scheduling changes;
- no cross-call random-distribution contamination;
- explicit and documented chain seed rules;
- no hidden use of global R, Armadillo, or static thread-local random-number state inside native chain execution.

### 2.5 Language-neutral native core

The native statistical engine should be usable through more than one language binding.

At this stage this means:

- the statistical core uses ordinary C++ types and exceptions;
- Rcpp remains in thin R-facing adapters;
- chain RNG does not depend on R, Armadillo global RNG, or language-runtime state;
- model specifications and result objects have typed C++ representations;
- core dimensions and semantics are independent of R list or matrix formatting;
- disk-backed LD and genotype formats remain language neutral.

A Python interface is not part of the current implementation scope. The framework should only preserve the architectural option of adding a future `pybind11` adapter over the same C++ engine.

The project does not need to promise a stable public C++ ABI during the refactor. The internal C++ API may evolve until the statistical framework is stable.

### 2.6 Backward-compatible migration

Existing public functions may remain as user-facing convenience wrappers:

- `sblr()`;
- `stblr_csr()`;
- `stblr_bed()`;
- `stblr_csr_annot()` and explicit annotation wrappers.

They should gradually become constructors of one resolved internal model specification.

Legacy native backends should remain available as validation references until the new implementation passes exact or tolerance-based reduction tests.

---

## 3. Current repository context

### 3.1 Single-trait framework

The single-trait side currently provides important reusable infrastructure:

- sparse-CSR summary-statistics backends;
- packed PLINK BED backends;
- block-eigen LD operators;
- BayesC, BayesR, SBayesRC, and individual-level BayesRC;
- fixed, learned, and grouped annotation-aware BayesC models;
- MAF-dependent `selection_s`;
- LD-swap support for selected CSR models;
- multiple-chain support;
- a named raw schema and common formatter;
- strict fixed-seed and thread-count validation for recent BED BayesRC work.

This is the strongest current foundation and should supply the canonical conventions for:

- chain ownership;
- RNG ownership;
- annotation preprocessing;
- raw result naming;
- formatted outputs;
- reduction testing.

### 3.2 Multivariate framework

The current multivariate code is older and structurally separate. It includes multiple parallel native implementations and an `sblr()` interface that handles model setup, native dispatch, and positional output formatting directly.

Before it becomes the basis of new model development, the multivariate code requires a statistical and architectural audit covering:

- trait-pattern state definitions;
- covariance updates;
- residual covariance assumptions;
- exact use of off-diagonal precision terms;
- RNG ownership;
- data-design assumptions;
- output semantics;
- duplicate and experimental implementations.

The existing multivariate code should initially be treated as a reference implementation, not as the final shared core.

### 3.3 New use cases

Two new model families motivate this design.

#### Hierarchical BayesC and hierarchical BayesRC

Marker-effect variances depend on multiple annotation layers, for example:

\[
\operatorname{Var}(\beta_j \mid d_j=1)
=
v_b
\prod_{h=1}^{H}
\theta_{h,g_h(j)}.
\]

A hybrid BayesRC model additionally uses annotation-dependent mixture membership and component scaling:

\[
\operatorname{Var}(\beta_j \mid c_j=k)
=
v_b \gamma_k q_j.
\]

#### Annotation-informed latent factor models

A gene- or SNP-by-trait evidence matrix is decomposed as:

\[
Y = U V^\top + E,
\]

with sparse annotation-informed row-factor loadings, trait-factor loadings, and later propagation from SNPs to genes and pathways.

This is related to factor-analytic multivariate BLR but is not the same likelihood. It should share a factor core while retaining a separate public interface and result class.

---

## 4. Unified statistical representation for regression models

Let:

\[
\mathbf b_j = (b_{j1},\ldots,b_{jT})^\top
\]

denote marker \(j\)'s effects on \(T\) traits.

Define a discrete marker state:

\[
z_j = (r_j,k_j),
\]

where:

- \(r_j\) is a trait-activity pattern;
- \(k_j\) is an effect-size component.

The general prior is:

\[
P(z_j \mid A_j,\Theta)
\]

and

\[
\mathbf b_j \mid z_j=(r,k)
\sim
N_T\left(
0,\;
q_j D_r \Sigma_k D_r
\right).
\]

Here:

- \(A_j\) is the marker annotation vector;
- \(D_r\) is a diagonal matrix selecting active traits;
- \(\Sigma_k\) is the component-specific trait covariance;
- \(q_j\) is a marker-specific scale multiplier.

This representation separates four concerns:

1. which traits are active;
2. which effect-size component is active;
3. how annotations affect state probabilities;
4. how annotations, MAF, groups, or hierarchies affect effect-size variance.

### 4.1 Model mappings

| Model | Marker state | State probabilities | Marker scale | Trait covariance |
|---|---|---|---|---|
| ST BayesC | null/active | global Beta/Bernoulli | \(1\) | scalar \(v_b\) |
| Hierarchical ST BayesC | null/active | global Beta/Bernoulli | \(q_j\) | scalar \(v_b\) |
| MT BayesC | trait pattern | global pattern probabilities | \(1\) | \(D_r\Sigma_BD_r\) |
| Hierarchical MT BayesC | trait pattern | global pattern probabilities | \(q_j\) | \(D_r\Sigma_BD_r\) |
| ST BayesR | mixture component | global Dirichlet | \(\gamma_k\) | scalar \(v_b\) |
| ST BayesRC | mixture component | annotation stick-breaking | \(\gamma_k\) | scalar \(v_b\) |
| Hierarchical BayesRC | mixture component | annotation stick-breaking | \(\gamma_k q_j\) | scalar \(v_b\) |
| MT BayesR/RC | trait pattern × component | global or annotation-informed | \(\gamma_k q_j\) | \(D_r\Sigma_BD_r\) |
| Factor MT-BLR | sparse factor state | factor-specific | factor scale | \(\Lambda\Lambda^\top+\Psi\) |

---

## 5. Architectural principles

### 5.1 Separate likelihood from prior

The likelihood operator must know how to:

- obtain marker diagonal information;
- calculate a marker score or partial residual;
- apply an effect change;
- rebuild the residual;
- expose any LD-neighbor information needed by optional proposals.

It must not know whether the model is BayesC, BayesR, BayesRC, hierarchical, or multivariate.

### 5.2 Separate discrete state from effect-scale priors

Trait-sharing patterns, mixture components, and annotation probabilities must not be embedded in the same class that computes hierarchy multipliers.

### 5.3 Separate covariance definitions from covariance updates

A covariance object defines the effective prior covariance.

A covariance update policy defines how its parameters are sampled.

For example, a full inverse-Wishart update and a separation strategy for variances and correlations are different statistical models and must have separate identifiers and tests.

### 5.4 Resolve all user input before C++

R should:

- align markers and annotations;
- validate dimensions and names;
- choose allowed marker states;
- initialize priors;
- choose operator and policy identifiers;
- construct an immutable resolved model specification.

C++ should execute the resolved model, not interpret a large ambiguous public argument list.

### 5.5 Keep the statistical core independent of language bindings

The model core must not construct or consume binding-specific objects such as:

- `Rcpp::List`;
- `Rcpp::NumericMatrix`;
- Python dictionaries or NumPy wrapper objects.

The preferred boundary is:

```text
R arguments
  -> R-side validation and preprocessing
  -> Rcpp adapter
  -> typed C++ model/data/control specification
  -> C++ statistical engine
  -> typed C++ result
  -> Rcpp result adapter
  -> R fit object
```

A future Python path should use the same middle layers:

```text
Python arguments
  -> Python validation
  -> pybind11 adapter
  -> the same typed C++ specification
  -> the same C++ statistical engine
  -> the same typed C++ result
  -> Python result adapter
```

The bindings must not contain separate statistical implementations.

### 5.6 Keep model-specific outputs namespaced

Common outputs should be standardized.

Model-specific outputs should be placed in namespaces such as:

- `raw$state`;
- `raw$component`;
- `raw$pattern`;
- `raw$annotation`;
- `raw$hierarchy`;
- `raw$factor`.

---

## 6. Framework layers

## 6.1 Data and likelihood operators

Conceptual interface:

```cpp
struct Operator {
  int n_markers() const;
  int n_traits() const;

  double diag(int marker, int trait) const;

  // Obtain the score/partial-residual contribution needed by a marker update.
  void score(int marker,
             const EffectState& effect,
             const ResidualState& residual,
             arma::vec& out) const;

  // Apply a marker effect difference to the residual state.
  void apply_effect_change(int marker,
                           const arma::vec& delta,
                           ResidualState& residual) const;

  // Reconstruct residual state from current effects.
  void rebuild(const EffectState& effect,
               ResidualState& residual) const;
};
```

Planned implementations:

- `DenseIndividualOperator`;
- `PackedBedOperator`;
- `DenseSummaryOperator`;
- `CsrOperator`;
- `BlockEigenOperator`.

The exact interface may differ, but the required operations and semantics must be explicit.

### Multivariate data contract

Every multivariate data specification must declare one of:

- `shared_design`: same individuals and genotype design across traits;
- `trait_specific_design`: different samples and/or genotype designs;
- `summary_shared_ld`: summary statistics with a shared LD operator;
- `evidence_matrix`: factor-model input, not a regression likelihood.

A full residual covariance matrix is only supported when the available sufficient statistics identify it.

Unsupported combinations must error during R-side specification resolution.

---

## 6.2 State-space policy

The state-space policy defines possible discrete states and their active trait/component structure.

Initial policies:

- `BinaryState`;
- `TraitPatternState`;
- `MixtureState`;
- `PatternMixtureState`.

A trait-pattern state consists of user-supplied binary trait masks.

Full \(2^T\) enumeration must not be the default for moderate or large \(T\). The public interface should support:

```r
trait_sharing = "independent"
trait_sharing = "shared"
trait_sharing = "patterns"
trait_sharing = "factor"
```

For `trait_sharing = "patterns"`, the allowed patterns are explicit and validated.

---

## 6.3 State-probability policy

The probability policy provides:

\[
P(z_j \mid A_j,\Theta).
\]

Initial policies:

- `GlobalBinaryProbability`;
- `GlobalDirichletProbability`;
- `GroupProbability`;
- `AnnotationLogitProbability`;
- `AnnotationProbitStickProbability`;
- future `AnnotationPatternProbability`.

The current BayesRC probit stick-breaking machinery should be generalized only at the utility level. A factor model's independent Bernoulli memberships are not the same as mutually exclusive BayesRC stick-breaking components.

---

## 6.4 Marker-scale policy

The scale policy provides \(q_j\), independently of state probabilities.

Initial policies:

- `UnitScale`;
- `FixedMarkerScale`;
- `MafSelectionScale`;
- `GroupScale`;
- `CategoricalHierarchyScale`;
- future `LogAdditiveAnnotationScale`;
- `CompositeScale`.

A composite scale may be:

\[
q_j =
q_j^{\mathrm{MAF}}(S)
\prod_{h=1}^{H}
\theta_{h,g_h(j)}.
\]

This supports the combination of `selection_s`, time-window effects, event-type effects, and other categorical variance hierarchies.

---

## 6.5 Trait-covariance policy

The trait-covariance policy provides \(\Sigma_k\).

Initial policies:

- `ScalarTraitVariance`;
- `DiagonalTraitCovariance`;
- `FullTraitCovariance`;
- future `FactorAnalyticTraitCovariance`.

For a trait pattern \(r\), the active covariance is:

\[
D_r \Sigma_k D_r.
\]

Small-\(T\) models may use explicit patterns and full covariance.

Moderate- or large-\(T\) models should use factor-analytic structure rather than enumerating all trait patterns.

---

## 6.6 Covariance-update policy

Initial policies:

- `FixedCovarianceUpdate`;
- `DiagonalInverseChiSquareUpdate`;
- `FullInverseWishartUpdate`;
- future `SeparationCovarianceUpdate`;
- future `FactorAnalyticCovarianceUpdate`.

Each policy must document:

- its prior;
- its conditional posterior or proposal;
- its valid data-design assumptions;
- its positive-definiteness guarantees;
- its effect on RNG ordering.

---

## 6.7 Chain runner and RNG ownership

A shared chain runner should own:

- chain seed;
- `std::mt19937`;
- stateful standard distributions;
- marker order;
- iteration loop;
- burn-in and thinning;
- state accumulation;
- residual rebuild schedule;
- chain diagnostics;
- failure capture.

No stateful random distribution may be:

- global;
- static;
- shared across chains;
- silently provided by Armadillo;
- silently provided by R's RNG inside an OpenMP chain;
- silently provided by a future Python/NumPy runtime.

The statistical engine should use chain-owned standard C++ RNG objects. This ensures that an R and a future Python binding can execute the same resolved model with the same seed semantics.

Parallelism should primarily distribute independent trait-chain jobs.

The statistical trajectory of a chain must not depend on which OpenMP thread executes it.

---

## 6.8 Aggregation and raw output

The chain aggregator must:

- preserve exact no-chain `NULL` conventions;
- combine posterior means;
- combine posterior traces using documented rules;
- retain compact chain summaries when requested;
- preserve trait and marker dimensions;
- validate identities such as component probabilities summing to one.

A future namespaced schema may be introduced as `blr_raw_v1` or `stblr_raw_v2`.

Proposed structure:

```r
raw$schema
raw$meta

raw$marker
raw$trace

raw$trait
raw$residual

raw$state
raw$component
raw$pattern

raw$scale
raw$hierarchy
raw$annotation

raw$factor
raw$diagnostics
raw$chains
```

Existing `stblr_raw_v1` should remain supported during migration.

## 6.9 Binding-neutral C++ specifications and results

The C++ engine should receive typed structures rather than R lists.

Illustrative internal structures:

```cpp
struct McmcControl {
  int nit;
  int nburn;
  int nthin;
  int nchains;
  int ncores;
  std::uint64_t seed;
  bool keep_chains;
};

struct ModelSpec {
  StateSpaceSpec state;
  ProbabilitySpec probability;
  ScaleSpec scale;
  TraitCovarianceSpec trait_covariance;
  ResidualCovarianceSpec residual_covariance;
};

struct BlrResult {
  MarkerResult marker;
  TraceResult trace;
  TraitResult trait;
  ResidualResult residual;
  StateResult state;
  ComponentResult component;
  PatternResult pattern;
  AnnotationResult annotation;
  HierarchyResult hierarchy;
  DiagnosticsResult diagnostics;
  std::vector<ChainResult> chains;
};
```

The exact structures may evolve. Their purpose is to ensure that:

- the sampler does not construct `Rcpp::List` objects;
- the R adapter converts typed results to the current R schema;
- a future Python adapter can expose the same numerical result;
- core output dimensions do not depend on R's dropping or naming rules.

Optional model-specific fields may use typed optionals, variants, or empty structures internally. Binding adapters should translate absence to the language-appropriate representation, such as R `NULL` or Python `None`.

## 6.10 Errors, progress, and logging

The statistical core should:

- throw standard C++ exceptions;
- avoid direct `Rcpp::stop()` calls;
- avoid direct `Rcpp::Rcout` calls;
- avoid direct Python logging in future;
- optionally report structured progress through a generic callback or diagnostics collector.

Binding adapters translate:

- C++ exceptions to R errors;
- C++ exceptions to Python exceptions in a future binding;
- structured progress to the language-specific user interface.

## 6.11 Data ownership and storage contracts

The core must clearly document:

- owned versus borrowed input data;
- read-only versus mutable memory;
- array orientation;
- row-major versus column-major conversions where relevant;
- lifetime of views passed from a binding;
- disk-backed operator ownership;
- thread-safety of shared read-only operators.

Large persisted data should use language-neutral formats. New core infrastructure should not require R serialization or Python pickle for LD, genotype, or model input resources.

The canonical internal dimension conventions are:

```text
marker effects:                 markers × traits
traces:                         samples × traits
trait covariance:               traits × traits
component probabilities:        markers × components × traits
pattern probabilities:          markers × patterns
annotation effects:             annotations × parameters × traits
hierarchy multipliers:          groups × traits within each layer
factor row loadings:            rows × factors
factor trait loadings:          traits × factors
```

Binding adapters may add names and classes, but they should not transpose or reinterpret these semantics silently.

---

## 7. Resolved R model specification

The long-term public construction should separate:

- data;
- model;
- MCMC control;
- output control.

Illustrative internal API:

```r
data_spec <- blr_data_spec(
  representation = "csr",
  stats = stats,
  ld_prefix = ld_prefix,
  design = "summary_shared_ld"
)

model_spec <- blr_model_spec(
  state = trait_pattern_state(
    family = "bayesc",
    patterns = patterns
  ),
  probability = global_pattern_probability(),
  scale = hierarchical_scale(
    window = window_id,
    event_type = event_type_id
  ),
  trait_covariance = full_trait_covariance(),
  residual_covariance = diagonal_residual_covariance()
)

control <- blr_mcmc_control(
  nit = 1000,
  nburn = 100,
  nthin = 1,
  nchains = 4,
  seed = 1
)

fit <- blr_fit(data_spec, model_spec, control)
```

A hierarchical BayesRC specification could be:

```r
model_spec <- blr_model_spec(
  state = mixture_state(
    gamma = c(0, 0.01, 0.1, 1)
  ),
  probability = annotation_stick_breaking(annotation),
  scale = hierarchical_scale(
    window = window_id,
    event_type = event_type_id
  ),
  trait_covariance = scalar_trait_variance()
)
```

These constructors may initially be internal. Existing public wrappers should resolve to the same structures.

---

## 8. Hierarchical marker-scale model

## 8.1 Initial categorical hierarchy

For layer \(h\), marker \(j\) belongs to group \(g_h(j)\).

Define:

\[
q_j =
\prod_{h=1}^{H}
\theta_{h,g_h(j)}.
\]

The effective active-effect prior is:

\[
b_j \mid d_j=1
\sim
N(0,v_b q_j)
\]

for a single trait, or:

\[
\mathbf b_j \mid r_j
\sim
N(0,q_jD_r\Sigma_BD_r)
\]

for the multivariate model.

The first version should support mutually exclusive categorical membership within each layer.

Example:

```r
hierarchy <- list(
  window = window_id,
  event_type = event_type,
  clinical_domain = domain_id
)
```

Each vector has one entry per marker.

## 8.2 Conditional multiplier update

For layer \(h\), group \(g\), define:

\[
q_j^{(-h)}
=
\prod_{\ell\ne h}
\theta_{\ell,g_\ell(j)}.
\]

For a single-trait active marker, a conjugate scaled inverse-chi-square-style update uses:

\[
S_{hg}
=
\sum_{\substack{j:g_h(j)=g\\d_j=1}}
\frac{b_j^2}{v_b q_j^{(-h)}}.
\]

The exact degrees of freedom and scale parameterization must follow one documented convention and match existing variance-prior conventions where possible.

For multivariate effects, the sufficient statistic must be defined from the active-trait covariance model, for example through a quadratic form involving the inverse active covariance.

## 8.3 Identifiability and normalization

The decomposition:

\[
v_b\prod_h\theta_{h,g_h(j)}
\]

is not identifiable without a scale constraint.

After updating layer \(h\), compute the marker-count-weighted geometric mean:

\[
c_h =
\exp\left[
\frac{\sum_g n_{hg}\log\theta_{hg}}
     {\sum_g n_{hg}}
\right].
\]

Then apply:

\[
\theta_{hg}\leftarrow\theta_{hg}/c_h,
\qquad
v_b\leftarrow v_b c_h.
\]

This transformation leaves every effective marker variance unchanged:

\[
v_bq_j.
\]

The scale policy must test this identity directly.

## 8.4 Overlapping annotations

Arbitrary overlaps should not be part of the first categorical implementation.

A later log-additive model may use:

\[
\log q_j = a_0 + A_j a,
\]

with shrinkage priors and Metropolis-Hastings, slice sampling, or another validated update.

This should be a separate `LogAdditiveAnnotationScale` rather than overloading the categorical hierarchy policy.

---

## 9. Hierarchical BayesRC hybrid

The first hybrid should use a marker-level hierarchy shared across non-null components:

\[
P(c_j=k)
=
\pi_{jk}(A_j),
\]

\[
b_j\mid c_j=k
\sim
N(0,v_b\gamma_k q_j).
\]

This separates:

- annotation effects on component membership;
- hierarchy effects on effect-size scale.

Initial reductions:

- \(q_j=1\) gives ordinary BayesRC;
- fixed global component probabilities give hierarchical BayesR;
- one non-null component gives hierarchical BayesC;
- annotation coefficients fixed at their intercept-only values give the matched global mixture model.

Component-specific hierarchy parameters:

\[
q_{jk}
=
\prod_h\theta_{h,g_h(j),k}
\]

should be deferred because they are more weakly identified and confounded with \(\gamma_k\).

---

## 10. Multivariate BLR strategy

## 10.1 Canonical first model

The first rebuilt multivariate model should be a small-\(T\), shared-design BayesC reference model with:

- explicit user-supplied trait patterns;
- one mathematically defined covariance prior;
- one residual covariance contract;
- chain-local RNG;
- named raw output;
- exact reduction to single-trait BayesC at \(T=1\).

The new framework should not initially migrate every legacy multivariate algorithm.

## 10.2 Trait-pattern scaling

The number of all binary trait patterns is \(2^T\).

Therefore:

- full enumeration may be supported only for small \(T\);
- `shared` and `independent` modes should have compact predefined states;
- `patterns` mode should require an explicit restricted set;
- factor models should handle larger trait collections.

## 10.3 Residual covariance

A full residual covariance is only valid when the likelihood and sufficient statistics identify cross-trait residual products.

The model specification must distinguish:

- same individuals and design across traits;
- overlapping samples with known cross-trait information;
- independent summary-statistic traits;
- evidence matrices.

Defaulting to a full residual covariance solely because \(T>1\) is not acceptable.

---

## 11. Factor framework

## 11.1 Shared factor core

A reusable factor core should eventually support:

- row-factor loadings;
- trait-factor loadings;
- sparse factor memberships;
- annotation-informed membership probabilities;
- factor-specific variance parameters;
- factor normalization;
- sign conventions;
- factor ordering;
- signed-permutation alignment across chains;
- posterior factor summaries.

## 11.2 Two separate likelihood frontends

### Factor-analytic MT-BLR

A regression model such as:

\[
Y = X F\Lambda^\top + E.
\]

Marker-to-factor effects are inferred directly from individual-level or summary-statistics regression data.

### Evidence-factor model

An evidence matrix model such as:

\[
Y^{(e)} = U V^\top + E^{(e)}.
\]

The input is a gene-by-trait or SNP-by-trait evidence matrix.

These models may share factor priors and alignment code but must not share a misleading common likelihood interface.

## 11.3 Recommended first evidence-factor implementation

Begin with a gene-by-trait matrix rather than genome-wide SNP-level evidence.

Advantages:

- avoids immediate LD modelling;
- directly targets gene and mechanism discovery;
- makes benchmarking easier;
- separates SNP-to-gene uncertainty from factor inference;
- allows a clean probabilistic factor model to be validated first.

Use signed approximately Gaussian evidence and, where available, observation variances:

\[
Y_{gt}
\sim
N\left(
z_g^\top v_t,
s_{gt}^2+\psi_t^2
\right).
\]

Raw PIPs should not be treated as ordinary Gaussian observations without an explicit transformed or alternative observation model.

## 11.4 Propagation

SNP-to-gene and gene-to-pathway propagation should initially be post-processing on posterior samples:

\[
Z^{(s)} = S^\top U^{(s)},
\]

\[
H^{(s)} = C^\top Z^{(s)}.
\]

This preserves uncertainty and prevents the first factor sampler from becoming a fully joint biological mapping model.

---

## 12. Proposed native and binding layout

The logical target organization is:

```text
src/
  core/
    blr_model_spec.h
    blr_state.h
    blr_chain_state.h
    blr_rng.h
    blr_result.h

    blr_state_space.h
    blr_probability_policy.h
    blr_scale_policy.h
    blr_trait_covariance.h
    blr_covariance_updates.h

    blr_marker_update.h
    blr_chain_runner.h
    blr_aggregation.h

    operators/
      dense_operator.h
      csr_operator.h
      bed_operator.h
      block_eigen_operator.h

    factor/
      factor_core.h
      factor_priors.h
      factor_alignment.h

  bindings/
    r/
      blr_spec_from_r.cpp
      blr_result_to_r.cpp
      blr_rcpp_entry.cpp

    python/                  # future; do not create now
      blr_pybind_entry.cpp
      blr_result_to_python.cpp
```

Because R package build systems may make deeply nested compilation inconvenient, the physical file layout may initially remain flatter. The required constraint is logical separation:

- core headers and implementation must not depend on Rcpp;
- Rcpp entry points and conversions must remain thin;
- no Python directory or dependency is required now;
- future Python support must not require reimplementing the statistical sampler.

This is a design target, not a requirement to create all files at once.

Extraction should occur incrementally and be protected by exact regression tests.

---

## 13. Migration strategy

### 13.1 Preserve validated behavior

Do not rewrite all existing single-trait samplers immediately.

First extract shared abstractions from validated paths without changing posterior output.

### 13.2 Keep reference backends

Legacy backends should remain until the new implementation passes:

- fixed-seed comparisons;
- reduction tests;
- schema comparisons;
- simulation recovery checks.

### 13.3 Route public wrappers gradually

A public wrapper may initially choose:

- legacy backend for established models;
- new framework backend for newly migrated models.

Once equivalence is established, the wrapper may route both to the new engine.

### 13.4 Remove code only after evidence

A legacy implementation is removable only when:

- all supported behavior exists in the new path;
- exact or documented-tolerance equivalence passes;
- documentation and examples have migrated;
- package-wide tests pass;
- a deprecation/migration note exists.

---

## 14. Implementation phases

## Phase 0 — architecture and audit

Deliverables:

- this implementation plan;
- model capability matrix;
- reduction-test matrix;
- audit-only Codex report of multivariate and shared infrastructure.

No sampler changes.

Readiness gate:

- model decomposition approved;
- current multivariate assumptions documented;
- canonical first covariance model selected;
- migration targets identified.

## Phase 1 — resolved model specification and binding boundary

Implement internal R constructors and validators for:

- data specification;
- state-space specification;
- probability specification;
- scale specification;
- covariance specification;
- MCMC control.

Existing public outputs remain unchanged.

Define the corresponding typed internal C++ specification and result contracts without adding a Python binding.

Readiness gate:

- existing wrappers can construct and print resolved specifications;
- Rcpp adapters convert to typed C++ structures;
- core specification/result headers do not include Rcpp;
- no Python dependency is added;
- no native statistical behavior changes;
- validation tests pass.

## Phase 2 — shared chain/RNG infrastructure

Extract or standardize:

- chain seed rules;
- chain-local random distributions;
- burn-in/thinning;
- chain retention;
- aggregation;
- failure reporting.

Readiness gate:

- exact fixed-seed equivalence for migrated backends;
- one-core/multi-core equality;
- no cross-call contamination.

## Phase 3 — canonical multivariate BayesC

Implement one reference backend with:

- shared design;
- explicit trait patterns;
- canonical covariance update;
- named raw output;
- \(T=1\) reduction.

Readiness gate:

- statistical audit complete;
- exact/tolerance-based reductions pass;
- simulated covariance and sharing recovery pass.

## Phase 4 — operator integration

Run the same validated marker prior through:

- dense;
- CSR;
- later packed BED.

Readiness gate:

- operator equivalence tests pass on exact fixtures;
- prior implementation is not duplicated.

## Phase 5 — categorical hierarchy scale

Implement:

- multiple categorical layers;
- conjugate updates;
- weighted-geometric normalization;
- hierarchy output and diagnostics.

Start in CSR.

Readiness gate:

- unit hierarchy reduces to BayesC;
- one layer reduces to group BayesC;
- fixed hierarchy reduces to fixed marker-scale BayesC;
- multivariate \(T=1\) reduction passes.

## Phase 6 — mixture and annotation composition

Add:

- hierarchical BayesR;
- hierarchical BayesRC;
- later multivariate mixture states.

Readiness gate:

- BayesRC and BayesR reductions pass;
- component and hierarchy identities pass;
- annotation and scale effects recover in simulation.

## Phase 7 — factor-analytic multivariate BLR

Implement:

- factor covariance;
- sparse marker-factor effects if required;
- normalization;
- chain factor alignment.

Readiness gate:

- low-rank covariance recovery;
- permutation/sign alignment;
- reduction to simpler covariance cases.

## Phase 8 — gene-level evidence-factor model

Implement a separate public model family:

- gene-by-trait evidence likelihood;
- sparse factor loadings;
- annotation-informed membership;
- posterior factor alignment;
- mechanism-informed outputs.

Readiness gate:

- simulated factor recovery;
- observation-uncertainty calibration;
- benchmark against simpler factor and target-prioritization methods.

## Phase 9 — SNP-to-gene/pathway propagation

Implement posterior-sample propagation with sparse mappings.

Readiness gate:

- deterministic mapping identities;
- uncertainty propagation tests;
- normalization options documented.

## Phase 10 — SNP-level evidence-factor extensions

Only after the gene-level model demonstrates value:

- fine-mapped SNP evidence;
- LD-aware residual modelling if necessary;
- uncertain SNP-to-gene mapping;
- multi-view molecular anchors.

---

## 15. Work split

### Architecture and statistical design

Primary responsibility: design/review work outside Codex.

Tasks:

- define priors and likelihood assumptions;
- select canonical covariance models;
- define identifiability constraints;
- define output semantics;
- design reduction and simulation tests;
- review Codex evidence and diffs;
- approve readiness.

### Repository implementation

Primary responsibility: Codex.

Tasks:

- inspect actual repository code;
- implement bounded architectural steps;
- compile native code;
- write tests;
- run focused and package-wide validation;
- document exact files and behavior changed;
- report unresolved blockers.

### Checkpoint decisions

Primary responsibility: package maintainer.

Tasks:

- approve model and API decisions;
- preserve clean Git checkpoints;
- commit validated phases;
- reject changes that weaken statistical invariants;
- decide when legacy code can be deprecated.

---

## 16. Codex task rules

Each Codex task should:

- have one principal architectural goal;
- state what must not change;
- name required reduction tests;
- require full reporting;
- prohibit unrelated refactoring;
- prohibit commit and push unless explicitly requested.

Examples of appropriately sized tasks:

- audit current multivariate BLR without code changes;
- extract chain/RNG infrastructure while preserving exact outputs;
- implement a hierarchy-scale policy against the existing group backend;
- migrate one backend to the named raw schema.

Avoid tasks that combine:

- architecture design;
- multivariate rewrite;
- hierarchy implementation;
- BayesRC implementation;
- public API redesign;
- performance optimization.

---

## 17. Required validation layers

Every implementation phase should use:

1. unit tests for one policy;
2. reduction tests to a simpler model;
3. deterministic seed/thread tests;
4. schema and field-identity tests;
5. simulation recovery tests;
6. regression tests for unaffected models;
7. package build/check;
8. benchmark tests only after correctness.

Performance optimization must not precede a validated reference implementation.

---

## 18. Decision log

The following decisions are adopted by this plan:

- models are compositions of policies;
- the statistical C++ core is language neutral;
- Rcpp is a thin current binding, not part of the model core;
- a future Python interface may use `pybind11` over the same engine;
- no Python package, stable C++ ABI, or Python dependency is required during the current refactor;
- data representation is separated from prior family;
- categorical hierarchy precedes overlapping log-additive hierarchy;
- hierarchy normalization uses a scale-preserving weighted geometric constraint;
- multivariate pattern enumeration is restricted to small trait sets;
- factor models handle larger trait collections;
- the evidence-factor model has a separate public interface and result class;
- the evidence-factor model begins at gene level;
- SNP-to-gene and pathway propagation begin as posterior post-processing;
- legacy multivariate code is audited before reuse;
- exact reproducibility is a framework invariant.

Open decisions to resolve during Phase 0:

- exact physical source layout compatible with the R package build;
- boundary between typed core results and R raw-schema conversion;
- canonical small-\(T\) covariance prior and update;
- exact supported residual covariance contracts;
- name and version of the future unified raw schema;
- public model-specification syntax;
- maximum default trait-pattern count;
- whether the first canonical MT model is dense individual-level or dense sufficient-statistics;
- migration and deprecation policy for `sblr()` legacy algorithms.

---

## 19. Immediate next actions

1. Add this document and its two companion matrices to `docs/dev`.
2. Run an audit-only Codex task over:
   - current multivariate native backends;
   - `sblr()` and related R wrappers;
   - chain/RNG infrastructure;
   - covariance updates;
   - raw return structures;
   - likelihood operators;
   - all Rcpp types, R RNG calls, Armadillo RNG calls, R console output, and R-specific errors inside statistical kernels;
   - candidate boundaries for typed C++ specifications and results.
3. Reconcile the audit with this plan.
4. Select the canonical multivariate BayesC likelihood and covariance model.
5. Create a staged implementation plan for Phases 1–3.
6. Do not implement hierarchical or factor models before the canonical multivariate and framework contracts are approved.
