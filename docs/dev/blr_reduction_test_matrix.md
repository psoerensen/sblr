# Unified BLR Framework: Reduction and Validation Test Matrix

**Status:** Draft validation specification  
**Date:** 2026-07-13  
**Target location:** `docs/dev/blr_reduction_test_matrix.md`

## 1. Purpose

This document defines the validation strategy for the unified BLR and factor framework.

The framework should be developed through a lattice of reductions:

> Every complex model must reduce to one or more simpler validated models under explicit parameter restrictions.

Reduction tests are stronger than isolated smoke tests because they verify:

- likelihood consistency;
- state-space consistency;
- prior-scale consistency;
- RNG ordering;
- chain aggregation;
- output semantics;
- operator equivalence.

---

## 2. Test categories

### Exact reduction

Expected to agree exactly or at machine-level tolerance under matched RNG trajectories.

Typical tolerance:

```r
tolerance <- 1e-12
```

Use exact equality when object types and floating-point operation order are identical.

### Statistical reduction

The same posterior target is represented through different computational paths, but floating-point ordering differs.

Use:

- tight numeric tolerance;
- matched posterior summaries;
- distributional validation over repeated simulations.

### Identity test

Tests an algebraic identity within one output, for example:

\[
v_{ld}=v_g-v_{le}.
\]

### Recovery test

Tests whether simulated architecture is recovered directionally or quantitatively.

### Regression test

Confirms an unaffected existing model remains unchanged.

### Contract test

Checks dimensions, names, `NULL` behavior, supported combinations, and informative errors.

---

## 3. Deterministic framework invariants

These apply to every native MCMC backend and must be independent of the language binding.

| ID | Invariant | Required comparison |
|---|---|---|
| RNG-01 | Repeated identical-seed calls are identical | all stochastic output |
| RNG-02 | One-core and multi-core execution are identical | all stochastic output |
| RNG-03 | Reversing thread-count call order has no effect | all stochastic output |
| RNG-04 | An unrelated intervening fit does not contaminate the next fit | reference fit equality |
| RNG-05 | A different seed changes at least one stochastic output | one or more posterior fields |
| RNG-06 | Chain seed rules are independent of thread assignment | retained chain equality |
| RNG-07 | No stateful distribution survives chain lifetime | code audit plus tests |
| RNG-08 | No Armadillo/R global RNG is used inside chain execution | code audit |
| RNG-09 | Binding conversion does not consume or alter sampler RNG state | adapter/core equality |
| RNG-10 | A future binding must preserve the same resolved seed and chain rules | cross-binding contract |

Required fields depend on model, but may include:

- posterior marker means;
- state probabilities;
- variance traces;
- component probabilities;
- hierarchy multipliers;
- annotation coefficients;
- covariance traces;
- factor loadings;
- CPO diagnostics.

---

## 4. Core regression reductions

| ID | Complex model | Restriction | Reference model | Required outputs | Type |
|---|---|---|---|---|---|
| RED-01 | MT BayesC | \(T=1\) | ST BayesC | marker effects, PIP, variances, state | Exact |
| RED-02 | MT diagonal BayesC | independent states/covariance | separate ST fits | trait-specific posterior fields | Exact/statistical |
| RED-03 | MT shared BayesC | only null/all-active patterns | shared-indicator reference | shared PIP, effect covariance | Exact |
| RED-04 | Trait-pattern model | patterns = null + one selected pattern | binary restricted model | state probabilities, effects | Exact |
| RED-05 | BayesR | one non-null component | BayesC with matched scale | effects, PIP, variances | Exact |
| RED-06 | BayesRC | fixed intercept-only annotation prior | fixed-pi BayesR | effects, `dm`, `comp_prob`, traces | Exact |
| RED-07 | MT BayesR | \(T=1\) | ST BayesR | component probabilities and effects | Exact |
| RED-08 | MT BayesRC | \(T=1\) | ST BayesRC | component and annotation outputs | Exact |

---

## 5. Hierarchy-scale reductions

| ID | Complex model | Restriction | Reference model | Required outputs | Type |
|---|---|---|---|---|---|
| HIER-01 | Hierarchical BayesC | all multipliers fixed at 1 | ordinary BayesC | all common posterior fields | Exact |
| HIER-02 | Hierarchical BayesC | one hierarchy layer | group BayesC | effects, PIP, group multiplier | Exact/statistical |
| HIER-03 | Hierarchical BayesC | fixed effective \(q_j\) | fixed marker-scale BayesC | effects, PIP, variances | Exact |
| HIER-04 | Multi-layer hierarchy | freeze all but one layer | one-layer hierarchy | selected layer outputs | Exact/statistical |
| HIER-05 | Composite hierarchy + MAF | hierarchy = 1 | `selection_s` model | common posterior fields | Exact |
| HIER-06 | Composite hierarchy + MAF | \(S=0\) or matched neutral setting | hierarchy-only model | common posterior fields | Exact |
| HIER-07 | Hierarchical MT BayesC | \(T=1\) | hierarchical ST BayesC | all common and hierarchy fields | Exact |
| HIER-08 | Hierarchical BayesR | hierarchy = 1 | BayesR | component outputs | Exact |
| HIER-09 | Hierarchical BayesRC | hierarchy = 1 | BayesRC | component and annotation outputs | Exact |
| HIER-10 | Hierarchical BayesRC | one non-null component | hierarchical BayesC | effects and hierarchy | Exact/statistical |

### Scale-preservation identity

For every hierarchy normalization step:

\[
v_b^{\text{before}}q_j^{\text{before}}
=
v_b^{\text{after}}q_j^{\text{after}}.
\]

Test for all markers within machine tolerance.

### Layer normalization identity

For each layer \(h\):

\[
\exp\left[
\frac{\sum_g n_{hg}\log\theta_{hg}}
     {\sum_g n_{hg}}
\right]
=1.
\]

---

## 6. Annotation-probability reductions

| ID | Complex model | Restriction | Reference | Required outputs | Type |
|---|---|---|---|---|---|
| ANN-01 | Learned BayesC probability | annotation coefficients = 0 | global BayesC | PIP and effects | Exact/statistical |
| ANN-02 | Group probability | one group | global BayesC | global/group probability | Exact |
| ANN-03 | BayesRC stick-breaking | intercept-only fixed coefficients | global BayesR | full component outputs | Exact |
| ANN-04 | Annotation preprocessing | shuffled rows with IDs | manually aligned annotations | all posterior fields | Exact |
| ANN-05 | Annotation preprocessing | extra annotation rows | selected aligned subset | all posterior fields | Exact |
| ANN-06 | Binary annotation scaling disabled | raw 0/1 input | manually prepared input | annotation matrix and posterior | Exact |
| ANN-07 | Annotation standardization | wrapper preprocessing | manually standardized input | native posterior | Exact |

---

## 7. Data-operator equivalence

Use tiny deterministic fixtures where the same \(X\), \(y\), marker order, scaling, and sufficient statistics are available.

| ID | Operator A | Operator B | Model | Required comparison | Type |
|---|---|---|---|---|---|
| OP-01 | dense individual | dense summary | ST BayesC | all common fields | Exact/statistical |
| OP-02 | dense individual | packed BED | ST BayesC | all common fields | Exact/statistical |
| OP-03 | dense summary | CSR | ST BayesC | effects, PIP, variance traces | Exact/statistical |
| OP-04 | CSR | block-eigen with no truncation/identity fixture | ST BayesC | all common fields | Exact/statistical |
| OP-05 | CSR | packed BED-derived sufficient statistics | ST BayesR | component outputs | Statistical |
| OP-06 | CSR | BED | BayesRC under matched prior | component and annotation outputs | Statistical |
| OP-07 | dense MT | dense MT sufficient statistics | canonical MT BayesC | pattern/effect/covariance outputs | Exact/statistical |
| OP-08 | dense MT summary | CSR MT | canonical MT BayesC | common fields | Statistical |

Each operator test must also verify:

- identical marker order;
- identical trait order;
- identical scaling convention;
- identical allele-frequency interpretation;
- exact residual rebuild identity where applicable.

---

## 8. Multivariate covariance tests

| ID | Test | Requirement |
|---|---|---|
| COV-01 | Positive definiteness | every sampled covariance is finite, symmetric, and SPD |
| COV-02 | Fixed covariance | disabling covariance updates preserves supplied covariance |
| COV-03 | Diagonal reduction | off-diagonal fixed at zero gives independent-trait reference |
| COV-04 | \(T=1\) reduction | covariance update reduces to scalar variance update |
| COV-05 | Same-design residual covariance | simulated residual covariance recovered |
| COV-06 | Unsupported design | full residual covariance errors without required cross-trait sufficient statistics |
| COV-07 | Covariance naming | rows/columns match trait names |
| COV-08 | Correlation aliases | `cov2cor` summaries agree with covariance matrices |
| COV-09 | Chain aggregation | aggregate covariance equals documented chain combination |
| COV-10 | Prior sensitivity | stronger prior shrinks estimates in expected direction |

A covariance implementation must not be validated only by checking that it is SPD.

---

## 9. State and component identities

### Binary state

\[
dm_j = E[d_j\mid y].
\]

Check:

- `dm` lies in \([0,1]\);
- final state is 0/1;
- posterior counts agree with retained samples when available.

### Mixture state

For each marker:

\[
\sum_k P(c_j=k\mid y)=1.
\]

If component zero is null:

\[
dm_j = 1 - P(c_j=0\mid y).
\]

Posterior mean component index:

\[
\bar c_j = \sum_k kP(c_j=k\mid y).
\]

Component counts:

\[
n_k = \sum_j P(c_j=k\mid y),
\qquad
\sum_k n_k = m.
\]

### Trait-pattern state

For each marker:

\[
\sum_r P(r_j=r\mid y)=1.
\]

Trait-specific PIP:

\[
dm_{jt}
=
\sum_r \mathbb{1}(r_t=1)P(r_j=r\mid y).
\]

Shared-pattern counts must agree with the posterior pattern-probability matrix.

---

## 10. Variance and diagnostic identities

For applicable ST-BLR outputs:

\[
v_{ld}=v_g-v_{le}.
\]

For CPO:

\[
\text{mean\_log\_cpo}
=
\frac{\text{log\_cpo}}{n_{\text{used}}}.
\]

For thinning:

\[
n_{\text{samples}}
=
\#\{i\ge nburn:(i-nburn)\bmod nthin=0\}.
\]

For prior component probabilities:

- `pi$final` must equal the documented marker-average final prior;
- `pim` must equal the posterior mean prior;
- `pis` must follow its documented per-iteration definition;
- prior probabilities must not be confused with realized state frequencies.

---

## 11. Chain aggregation tests

When `keep_chains = FALSE`:

```r
expect_null(raw$chains)
```

When `keep_chains = TRUE`:

- one list per trait or model-defined top-level unit;
- one list per chain;
- dimensions documented and tested;
- top-level posterior means equal the documented average of chain summaries;
- matrices remain matrices where semantically required;
- vectors may be compared after explicit `drop()` only in tests where shape is not semantic.

Required aggregation checks may include:

- `bm`;
- `dm`;
- component probabilities;
- pattern probabilities;
- annotation coefficients;
- hierarchy multipliers;
- covariance summaries;
- factor loadings after alignment.

---

## 12. Factor-model reductions

| ID | Complex model | Restriction | Reference | Required outputs | Type |
|---|---|---|---|---|---|
| FAC-01 | Gaussian factor model | \(K=1\) | rank-one Bayesian regression | reconstruction and loadings | Statistical |
| FAC-02 | Sparse factor model | inclusion fixed to 1 | dense Gaussian factor model | loadings and reconstruction | Exact/statistical |
| FAC-03 | Annotation-informed membership | coefficients fixed to 0 | global-membership factor model | membership probabilities | Exact/statistical |
| FAC-04 | Multi-view factor model | one view | single-view factor model | common factor outputs | Exact |
| FAC-05 | Factor-analytic MT-BLR | no idiosyncratic factor contribution | diagonal/full covariance reference | covariance/effect outputs | Statistical |
| FAC-06 | Gene factor model | identity row mapping | direct factor outputs | loadings and reconstruction | Exact |
| FAC-07 | SNP-to-gene propagation | \(S=I\) | SNP-factor loadings | propagated loadings | Exact |
| FAC-08 | Gene-to-pathway propagation | \(C=I\) | gene-factor loadings | pathway loadings | Exact |
| FAC-09 | Fixed uncertain maps | one map with probability 1 | fixed mapping | posterior propagation | Exact |

---

## 13. Factor identifiability tests

### Scale normalization

For each factor \(k\), apply the selected convention, for example:

\[
\|V_{\cdot k}\|_2=1.
\]

Verify that normalization leaves:

\[
UV^\top
\]

unchanged.

### Sign convention

The selected anchor loading, such as the largest absolute trait loading, must be positive.

Verify reconstruction remains unchanged after sign normalization.

### Factor ordering

Factors should be sorted by a documented quantity, such as:

- explained variance;
- factor scale;
- posterior expected contribution.

### Cross-chain alignment

Construct chains with known:

- factor permutations;
- sign flips;
- small noise.

Verify the alignment procedure recovers the known matching.

Do not average unaligned factor matrices.

---

## 14. Simulation recovery suite

## 14.1 Hierarchical BayesC

Simulate:

- multiple categorical layers;
- known enriched and depleted groups;
- enough markers per group;
- realistic sparsity;
- moderate LD.

Metrics:

- rank correlation between true and estimated multipliers;
- coverage of posterior intervals;
- directional enrichment recovery;
- calibrated marker PIP;
- prediction/effect recovery;
- invariance to layer label ordering.

## 14.2 Multivariate BayesC

Simulate:

- independent trait effects;
- fully shared effects;
- selected partial-sharing patterns;
- positive and negative effect correlations;
- diagonal and full residual covariance where supported.

Metrics:

- trait-pattern classification;
- trait-specific PIP calibration;
- effect RMSE;
- covariance recovery;
- credible interval coverage.

## 14.3 Hierarchical BayesRC

Simulate:

- annotation-driven component membership;
- hierarchy-driven scale;
- both signals together;
- null annotation signal;
- null hierarchy signal.

Metrics:

- component probability calibration;
- hierarchy multiplier recovery;
- separation of membership and scale effects;
- BayesRC/BayesR reduction under null restrictions.

## 14.4 Gene evidence-factor model

Simulate:

- known factor number;
- sparse gene-factor loadings;
- trait-factor loadings;
- annotation enrichment;
- heteroscedastic observation uncertainty;
- optional molecular anchor views.

Metrics:

- factor subspace recovery;
- aligned loading correlation;
- membership PIP calibration;
- reconstruction error;
- annotation-effect recovery;
- factor-number sensitivity.

---

## 15. Input and contract tests

Every public/internal resolver should test:

- missing required inputs;
- duplicate marker or row IDs;
- missing selected IDs;
- invalid trait names;
- non-finite values;
- invalid dimensions;
- invalid state patterns;
- duplicate state patterns;
- no null state where required;
- invalid component scales;
- invalid covariance matrices;
- invalid hierarchy groups;
- empty hierarchy groups;
- invalid prior parameters;
- unsupported data/covariance combinations;
- unsupported operator/model combinations;
- invalid MCMC controls;
- invalid chain seeds;
- unsupported public arguments.

Errors should identify:

- the failing object;
- the expected dimension or contract;
- a small example of offending IDs where relevant.

---

## 16. Raw-schema tests

Every raw result must include:

```r
raw$schema
raw$meta
raw$marker
raw$trace
raw$diagnostics
```

Additional required namespaces depend on the model:

| Model | Required namespace |
|---|---|
| BayesC | `state` or binary marker state representation |
| BayesR/RC | `component` |
| Trait-pattern model | `pattern` |
| Annotation model | `annotation` |
| Hierarchy model | `hierarchy` |
| Multivariate model | `trait`, optionally `residual` |
| Factor model | `factor` |
| Multiple chains | `chains` or explicit `NULL` |

Schema tests must verify:

- canonical dimension orientation;
- names;
- no ambiguous positional interpretation;
- optional fields are actual `NULL`;
- formatted aliases map to authoritative raw fields.

---

## 17. Regression-test requirements

A framework change must run tests for:

- the directly changed model;
- the simpler reference model used in reduction;
- all operators using the changed shared code;
- chain/RNG tests;
- raw schema;
- formatter;
- public interface;
- selection and annotation models when shared utilities change;
- block-eigen and LD-swap when operator code changes;
- full package suite.

No tolerance may be weakened solely to make a refactor pass.

Any changed tolerance must have a documented numerical reason and a separate validation of the posterior target.

---

## 18. Performance validation

Performance tests begin only after correctness.

Track:

- wall time;
- peak memory;
- residual rebuild count;
- marker updates per second;
- covariance-update time;
- state-evaluation count;
- factor alignment time;
- scaling with markers, traits, states, and chains.

A performance optimization is acceptable only when:

- all statistical reductions still pass;
- RNG behavior remains documented;
- output schema remains unchanged;
- benchmark improvement is material.

---

## 19. CI tiers

### Tier 1 — fast pull-request tests

- policy unit tests;
- tiny exact reductions;
- RNG one-core/two-core tests;
- schema and input contracts;
- short deterministic simulations.

### Tier 2 — package regression

- full `testthat` suite;
- package build;
- package check without external-data workflows;
- rendered developer documentation where applicable.

### Tier 3 — extended validation

- larger simulation recovery;
- operator equivalence on realistic LD;
- performance benchmarks;
- longer multi-chain factor alignment;
- external example workflows.

---

## 20. Phase readiness gates

### Architecture-ready

- capability matrix approved;
- reduction matrix approved;
- multivariate audit completed;
- canonical covariance model selected.

### Framework-core-ready

- model specifications resolve correctly;
- chain/RNG infrastructure passes deterministic tests;
- named raw output is stable.

### Canonical MT-ready

- \(T=1\) reduction passes;
- diagonal/reference reductions pass;
- covariance recovery passes;
- unsupported residual-covariance designs error.

### Hierarchy-ready

- HIER-01 through HIER-07 pass;
- normalization identities pass;
- simulation recovers multiple layers.

### Hierarchical BayesRC-ready

- HIER-08 through HIER-10 pass;
- ANN-03 passes;
- membership and scale signals separate in simulation.

### Factor-core-ready

- FAC-01 through FAC-05 pass;
- sign/permutation alignment passes;
- reconstruction is invariant to normalization.

### Evidence-factor-ready

- calibrated simulated recovery;
- observation uncertainty handled;
- benchmark against simpler alternatives;
- mechanism-informed interpretation documented without causal overclaiming.

---

## 21. Required implementation report

Every Codex implementation report should include:

- branch and initial Git status;
- files inspected;
- files changed;
- model/policy implemented;
- mathematical formulas changed;
- RNG effects;
- output-schema effects;
- exact reduction tests run;
- simulation tests run;
- pass/fail/skip counts;
- package build/check status;
- unresolved blockers;
- readiness status;
- confirmation that no unrelated model was changed.

---

## 22. Language-binding and core-boundary tests

These tests protect the option of a future Python interface without requiring Python to be implemented now.

### Current required tests

| ID | Test | Requirement |
|---|---|---|
| BIND-01 | Core header audit | Core statistical headers do not include Rcpp headers |
| BIND-02 | Core result audit | Migrated samplers return typed C++ results rather than `Rcpp::List` |
| BIND-03 | Exception audit | Core throws standard C++ exceptions |
| BIND-04 | Console audit | Core does not call `Rcpp::Rcout`, `Rprintf`, or Python logging |
| BIND-05 | RNG audit | Core does not call R, Armadillo global, or NumPy runtime RNG |
| BIND-06 | Dimension contract | Core dimensions match canonical orientation before R naming/formatting |
| BIND-07 | Adapter round-trip | R specification conversion preserves all typed specification fields |
| BIND-08 | Result conversion | R result conversion preserves all values and dimensions |
| BIND-09 | Optional-field conversion | Missing typed fields become actual R `NULL` |
| BIND-10 | Error translation | Standard C++ validation errors become informative R errors |
| BIND-11 | Data ownership | Borrowed views remain valid for the complete native call |
| BIND-12 | Disk format | New persisted core data does not require R serialization |

These tests should be introduced as each backend migrates. Legacy backends are not required to satisfy the complete boundary immediately, but the audit must identify every violation.

### Future cross-binding tests

When a Python binding is eventually added, require:

| ID | Test | Requirement |
|---|---|---|
| PY-01 | Resolved specification equality | R and Python adapters create equivalent typed C++ specs |
| PY-02 | Fixed-seed equality | Same core call produces identical numerical results |
| PY-03 | Dimension equality | R arrays and NumPy arrays represent the same canonical dimensions |
| PY-04 | Optional values | R `NULL` and Python `None` correspond to the same absent core field |
| PY-05 | Exception semantics | Both bindings report the same core validation failure |
| PY-06 | Disk resource equality | Both bindings read the same CSR/BED resource metadata |
| PY-07 | Serialization independence | Core computation does not depend on RDS or pickle |
| PY-08 | Performance sanity | Binding conversion is not the dominant runtime for realistic analyses |

Cross-binding equality should be tested at the typed C++ result level where practical, before language-specific names, classes, or presentation metadata are added.

---

## 23. Binding readiness gates

### R-binding boundary ready

- typed C++ specification exists;
- typed C++ result exists for the migrated backend;
- Rcpp is confined to adapter files;
- deterministic and reduction tests still pass;
- R output remains backward compatible where required.

### Python option preserved

- no R-specific types occur in the migrated core;
- no R/Armadillo global RNG occurs in the migrated core;
- standard C++ exceptions are used;
- canonical dimension conventions are documented;
- disk-backed input formats are language neutral;
- no Python dependency or ABI commitment has been introduced.

### Future Python binding ready

- all `PY-*` tests pass;
- the Python binding calls the same core entry points;
- no statistical code is duplicated;
- R and Python documentation describe the same model semantics.
