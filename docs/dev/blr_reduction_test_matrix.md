# Unified BLR Framework: Reduction, Regression, and Validation Test Matrix

Phase 16A disposition note: the experimental packed-BED scheduled single-chain
route is not an exact reduction of canonical one-chain BayesC because its
historical logical-chain seed mapping differs; Phase 11B permanently records the
first numerical nonidentity. The experimental sparse route uses a distinct
stochastic null-marker scheduling policy and is not required to reduce exactly
to canonical scheduled BayesC. Both remain noncanonical research routes.

**Status:** Phase 12A hardened validation specification
**Date:** 2026-07-18
**Target location:** `docs/dev/blr_reduction_test_matrix.md`

## 1. Purpose

This document defines the tests required for an incremental architectural refactor of `sblr`.

The test strategy has two complementary goals:

1. preserve the behavior, speed, and memory efficiency of validated current kernels during migration;
2. validate newly implemented models through mathematical reductions and simulation recovery.

A refactor is not complete merely because it compiles or because posterior means look plausible.

---

## 2. Test classes

### 2.1 Exact regression test

Used when mathematics, operation order, RNG, and traversal are intentionally unchanged.

Expected:

```r
expect_identical(...)
```

or machine-level equality.

### 2.2 Numerical regression test

Used when representation or floating-point order changes without changing the target.

Expected:

- tight documented tolerance;
- dimension/name equality;
- stable derived identities.

### 2.3 Statistical equivalence test

Used for newly implemented samplers or intentionally changed RNG order.

Expected:

- repeated-simulation agreement;
- posterior moment agreement;
- interval coverage;
- calibrated state probabilities;
- no systematic bias.

### 2.4 Reduction test

A complex model is restricted to a simpler validated model.

### 2.5 Recovery test

Known simulated architecture is recovered.

### 2.6 Performance regression test

Wall time and throughput are compared against the pre-migration backend.

### 2.7 Memory regression test

Peak memory and allocation behavior are compared against the pre-migration backend.

### 2.8 Contract test

Checks dimensions, names, supported combinations, ownership, and informative failures.

---

## 3. General migration invariants

For a preserved existing kernel:

| ID | Invariant |
|---|---|
| MIG-01 | Public statistical target is unchanged |
| MIG-02 | Marker traversal is unchanged unless explicitly approved |
| MIG-03 | Burn-in and thinning semantics are unchanged |
| MIG-04 | Chain seed formulas are unchanged unless corrected in a separate task |
| MIG-05 | Residual update mathematics are unchanged |
| MIG-06 | Variance update mathematics are unchanged |
| MIG-07 | Output values and dimensions are unchanged during boundary migration |
| MIG-08 | No operator data are newly copied per chain |
| MIG-09 | No full marker-sample storage is introduced |
| MIG-10 | Runtime and memory regressions are measured |

---

## 4. Deterministic RNG tests

| ID | Test | Requirement |
|---|---|---|
| RNG-01 | Repeated fixed-seed call | identical stochastic output |
| RNG-02 | One-core versus multi-core | identical output for supported deterministic task scheduling |
| RNG-03 | Reversed core-order calls | no call-order effect |
| RNG-04 | Intervening unrelated fit | no contamination |
| RNG-05 | Different seed | changes at least one stochastic output |
| RNG-06 | Explicit chain seeds | exact documented mapping |
| RNG-07 | Thread assignment | chain trajectory independent of worker |
| RNG-08 | Distribution ownership | no global/static/persistent thread-local stochastic state |
| RNG-09 | Core RNG audit | no R, NumPy, or Armadillo global RNG |
| RNG-10 | Binding conversion | does not consume sampler RNG state |

For an RNG correction, compare posterior target statistically if exact trajectory changes.

---

## 5. Binding and language-neutral boundary tests

| ID | Test | Requirement |
|---|---|---|
| BIND-01 | Core headers | no Rcpp or Python headers |
| BIND-02 | Standard exceptions | core uses C++ exceptions |
| BIND-03 | Console separation | no R/Python output inside reusable core |
| BIND-04 | Typed spec conversion | R conversion preserves all values |
| BIND-05 | Typed result conversion | values and dimensions preserved |
| BIND-06 | Optional fields | absent core values map to actual R `NULL` |
| BIND-07 | Canonical orientation | core dimensions are independent of R dropping |
| BIND-08 | Ownership | borrowed input remains valid for native call |
| BIND-09 | Disk formats | new core resources do not require R serialization |
| BIND-10 | Future cross-binding | R and Python must call the same engine |

Future Python tests should compare typed specifications and typed results before presentation metadata.

---

## 6. CSR BayesC migration tests

The first architecture migration should use unscheduled CSR BayesC.

### 6.1 Exact old-versus-migrated comparison

Compare:

- marker means;
- PIP;
- final effects and states;
- residual scores;
- variance traces;
- covariance summaries;
- pi traces/final values;
- diagnostics;
- chains;
- metadata;
- formatted fit.

### 6.2 Configurations

Test:

- one trait;
- multiple independent traits;
- one chain;
- multiple chains;
- explicit chain seeds;
- `keep_chains = FALSE`;
- `keep_chains = TRUE`;
- `selection_s` disabled;
- fixed `selection_s`;
- LD-swap disabled;
- one-core and two-core execution.

### 6.3 Shape tests

Include:

- one marker;
- one trait;
- one retained sample;
- matrix preservation;
- actual `NULL` behavior.

---

## 7. Scalar model reductions

| ID | Complex model | Restriction | Reference |
|---|---|---|---|
| RED-ST-01 | BayesR | one non-null component | BayesC |
| RED-ST-02 | BayesRC | fixed intercept-only probabilities | BayesR |
| RED-ST-03 | Fixed marker prior | constant pi and scale | BayesC |
| RED-ST-04 | Group probability | one group | global BayesC |
| RED-ST-05 | Group scale | one group with multiplier 1 | BayesC |
| RED-ST-06 | Learned annotation probability | coefficients zero | global BayesC |
| RED-ST-07 | Learned annotation scale | coefficients zero | unit/fixed scale reference |
| RED-ST-08 | Hierarchical BayesC | all multipliers 1 | BayesC |
| RED-ST-09 | Hierarchical BayesC | one layer | approved group/hierarchy reference |
| RED-ST-10 | Hierarchical BayesR | hierarchy 1 | BayesR |
| RED-ST-11 | Hierarchical BayesRC | hierarchy 1 | BayesRC |
| RED-ST-12 | Hierarchical BayesRC | one non-null component | hierarchical BayesC |

---

## 8. Scale-policy identities

### 8.1 Fixed scale

For fixed \(q_j\), migrated policy output must reproduce the current fixed-prior backend.

### 8.2 `selection_s`

Preserve the current scale definition exactly during migration.

Test:

- fixed \(S\);
- sampled \(S\);
- bounds;
- acceptance diagnostics;
- MAF alignment;
- old-versus-migrated equality.

Do not assume \(S=0\) means unit scale unless that is the documented current formula.

### 8.3 Hierarchy normalization

For each hierarchy layer:

\[
\exp\left[
\frac{\sum_g n_{hg}\log\theta_{hg}}
{\sum_g n_{hg}}
ight]=1.
\]

Scale preservation:

\[
v_b^{\mathrm{before}}q_j^{\mathrm{before}}
=
v_b^{\mathrm{after}}q_j^{\mathrm{after}}.
\]

Test all markers within numerical tolerance.

---

## 9. Probability and state identities

### Binary state

\[
0\le dm_j\le 1.
\]

When retained state samples are available:

\[
dm_j = \frac{1}{S}\sum_s d_j^{(s)}.
\]

### Mixture state

\[
\sum_k P(c_j=k\mid y)=1.
\]

If component 0 is null:

\[
dm_j=1-P(c_j=0\mid y).
\]

Component counts:

\[
\sum_k n_k=m.
\]

### Trait patterns

\[
\sum_r P(r_j=r\mid y)=1.
\]

Trait-specific PIP:

\[
dm_{jt}
=
\sum_r I(r_t=1)P(r_j=r\mid y).
\]

---

## 10. Operator equivalence tests

Use tiny exact fixtures where the same model can be represented through multiple operators.

| ID | Operator A | Operator B | Model |
|---|---|---|---|
| OP-01 | dense individual | dense sufficient statistics | ST BayesC |
| OP-02 | dense individual | BED | ST BayesC |
| OP-03 | dense summary | CSR | ST BayesC |
| OP-04 | CSR | exact/no-truncation block eigen fixture | ST BayesC |
| OP-05 | CSR | BED-derived sufficient statistics | BayesR |
| OP-06 | CSR | BED | BayesRC under matched prior |
| OP-07 | dense MT individual | dense MT sufficient statistics | MT BayesC |
| OP-08 | dense MT | later CSR MT | MT BayesC |

Every operator test must verify:

- marker order;
- trait order;
- scaling;
- sample-size convention;
- allele-frequency convention;
- residual rebuild identity.

---

## 11. Multivariate tests

The new multivariate engine requires a complete suite.

### 11.1 Reductions

| ID | Complex model | Restriction | Reference |
|---|---|---|---|
| MT-RED-01 | MT BayesC | \(T=1\) | scalar BayesC |
| MT-RED-02 | MT diagonal covariance | independent states/covariance | separate scalar fits |
| MT-RED-03 | MT shared pattern | null/all-active only | shared-indicator reference |
| MT-RED-04 | Explicit pattern set | one active pattern plus null | restricted binary reference |
| MT-RED-05 | MT BayesR | \(T=1\) | scalar BayesR |
| MT-RED-06 | MT BayesRC | \(T=1\) | scalar BayesRC |
| MT-RED-07 | Fixed covariance | updates disabled | supplied covariance preserved |
| MT-RED-08 | Full covariance | off-diagonal zero | diagonal reference |

### 11.2 Covariance tests

- positive definiteness;
- correct full residual score use;
- exact shared-design residual cross-product;
- fixed/update flag behavior;
- covariance naming;
- simulation recovery;
- prior sensitivity;
- unsupported overlap/design errors.

### 11.3 State-space scaling

- no automatic impractical \(2^T\) enumeration;
- duplicate pattern rejection;
- null pattern requirement;
- restricted pattern validation;
- explicit maximum/default behavior.

---

## 12. Factor-model tests

| ID | Complex model | Restriction | Reference |
|---|---|---|---|
| FAC-01 | Sparse factor model | inclusion fixed to 1 | dense factor model |
| FAC-02 | Annotation membership | coefficients zero | global membership |
| FAC-03 | Multi-view model | one view | single-view model |
| FAC-04 | SNP-to-gene mapping | identity matrix | original loadings |
| FAC-05 | Gene-to-pathway mapping | identity matrix | gene loadings |

Additional tests:

- reconstruction invariance under sign normalization;
- reconstruction invariance under scale normalization;
- factor ordering;
- signed-permutation alignment;
- known-factor simulation recovery;
- observation-uncertainty calibration.

---

## 13. Raw result and formatting tests

During migration of existing ST models:

- `stblr_raw_v1` remains valid;
- canonical dimensions remain unchanged;
- names remain unchanged;
- optional fields retain current behavior;
- formatters preserve matrix dimensions;
- no positional output is introduced.

For a future joint MT or factor schema:

- define named namespaces;
- document dimensions before implementation;
- test version detection;
- test conversion from typed core result;
- do not overload historical positional MT slots.

---

## 14. Performance regression tests

Performance is a preservation gate.

### 14.1 Benchmarks

For each migrated backend, compare:

- pre-migration implementation;
- migrated implementation.

Record:

- wall time;
- initialization time;
- marker updates per second where available;
- one-core and multiple-core behavior;
- one-chain and multi-chain behavior.

### 14.2 Acceptance

A practical default:

- differences within approximately 5–10% may be benchmark noise;
- larger regression requires repeated measurement and explanation;
- public replacement should not occur with an unexplained material slowdown.

Do not optimize from one tiny benchmark.

### 14.3 Benchmark fixtures

Include:

- tiny correctness fixture;
- moderate representative fixture;
- large realistic fixture where available;
- minimal-output run;
- rich-output run.

Benchmarks are not required in fast CI.

---

## 15. Memory regression tests

Measure peak resident memory for:

- legacy backend;
- migrated backend;
- one chain;
- multiple chains;
- minimal output;
- rich output.

Required invariants:

| ID | Invariant |
|---|---|
| MEM-01 | No full immutable operator copy per chain |
| MEM-02 | No full marker-by-sample storage by default |
| MEM-03 | Marker-summary memory independent of retained sample count |
| MEM-04 | Optional component/factor output measured separately |
| MEM-05 | No meaningful unexplained peak-memory regression |
| MEM-06 | One-marker/one-trait shapes do not trigger hidden copies |

Allocation instrumentation may be added for newly written hot loops, but existing efficient kernels need not be rewritten solely to satisfy an allocation-counter design.

---

## 16. Chain aggregation tests

When chains are not retained:

- raw and formatted absence follows the documented contract.

When chains are retained:

- dimensions are documented;
- top-level posterior means equal documented chain aggregation;
- one-chain aggregation is a degenerate special case;
- matrices remain matrices;
- chain order follows explicit seed/task order.

---

## 17. Diagnostic identities

For applicable ST outputs:

\[
v_{ld}=v_g-v_{le}.
\]

For CPO:

\[
\mathrm{mean\_log\_cpo}
=
\frac{\mathrm{log\_cpo}}{n_{\mathrm{used}}}.
\]

For thinning:

\[
n_{\mathrm{samples}}
=
\#\{i\ge nburn:(i-nburn)\bmod nthin=0\}.
\]

Prior probabilities and realized state frequencies must remain distinct.

---

## 18. Input and contract tests

Test:

- missing inputs;
- invalid dimensions;
- duplicate IDs;
- missing IDs;
- non-finite values;
- invalid patterns;
- invalid covariance;
- invalid hierarchy groups;
- invalid component scales;
- invalid seeds;
- invalid MCMC control;
- unsupported operator/model combination;
- unsupported residual covariance/design combination;
- malformed typed specifications.

Errors should identify the failing field and expected contract.

---

## 19. CI tiers

### Tier 1 — fast

- typed specification validation;
- exact regression fixtures;
- RNG one-core/two-core;
- raw schema;
- small policy reductions.

### Tier 2 — full package

- complete `testthat`;
- package build/check;
- migrated operator/model regression;
- documentation consistency where changed.

### Tier 3 — extended

- simulation recovery;
- realistic operator equivalence;
- performance benchmarks;
- peak-memory benchmarks;
- long-chain factor alignment.

---

## 20. Migration readiness gates

### Boundary ready

- typed specs/results exist;
- binding conversion tested;
- no Rcpp in reusable core;
- baseline fixtures recorded.

### Existing kernel migrated

- exact/statistical behavior preserved;
- reproducibility preserved;
- runtime has no unexplained material regression;
- peak memory has no meaningful regression;
- responsibilities are clearer.

### Shared utility ready

- used by at least two migrated models;
- reductions protect both;
- no model-specific semantics hidden in generic names.

### Legacy removal ready

- public route has stabilized;
- full tests pass;
- benchmark and memory gates pass;
- no unsupported feature remains in legacy only;
- removal is explicitly approved.

### New multivariate ready

- all MT reductions pass;
- covariance contract approved;
- simulated recovery passes;
- exact shared-design reference validated.

### Factor ready

- alignment and identifiability tests pass;
- simulation recovery passes;
- output uncertainty is documented.
