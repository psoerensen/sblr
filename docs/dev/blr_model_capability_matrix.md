# Unified BLR Framework: Model Capability Matrix

**Status:** Draft working matrix  
**Date:** 2026-07-13  
**Target location:** `docs/dev/blr_model_capability_matrix.md`

## 1. Purpose

This matrix defines how current and planned models decompose into reusable framework components.

It is intended to answer:

- which model capabilities currently exist;
- which components are reusable;
- which combinations are planned;
- which combinations are invalid or intentionally deferred;
- which reductions should hold.

Status labels:

- **Current** — implemented and validated in the current package.
- **Partial** — some required components exist, but not through the proposed unified framework.
- **Planned** — accepted target for staged implementation.
- **Deferred** — technically possible but not an initial priority.
- **Unsupported** — mathematically invalid under the stated data contract or outside package scope.

---

## 2. Component vocabulary

### State spaces

| ID | Meaning |
|---|---|
| `binary` | Null versus active marker |
| `trait_pattern` | One binary activity mask across traits |
| `mixture` | One BayesR/RC effect-size component |
| `pattern_mixture` | Trait pattern × effect-size component |
| `factor` | Sparse marker/gene membership in one or more latent factors |

### Probability policies

| ID | Meaning |
|---|---|
| `global_binary` | Global BayesC inclusion probability |
| `global_dirichlet` | Global multinomial mixture probabilities |
| `group` | Group-specific probabilities |
| `annot_logit` | Annotation-dependent Bernoulli/logit probabilities |
| `annot_probit_stick` | Annotation-dependent mutually exclusive component probabilities |
| `factor_annot` | Independent annotation-informed factor membership probabilities |

### Scale policies

| ID | Meaning |
|---|---|
| `unit` | No marker-specific scaling |
| `fixed_marker` | User-supplied marker multiplier |
| `maf_s` | MAF-dependent `selection_s` multiplier |
| `group` | One categorical group layer |
| `hierarchy` | Multiple multiplicative categorical layers |
| `log_additive` | Overlapping annotation incidence on log scale |
| `component` | BayesR/RC fixed component multipliers \(\gamma_k\) |
| `composite` | Product of two or more scale policies |

### Trait covariance policies

| ID | Meaning |
|---|---|
| `scalar` | Single trait |
| `diagonal` | Independent trait-specific variances |
| `full` | Full trait covariance |
| `factor_analytic` | Low-rank plus diagonal covariance |

### Data operators

| ID | Meaning |
|---|---|
| `dense_individual` | Dense individual-level \(X,y\) |
| `bed` | Packed PLINK BED |
| `dense_summary` | Dense sufficient statistics |
| `csr` | Sparse-CSR LD sufficient statistics |
| `block_eigen` | Block-eigen LD approximation |
| `evidence` | Gene/SNP-by-trait evidence matrix |

---

## 3. Regression model matrix

| Model | State | Probability | Scale | Trait covariance | Data | Status | Primary output |
|---|---|---|---|---|---|---|---|
| ST BayesC | binary | global_binary | unit | scalar | CSR, BED | Current | `bm`, `dm`, variance traces |
| ST BayesC + fixed marker priors | binary | fixed/global | fixed_marker | scalar | CSR | Current | effects and resolved marker priors |
| ST learned annotation BayesC | binary | annot_logit/global | optional learned scale | scalar | CSR | Current | annotation effects and PIPs |
| ST group BayesC | binary | group | group | scalar | CSR | Current | group probabilities and multipliers |
| ST BayesC + `selection_s` | binary | global_binary | maf_s | scalar | CSR | Current | effects and sampled/fixed \(S\) |
| ST BayesR | mixture | global_dirichlet | component | scalar | CSR, BED | Current | component probabilities |
| ST SBayesRC | mixture | annot_probit_stick | component | scalar | CSR | Current | component probabilities and annotation effects |
| ST BED BayesRC | mixture | annot_probit_stick | component | scalar | BED | Current | individual-level BayesRC outputs |
| ST hierarchical BayesC | binary | global_binary | hierarchy | scalar | CSR first | Planned | layer multipliers and effective scale |
| ST hierarchical BayesR | mixture | global_dirichlet | component × hierarchy | scalar | CSR first | Planned | components plus hierarchy |
| ST hierarchical BayesRC | mixture | annot_probit_stick | component × hierarchy | scalar | CSR first | Planned | annotation membership plus hierarchy |
| ST overlapping hierarchy | binary/mixture | compatible | log_additive | scalar | CSR | Deferred | log-scale annotation coefficients |
| MT BayesC | trait_pattern | global_dirichlet | unit | full/diagonal | legacy dense/summary | Partial | pattern probabilities and covariance |
| Canonical MT BayesC | trait_pattern | global_dirichlet | unit | selected canonical policy | dense then CSR | Planned | named multivariate raw output |
| MT independent BayesC | trait_pattern or separate binary | global_binary | unit | diagonal | dense/CSR | Planned | trait-specific effects |
| MT shared BayesC | shared/null pattern | global_binary | unit | full | shared design | Planned | shared PIP and correlated effects |
| MT hierarchical BayesC | trait_pattern | global_dirichlet | hierarchy | full/diagonal | dense/CSR | Planned | patterns, covariance, hierarchy |
| MT BayesR | pattern_mixture | global_dirichlet | component | full/diagonal | dense/CSR | Planned | pattern-component probabilities |
| MT BayesRC | pattern_mixture | annotation-informed | component | full/diagonal | dense/CSR | Deferred | marker pattern-component probabilities |
| MT hierarchical BayesRC | pattern_mixture | annotation-informed | component × hierarchy | full/diagonal | dense/CSR | Deferred | full composed prior |
| Factor-analytic MT-BLR | factor | factor prior | factor/composite | factor_analytic | dense/CSR | Planned later | marker-factor and trait loadings |

---

## 4. Evidence-factor model matrix

| Model | Row unit | Observation model | Membership | Annotation role | Mapping | Status |
|---|---|---|---|---|---|---|
| Gene-level Gaussian factor model | gene | Gaussian with trait residual variance | dense Gaussian | none | none | First planned factor model |
| Sparse gene factor model | gene | Gaussian | spike-and-slab | none | none | Planned |
| Annotation-informed gene factor model | gene | Gaussian | Bernoulli/probit or logit | gene-factor membership | none | Planned |
| Multi-view gene factor model | gene | view-specific Gaussian | shared factors | gene annotations | molecular anchors | Planned later |
| SNP-level factor model | SNP/fine-mapped unit | Gaussian or transformed evidence | sparse | SNP-factor membership | fixed \(S\) post hoc | Deferred until gene model validated |
| LD-aware SNP factor model | SNP | correlated residual/LD-aware | sparse | SNP-factor membership | fixed/uncertain \(S\) | Deferred research extension |
| Joint uncertain SNP-to-gene model | SNP and gene | joint hierarchical | sparse | SNP and gene priors | sampled \(S\) | Deferred |
| Pathway factor model | pathway | derived/posterior mapping | inherited | pathway annotations | \(C^\top Z\) | Post-processing first |

---

## 5. Data-operator compatibility

Legend:

- **Yes** — intended supported combination.
- **Reference** — current/legacy path used for validation.
- **Later** — planned after a reference implementation.
- **No** — invalid or intentionally out of scope.

| Model family | dense individual | BED | dense summary | CSR | block-eigen | evidence |
|---|---:|---:|---:|---:|---:|---:|
| ST BayesC | Later/reference | Yes | Reference | Yes | Yes/internal | No |
| ST BayesR | Later/reference | Yes | Reference | Yes | Yes/internal | No |
| ST BayesRC | Later | Yes | Later | Yes | Internal/partial | No |
| ST hierarchy | Later | Later | Later | First target | Later | No |
| Canonical MT BayesC | First/reference | Later | First/reference | Later | Later | No |
| MT mixture models | Later | Later | Later | Later | Later | No |
| Factor-analytic MT-BLR | First | Later | Later | Later | Later | No |
| Gene evidence factor | No | No | No | No | No | Yes |
| SNP evidence factor | No | No | No | optional LD metadata | optional | Yes |

---

## 6. Trait-design contracts

| Data/design case | Full effect covariance | Full residual covariance | Notes |
|---|---:|---:|---|
| Same individuals and same genotype matrix | Yes | Yes | Cross-trait phenotype products are identifiable |
| Same genotype matrix, different phenotypes with missingness | Conditional | Conditional | Requires explicit missingness/sufficient-statistic contract |
| Partially overlapping GWAS samples | Yes for prior | Only with overlap information | Must not infer residual covariance from marginal summaries alone |
| Independent GWAS samples | Yes for prior | Usually diagonal | Genetic covariance and residual covariance are distinct |
| Gene/SNP evidence matrix | Factor covariance | Trait residual covariance by observation model | Separate factor likelihood |

---

## 7. Hierarchy capability levels

| Capability | Version 1 | Later |
|---|---:|---:|
| Multiple categorical layers | Yes | Yes |
| One membership per marker per layer | Yes | Yes |
| Fixed layer multipliers | Yes | Yes |
| Learned conjugate layer multipliers | Yes | Yes |
| Weighted geometric normalization | Yes | Yes |
| Trait-specific multipliers | Optional after ST | Yes |
| Component-independent BayesRC hierarchy | Yes after BayesC | Yes |
| Component-specific hierarchy | No | Deferred |
| Overlapping incidence matrices | No | Yes through log-additive policy |
| Nonlinear annotation encoder | No | Research extension |
| Combination with `selection_s` | Planned | Yes |

---

## 8. Multivariate sharing modes

| Mode | State representation | Intended trait count | Comments |
|---|---|---:|---|
| `independent` | independent binary states or restricted patterns | small to moderate | Diagonal prior covariance or separate fits |
| `shared` | null versus all-traits-active | small to moderate | One shared marker inclusion |
| `patterns` | explicit user patterns | small | Avoid automatic \(2^T\) explosion |
| `factor` | sparse latent factors | moderate to large | Preferred scalable multivariate strategy |

The framework should never silently enumerate an impractical number of states.

---

## 9. Covariance policy matrix

| Policy | Prior/update | Small \(T\) | Large \(T\) | Initial status |
|---|---|---:|---:|---|
| Scalar | inverse-chi-square | Yes | N/A | Current |
| Diagonal | independent inverse-chi-square | Yes | Yes | Planned canonical option |
| Full inverse-Wishart | inverse-Wishart | Yes | No | Candidate after audit |
| Separation strategy | variance + correlation priors | Yes | Limited | Deferred until formally specified |
| Factor analytic | loading/shrinkage priors | Yes | Yes | Planned later |
| Fixed covariance | user supplied | Yes | Yes | Planned validation mode |

No covariance policy should be called simply `sampleB()` in the new core.

---

## 10. Public interface mapping

Existing public wrappers should map to resolved specifications.

| Public call | Resolved state | Probability | Scale | Operator |
|---|---|---|---|---|
| `stblr_csr(method="bayesc")` | binary | global_binary | unit/MAF | CSR |
| `stblr_bed(method="bayesc")` | binary | global_binary | unit | BED |
| `stblr_csr(method="bayesr")` | mixture | global_dirichlet | component/MAF | CSR |
| `stblr_bed(method="bayesr")` | mixture | global_dirichlet | component | BED |
| `stblr_csr_annot(annotation_model="sbayesrc")` | mixture | annot_probit_stick | component/MAF | CSR |
| `stblr_bed(method="bayesrc")` | mixture | annot_probit_stick | component | BED |
| future hierarchical wrapper | binary/mixture | compatible | hierarchy/composite | CSR first |
| future canonical MT wrapper | trait_pattern | global/annotation | compatible | dense/CSR |
| future factor wrapper | factor | factor_annot | factor | evidence or regression |

---

## 11. Output capability matrix

| Output namespace | BayesC | BayesR/RC | Hierarchy | MT patterns | Factor model |
|---|---:|---:|---:|---:|---:|
| `marker` | Yes | Yes | Yes | Yes | row loadings |
| `trace` | Yes | Yes | Yes | Yes | Yes |
| `state` | binary | component | binary/component | pattern | factor membership |
| `component` | No | Yes | optional | optional | No |
| `pattern` | No | No | No | Yes | No |
| `annotation` | optional | BayesRC | optional | optional | membership effects |
| `hierarchy` | No | No | Yes | Yes | optional |
| `trait` | variance | variance | variance | covariance | trait loadings |
| `residual` | variance | variance | variance | covariance | observation residual |
| `factor` | No | No | No | factor-analytic only | Yes |
| `chains` | optional | optional | optional | optional | required for alignment |

Common formatted aliases may remain at the top level, but namespaced output is authoritative.

---

## 12. Explicitly deferred combinations

The following should not be implemented during the initial framework migration:

- arbitrary overlapping hierarchy layers;
- component-specific hierarchy parameters for every BayesRC component;
- automatic \(K^T\) trait-component combinations;
- genome-wide LD-aware SNP factorization;
- joint sampling of uncertain SNP-to-gene links;
- nonlinear/deep annotation encoders;
- simultaneous migration of all legacy multivariate algorithms;
- performance optimization before a validated reference implementation.

---

## 13. Framework readiness checklist

A new model/data combination is supported only when:

- the state, probability, scale, covariance, and operator policies are identified;
- the data-design contract is valid;
- a reduction path exists;
- deterministic RNG tests pass;
- output namespaces are defined;
- simulation recovery tests pass;
- unaffected model regressions pass;
- documentation states current limitations.

---

## 14. Language-binding capability matrix

The statistical engine should support a current R binding and preserve the option of a future Python binding.

| Layer | R now | Python now | Python later | Core requirement |
|---|---:|---:|---:|---|
| Public package interface | Yes | No | Yes | Binding-specific |
| Input validation and friendly errors | Yes | No | Yes | Mostly binding-specific |
| Marker/trait name alignment | Yes | No | Yes | Shared semantics, separate adapters |
| Typed model specification | Internal | No | Shared | Ordinary C++ |
| Statistical sampler | Yes through R | No | Same sampler | Ordinary C++ only |
| Chain RNG | Yes | No | Same behavior | Standard C++ chain-owned RNG |
| Likelihood operators | Yes | No | Shared | Ordinary C++ |
| Typed result object | Internal target | No | Shared | Ordinary C++ |
| R raw-list conversion | Yes | N/A | N/A | Rcpp adapter |
| NumPy/dict conversion | No | No | Yes | Future pybind11 adapter |
| Progress display | R adapter | No | Python adapter | Generic callback or diagnostics |
| Exception translation | R error | No | Python exception | Standard C++ exceptions in core |
| Disk-backed LD/BED access | Yes | No | Shared | Language-neutral file formats |

### Required current-stage constraints

- Core model and sampler headers must not require Rcpp.
- Core stochastic code must not use `R::rnorm`, `R::rchisq`, `arma::randn`, or equivalent runtime-global RNG.
- Core results must not be created as `Rcpp::List`.
- Core errors must use standard C++ exceptions.
- Core dimensions must follow the canonical conventions in the implementation plan.
- No Python dependency is added at this stage.
- No stable public C++ ABI is promised at this stage.

### Future Python scope

A future Python interface may provide:

- NumPy input/output conversion;
- Python model-specification classes or dataclasses;
- pybind11 exception translation;
- Python-native progress and diagnostics;
- access to the same language-neutral CSR and BED resources.

It must not contain a second implementation of the statistical methods.
