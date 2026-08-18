# sblr

**Scalable Bayesian Linear Regression**

`sblr` is an open-source R package for Bayesian genomic prediction, variable
selection, fine-mapping, and parameter estimation. It provides scalar-trait and
joint multi-trait models for individual-level genotypes and GWAS summary
statistics.

The package combines:

- BayesC, BayesR, and annotation-informed prior models;
- packed PLINK BED genotypes;
- sparse CSR LD matrices;
- retained low-rank and reconstructed block-eigen LD operators;
- deterministic multichain execution with optional OpenMP scheduling;
- validated posterior outputs, convergence diagnostics, and canonical
  extraction and summarization functions.

The statistical models, sampler transitions, operator contracts, fit schema,
and seven maintained fitting interfaces have been extensively checked through
independent reference calculations, mathematical reductions, serial/parallel
comparisons, and benchmark studies.

> **Scientific software**
>
> `sblr` is actively developed and distributed from GitHub. Its maintained
> interfaces and output contracts are documented and tested. Scientific
> analyses still require suitable chain lengths, convergence assessment,
> sensitivity analysis, and validation appropriate to the application.

## Installation

The package contains Rcpp, Armadillo, and OpenMP-aware native code. Install a
compiler toolchain compatible with your R release—Rtools on Windows—and then
install `sblr` from GitHub:

```r
install.packages("remotes")
remotes::install_github("psoerensen/sblr")
```

## What `sblr` supports

| Capability | Available implementation |
|---|---|
| Individual-level analysis | Packed PLINK BED genotypes |
| Summary-statistics analysis | Sparse CSR and block-eigen LD operators |
| Scalar-trait models | BayesC, BayesR, and annotation-informed models |
| Multi-trait models | MT-BayesC and MT-BayesR |
| Variable selection | Marker inclusion probabilities and mixture states |
| Fine-mapping | Posterior inclusion probabilities and credible-set helpers |
| Parameter estimation | Genetic, residual, effect, and applicable covariance quantities |
| Parallel execution | Deterministic chain-level OpenMP scheduling |
| Diagnostics | Convergence traces, chain summaries, seeds, workers, and provenance |

## Canonical fitting interfaces

```r
# Scalar-trait models
stblr_bed()
stblr_csr()
stblr_csr_annot()
stblr_block_eigen()

# Joint multi-trait models
mtblr_bed()
mtblr_csr()
mtblr_block_eigen()
```

The scalar-trait functions fit one trait at a time. The multi-trait functions
fit one joint model across traits within each logical chain.

## Model and data support

| Interface | Data representation | Maintained models |
|---|---|---|
| `stblr_bed()` | Individual-level packed BED | BayesC, BayesR, annotation-informed BayesR |
| `stblr_csr()` | Summary statistics with sparse CSR LD | SBayesC, SBayesR |
| `stblr_csr_annot()` | Summary statistics with sparse CSR LD and annotations | Fixed-marker, grouped, learned-logistic, log-variance, and annotation-informed mixture models |
| `stblr_block_eigen()` | Summary statistics with block-eigen LD | SBayesC, SBayesR, and supported annotation-informed models |
| `mtblr_bed()` | Common-sample packed BED | MT-BayesC and MT-BayesR |
| `mtblr_csr()` | Independent summary providers with sparse CSR LD | MT-BayesC and MT-BayesR |
| `mtblr_block_eigen()` | Independent summary providers with block-eigen LD | MT-BayesC and MT-BayesR |

Unsupported model, data, or operator combinations fail before numerical
sampling.

Current multi-trait summary-statistics models require providers declared
independent. Overlap-aware multi-trait likelihoods are not yet implemented.

## Model names

For scalar-trait models, the `s` prefix denotes a summary-statistics
likelihood:

| Individual-level model | Summary-statistics model |
|---|---|
| `bayesc` | `sbayesc` |
| `bayesr` | `sbayesr` |
| `bayesrc` | `sbayesrc` |

The `s` prefix does not activate MAF-dependent effect scaling.
`maf_effect_s` is a separate optional effect-variance control on interfaces
that explicitly support it.

The multi-trait interfaces use `method = "bayesc"` or `method = "bayesr"`.
Their data representation is determined by the selected interface:

- `mtblr_bed()` for individual-level genotypes;
- `mtblr_csr()` for sparse-LD summary statistics;
- `mtblr_block_eigen()` for block-eigen summary statistics.

MT-BayesC uses joint marker-activity patterns. MT-BayesR extends these joint
activity patterns with multiple effect-size scales. The detailed prior
factorization and covariance transitions are documented in the statistical
methods.

## Posterior access

Maintained scalar- and multi-trait fits use the same access pattern:

```r
fit <- stblr_bed(...)

effects <- extract_posterior(fit, "effects")
pips <- extract_posterior(fit, "pips")

effect_summary <- summarise_posterior(
  fit,
  quantity = "effects"
)

diagnostics <- extract_diagnostics(fit)
```

The responsibilities are intentionally separate:

- `extract_posterior()` retrieves stored posterior quantities without
  summarizing them;
- `summarise_posterior()` summarizes retained draws and chains;
- `extract_diagnostics()` returns sampler, execution, convergence, provider,
  and provenance information.

Unavailable scientific quantities return `NULL`; they are not reconstructed
from unrelated fields. Trait, marker, chain, and covariance dimensions are
preserved where applicable.

## LD representations

`sblr` supports three principal data representations:

- packed BED for individual-level analysis;
- sparse CSR LD for summary-statistics analysis;
- block-eigen LD for compressed or reconstructed summary-statistics analysis.

Operator provenance, marker alignment, approximation metadata, and retained
rank are recorded with the fitted object. The
[data-representation guide](https://psoerensen.github.io/sblr/notes/data_representations.html)
describes the corresponding likelihoods and alignment requirements.

## Diagnostics and reproducibility

Multichain execution uses deterministic task seeds and one logical task per
chain. When OpenMP is available, chains can be scheduled across workers
without changing the scientific result.

Depending on the fitted model, diagnostics can include:

- variance and covariance traces;
- activity- and scale-probability traces;
- marker-state diagnostics;
- rank-normalized R-hat;
- effective sample sizes and Monte Carlo standard errors;
- retained and convergence indices;
- task seeds, worker assignments, and provider provenance.

Diagnostics assess MCMC behavior. They do not by themselves establish model
correctness or guarantee convergence.

## Documentation

- [Documentation website](https://psoerensen.github.io/sblr/)
- [Practical notes](https://psoerensen.github.io/sblr/notes/)
- [Statistical methods](https://psoerensen.github.io/sblr/methods/)
- [Annotation-informed priors](https://psoerensen.github.io/sblr/methods/annotation_priors.html)
- [Data representations](https://psoerensen.github.io/sblr/notes/data_representations.html)
- [Credible sets and fine-mapping](https://psoerensen.github.io/sblr/notes/credible_sets_and_finemapping.html)
- [Maintained workflows](examples/workflows/README.md)
- [Benchmark and validation studies](https://psoerensen.github.io/sblrbench/)
- [Development contracts](docs/dev/blr_architecture.md)

Start with the
[maintained workflow examples](examples/workflows/README.md). They use small
deterministic fixtures or explicitly identify required external BED, summary
statistics, and LD inputs.
