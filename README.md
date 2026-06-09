# sblr

`sblr` stands for **Scalable Bayesian Linear Regression**.

## Overview

`sblr` is an experimental R package for Bayesian genomic prediction and
variable selection. It combines R interfaces with C++/Rcpp implementations of
scalable Bayesian linear regression models.

This repository is a research snapshot. It is not yet production-ready,
CRAN-ready, or intended as a stable analysis platform.

## Current status

- The package is experimental and APIs may change.
- Scripts in `examples/` are exploratory development workflows, not formal
  automated tests.
- Documentation, tests, and the intended public API are still being developed.
- Some workflows require external data, local path configuration, and long
  runtimes.

## Main ideas

The package explores:

- Bayesian linear regression for genomic prediction and variable selection
- Sparse mixture priors, including BayesC, BayesR, and SBayesRC-style models
- Models based on GWAS summary statistics and sparse linkage disequilibrium
  matrices
- Workflows using PLINK BED genotype files
- Annotation-informed prior models
- C++/Rcpp MCMC samplers with OpenMP parallelization

## Typical workflow

```text
Genotypes / phenotypes
  -> summary statistics and LD
  -> R wrapper initializes priors and MCMC state
  -> C++ Gibbs sampler
  -> posterior summaries
```

Posterior summaries can include marker effects, posterior inclusion
probabilities, variance components, and covariance estimates.

## Installation

Installation from GitHub requires an R package compilation toolchain because
the package contains C++ code. On Windows, install a version of Rtools
compatible with your R version.

```r
install.packages("remotes")
remotes::install_github("psoerensen/sblr")
```

The package is under active research development, so installation or workflows
may require additional setup.

## Examples

See [`examples/README.md`](examples/README.md) before running scripts in
`examples/`. These scripts include development notebooks, benchmarks, and
manual testing workflows, and some are machine-specific or computationally
expensive.
