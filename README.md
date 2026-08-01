# sblr

`sblr` is an experimental R package for scalable Bayesian genomic prediction
and variable selection. It implements scalar-trait BLR (STBLR) and joint
small-trait-count BLR (MTBLR) samplers in R/C++ with deterministic multichain
execution and optional OpenMP chain dispatch.

The package is a research implementation and is not yet CRAN-ready. Its model
semantics, operator contracts, fit schema, and seven fitting interfaces are,
however, deliberately stabilized and tested.

## Installation

The package contains Rcpp, Armadillo, and OpenMP-aware native code. Install a
compiler toolchain compatible with your R release (Rtools on Windows), then:

```r
install.packages("remotes")
remotes::install_github("psoerensen/sblr")
```

R 3.5.0 or newer is declared; current development and validation use a modern
R toolchain. Some genotype workflows also use qgg-compatible PLINK BED
metadata.

## Canonical fitting interfaces

```r
stblr_csr()          # scalar-trait summary statistics + sparse LD
stblr_csr_annot()    # scalar-trait annotation-aware sparse-LD models
stblr_block_eigen()  # scalar summary statistics + retained low-rank block eigen
stblr_bed()          # scalar-trait individual-level packed BED

mtblr_csr()          # joint multi-trait summary statistics + sparse LD
mtblr_block_eigen()  # joint multi-trait summary statistics + block-eigen LD
mtblr_bed()          # joint multi-trait individual-level packed BED
```

Public model names encode data level:

| Model | Prior kernel | Data level |
|---|---|---|
| `bayesc` | BayesC | individual level |
| `sbayesc` | BayesC | summary statistics |
| `bayesr` | BayesR | individual level |
| `sbayesr` | BayesR | summary statistics |
| `bayesrc` | annotation-informed BayesR | individual level |
| `sbayesrc` | annotation-informed BayesR | summary statistics |

The `s` prefix means **summary statistics**. It never activates
`maf_effect_s`, which is an independent optional MAF-dependent effect-variance
scale.

## Model and operator support

| Interface | Operator | Supported models |
|---|---|---|
| `stblr_bed()` | packed BED | `bayesc`, `bayesr`, `bayesrc` |
| `stblr_csr()` | CSR sparse LD | `sbayesc`, `sbayesr` |
| `stblr_block_eigen()` | retained low-rank block eigen (default); reconstructed dense by request | `sbayesc`, `sbayesr`, `sbayesrc` |
| `stblr_csr_annot()` | CSR sparse LD | `sbayesc` with `fixed_marker`, `group`, or `learned_logistic`; `sbayesrc` with `annotation_probit_stick` |
| `mtblr_bed()` | packed BED | `bayesc`, `bayesr`, `bayesrc` |
| `mtblr_csr()` | CSR sparse LD | `sbayesc`, `sbayesr`, `sbayesrc` |
| `mtblr_block_eigen()` | block eigen | `sbayesc`, `sbayesr`, `sbayesrc` |

Unsupported combinations fail before numerical execution.

The retained low-rank operator follows the GCTB/SBayesRC eigenspace likelihood
strategy, represented in `sblr` cross-product units with a global projected
residual-variance contract. Retained rank, transformed scores, and marker
conditionals are crosswalked to the pinned GCTB strategy; `sblr` does not
reproduce GCTB's block-specific residual-variance procedure.

## Diagnostics

`convergence = "auto"` preserves the core five trait-level diagnostics when
multiple chains are available. `"core"` requests `vbs`, `vgs`, `ves`, `vle`,
and `vld`; `"extended"` can additionally diagnose covariance, probability,
sampled `maf_effect_s`, and annotation parameters. Explicitly selected marker
IDs or indices can retain mixing diagnostics for effective effects (`b`),
binary activity (`d`), and mixture components where defined. Diagnostic trace
capture is post-burn, unthinned, chain-private, RNG-neutral, and protected by
a pre-execution memory guard.

Diagnostics assess MCMC mixing; they do not establish model correctness or
absolute convergence.

## Documentation

- [GitHub Pages site](https://psoerensen.github.io/sblr/)
- [Practical Notes](https://psoerensen.github.io/sblr/notes/)
- [Statistical Methods](https://psoerensen.github.io/sblr/methods/)
- [Credible sets and fine-mapping](https://psoerensen.github.io/sblr/notes/credible_sets_and_finemapping.html)
- [Maintained workflows](examples/workflows/README.md)
- [Generated R help sources](man/) (`help(package = "sblr")` after installation)
- [Current implementation contracts](docs/dev/blr_architecture.md)

Start with [the workflow index](examples/workflows/README.md). The scripts use
small deterministic fixtures or clearly identify required external BED/LD
inputs; production analyses require longer runs and sensitivity checks.
