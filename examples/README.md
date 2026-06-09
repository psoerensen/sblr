# Development Examples

The scripts in this directory are exploratory workflows, development notebooks,
benchmarks, and manual testing utilities used while developing `sblr`. They are
not formal automated tests or consistently polished user-facing examples.

Some scripts may download external data, require local paths to be edited, run
long benchmarks, call experimental functions, or contain development-only
diagnostics. Review a script and its dependencies before running it.

## Clean Workflows

- `workflows/basic_sblr_summary_stats.R`: Concise template showing the intended
  high-level `sblr()` workflow using summary statistics and LD-derived matrices.

## Current Scripts

- `Evaluation_mtblr_revised.R`: Exploratory workflow and benchmark for
  multi-trait Bayesian linear regression and sparse linkage disequilibrium
  utilities.
- `Test_sparse_ld.R`: Manual development checks for PLINK BED readers, sparse
  linkage disequilibrium utilities, and sampler comparisons.
- `mtsim_sparse_ld.R`: Development notebook containing simulation workflows,
  candidate sparse-LD and BED wrapper functions, benchmarks, and recovery
  diagnostics.
- `mtsim_prior_sparse_ld.R`: Development notebook for annotation-informed
  priors, annotation-aware sampler wrappers, simulations, benchmarks, and
  diagnostics.

## Future Cleanup

- Reusable package functions should move into `R/` with documentation and
  focused tests.
- Small stable checks should move into `tests/testthat/` once the test suite is
  added.
- Polished user-facing workflows should move into `vignettes/` or concise
  examples in the package README.
- Obsolete experiments can later move into `examples/archive/`.

Automated checks in `tests/testthat/` should remain small, deterministic, and
independent of machine-specific paths or external data downloads.
