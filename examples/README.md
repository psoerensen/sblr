# Examples

## Maintained workflows

The supported, reviewed examples live in [`workflows/`](workflows/README.md).
They cover the canonical STBLR and MTBLR interfaces, all three operator
families, annotation models, `maf_effect_s`, convergence diagnostics, and
operator comparisons. These are the examples promised by the package.

## Exploratory scripts

Other scripts directly under `examples/` are research or development
experiments. They may require external data, local configuration, substantial
runtime, or internal functions, and are not stable public workflows:

- `Evaluation_mtblr_revised.R`
- `mtsim_prior_sparse_ld.R`
- `mtsim_sparse_ld.R`

## Benchmarks and tests

- Manual performance and memory studies: [`tools/benchmarks/`](../tools/benchmarks/)
- Permanent automated tests: [`tests/testthat/`](../tests/testthat/)

Do not use exploratory scripts as API documentation; use the maintained
workflow index and generated R help instead.
