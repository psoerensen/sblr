# Cleaned `sblr` example workflow suite

This folder contains a curated set of example workflows derived from the previous
example scripts:

- `basic_stblr_workflow.R`
- `sparse_ld_bed_workflow.R`
- `annotation_based_models.R`

The goal is to demonstrate package capabilities without mixing canonical examples,
heavy benchmarks, debugging snippets, and exploratory analyses in the same file.

## Files

| File | Purpose |
|---|---|
| `00_workflow_helpers.R` | Shared base-R helper functions used by the workflow scripts. |
| `01_basic_csr_workflow.R` | Canonical summary-statistics CSR workflow: summary stats, sparse LD, BayesC/BayesR, summaries, architecture, credible sets. |
| `02_sparse_ld_and_bed_workflow.R` | Parallel CSR and BED examples using `stblr_csr()` and `stblr_bed()`. |
| `03_annotation_models_workflow.R` | Annotation-unaware and annotation-aware CSR models: prior, learned, group, and SBayesRC. |
| `04_ld_swap_and_finemapping_workflow.R` | LD-swap diagnostics and credible-set/fine-mapping post-processing. |
| `05_selection_s_workflow.R` | Fixed and sampled `selection_s` examples for CSR BayesC/BayesR. |
| `06_multichain_diagnostics_workflow.R` | Single-chain and multi-chain diagnostics, chain summaries, and LD-swap chain summaries. |

## Running the workflows

Set the example data directory before running scripts that load example data:

```r
Sys.setenv(SBLR_EXAMPLE_DATA_DIR = "path/to/example/data")
```

The directory is expected to contain a file such as:

```text
Glist_sparseLD_1k.RDS
```

Heavy or benchmark-like sections are guarded by:

```r
Sys.setenv(SBLR_RUN_HEAVY_EXAMPLES = "true")
```

By default, those sections are skipped.

## Notes

These scripts are demonstration templates. The MCMC settings are intentionally short.
Real analyses need longer chains, convergence checks, posterior predictive checks, and
careful validation of LD, marker order, allele alignment, and phenotype preprocessing.
