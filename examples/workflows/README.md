# Maintained BLR workflows

These scripts use the current public fitting API and canonical model names.
Source `workflow_helpers.R` when a script requests its deterministic fixture
helpers.

| Workflow | Current contract |
|---|---|
| [`st_summary_workflow.R`](st_summary_workflow.R) | `stblr_csr()` with `sbayesc`/`sbayesr` and sparse summary operators |
| [`st_individual_workflow.R`](st_individual_workflow.R) | `stblr_bed()` with `bayesc`/`bayesr`, plus matched summary comparisons |
| [`st_annotation_workflow.R`](st_annotation_workflow.R) | `stblr_csr_annot()` with all four canonical annotation policies |
| [`mt_models_workflow.R`](mt_models_workflow.R) | Corrected Cheng MT-BayesC-Pi with common-sample packed BED data |
| [`maf_effect_s_workflow.R`](maf_effect_s_workflow.R) | Independent fixed and supported sampled MAF-dependent scaling |
| [`operator_comparison_workflow.R`](operator_comparison_workflow.R) | Construction and comparison of CSR, retained low-rank block eigen, reconstructed-dense reference, and packed-BED routes |
| [`workflow_helpers.R`](workflow_helpers.R) | Small deterministic fixture and output-inspection helpers |

The public `s` model prefix denotes summary-statistics data. It is unrelated
to whether `maf_effect_s` is active. Demonstration MCMC lengths are intentionally
small and are not recommendations for scientific analyses.
