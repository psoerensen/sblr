# Maintained BLR workflows

These scripts use the current seven-function API and canonical model names.
Source `workflow_helpers.R` when a script requests its deterministic fixture
helpers.

| Workflow | Current contract |
|---|---|
| [`st_summary_workflow.R`](st_summary_workflow.R) | `stblr_csr()` with `sbayesc`/`sbayesr` and sparse summary operators |
| [`st_individual_workflow.R`](st_individual_workflow.R) | `stblr_bed()` with `bayesc`/`bayesr`, plus matched summary comparisons |
| [`st_annotation_workflow.R`](st_annotation_workflow.R) | `stblr_csr_annot()` with all four canonical annotation policies |
| [`mt_models_workflow.R`](mt_models_workflow.R) | MT BayesC/BayesR across packed BED, CSR, and block eigen |
| [`mt_bayesrc_workflow.R`](mt_bayesrc_workflow.R) | MT BayesRC/SBayesRC annotation preprocessing and probability outputs |
| [`maf_effect_s_workflow.R`](maf_effect_s_workflow.R) | Independent fixed and supported sampled MAF-dependent scaling |
| [`convergence_workflow.R`](convergence_workflow.R) | Core/extended convergence, selected markers, retained traces, and memory preflight |
| [`operator_comparison_workflow.R`](operator_comparison_workflow.R) | Construction and comparison of CSR, retained low-rank block eigen, reconstructed-dense reference, and packed-BED routes |
| [`workflow_helpers.R`](workflow_helpers.R) | Small deterministic fixture and output-inspection helpers |

The public `s` model prefix denotes summary-statistics data. It is unrelated
to whether `maf_effect_s` is active. Demonstration MCMC lengths are intentionally
small and are not recommendations for scientific analyses.
