# Maintained workflows

These small scripts demonstrate the canonical APIs. They use simulated or
explicitly identified external inputs and are not substitutes for a complete
analysis plan.

| Script | Contract |
|---|---|
| `workflow_helpers.R` | shared input and compact reporting helpers |
| `st_summary_workflow.R` | ST summary-statistics CSR workflow |
| `st_individual_workflow.R` | ST individual-level packed-BED workflow |
| `st_annotation_workflow.R` | ST annotation policies |
| `mt_models_workflow.R` | MT BayesC/BayesR summary and individual routes |
| `mt_bayesrc_workflow.R` | MT BayesRC/SBayesRC annotations |
| `selection_s_workflow.R` | fixed and sampled supported `selection_s` cases |
| `convergence_workflow.R` | core/extended and selected-marker diagnostics |
| `operator_comparison_workflow.R` | CSR, block-eigen, and packed-BED reductions |

Heavy or external-data sections must be explicitly enabled and must never embed
machine-specific paths.
