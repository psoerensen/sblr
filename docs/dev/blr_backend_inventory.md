# BLR backend inventory

| Family | Operator | Public models | Native disposition |
|---|---|---|---|
| ST | packed BED | `bayesc`, `bayesr`, `bayesrc` | canonical sample-space kernels |
| ST | CSR | `sbayesc`, `sbayesr`; fixed-marker/group/learned-logistic BayesC, probit-stick SBayesRC, and BayesC-LV/BayesR-LV through `stblr_csr_annot()` | canonical summary kernels; LV has separate policy-aware scientific backends |
| ST | block eigen | `sbayesc`, `sbayesr`, `sbayesrc`; BayesC-LV/BayesR-LV through `annotation_model = "log_variance"` | canonical retained-factor kernels; reconstructed blocks are explicit references |
| MT | packed BED | `bayesc`, `bayesr`, `bayesrc` | canonical joint sample-space kernels |
| MT | CSR | `sbayesc`, `sbayesr`, `sbayesrc` | canonical joint summary kernels |
| MT | block eigen | `sbayesc`, `sbayesr`, `sbayesrc` | canonical joint reconstructed-block kernels |

Dense and low-level routes are internal references only when a permanent
scientific reduction owns them. Uncalled and unregistered experimental routes
are retired rather than retained for hypothetical future use.

Current LV support is ST CSR and retained block eigen only. BED-LV and MT-LV
are not implemented. For the complete current/proposed annotation-provider
matrix, see
[the annotation-prior architecture matrix](annotation_prior_architecture_matrix.md).
