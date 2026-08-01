# BLR backend inventory

| Family | Operator | Public models | Native disposition |
|---|---|---|---|
| ST | packed BED | `bayesc`, `bayesr`, `bayesrc` | canonical sample-space kernels |
| ST | CSR | `sbayesc`, `sbayesr`; annotation policies through `stblr_csr_annot()` | canonical summary kernels |
| ST | block eigen | `sbayesc`, `sbayesr`, `sbayesrc` | canonical retained-factor kernels; reconstructed blocks are explicit references |
| MT | packed BED | `bayesc`, `bayesr`, `bayesrc` | canonical joint sample-space kernels |
| MT | CSR | `sbayesc`, `sbayesr`, `sbayesrc` | canonical joint summary kernels |
| MT | block eigen | `sbayesc`, `sbayesr`, `sbayesrc` | canonical joint reconstructed-block kernels |

Dense and low-level routes are internal references only when a permanent
scientific reduction owns them. Uncalled and unregistered experimental routes
are retired rather than retained for hypothetical future use.
