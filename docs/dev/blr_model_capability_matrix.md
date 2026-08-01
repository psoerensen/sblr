# Model capability matrix

| Public route | Supported models | Annotation policy | Fixed `maf_effect_s` | Sampled `maf_effect_s` |
|---|---|---|---|---|
| `stblr_bed()` | `bayesc`, `bayesr`, `bayesrc` | global or probit-stick BayesRC | model-dependent | supported where the scalar kernel implements it |
| `stblr_csr()` | `sbayesc`, `sbayesr` | global | supported | supported by validated scalar routes |
| `stblr_csr_annot()` | `sbayesc`, `sbayesrc` | fixed-marker, group, learned-logistic, probit-stick | supported where validated | supported only where implemented |
| `stblr_block_eigen()` | `sbayesc`, `sbayesr`, `sbayesrc` | global or probit-stick | supported; retained low rank is default and keeps priors in marker space | supported where the shared scalar kernel implements it |
| `mtblr_bed()` | `bayesc`, `bayesr`, `bayesrc` | global or probit-stick BayesRC | supported for mixture models | unsupported |
| `mtblr_csr()` | `sbayesc`, `sbayesr`, `sbayesrc` | global or probit-stick BayesRC | supported for mixture models | unsupported |
| `mtblr_block_eigen()` | `sbayesc`, `sbayesr`, `sbayesrc` | global or probit-stick BayesRC | supported for mixture models | unsupported |

Unsupported combinations fail before native execution. The runtime matrix is
owned by `.blr_model_capability_matrix()` and tested against this public
contract.
