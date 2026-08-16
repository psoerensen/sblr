# Model capability matrix

> **Status: historical and superseded for annotation/prior capabilities.**
> This table predates the qualified ST CSR/block-eigen log-variance routes and
> also reflects an executable capability resolver that is known to lag tested
> public MT BayesRC/SBayesRC routes. Do not use it as current annotation-model
> capability truth. Use
> [the audited annotation-prior capability matrix](annotation_prior_architecture_matrix.md)
> for annotation and prior architectures, and
> [the backend inventory](blr_backend_inventory.md) for the broader current
> ST/MT/operator route inventory. The table below is retained as a historical
> snapshot until Phase 1 replaces the duplicated capability declarations.

## Historical snapshot

| Public route | Supported models | Annotation policy | Fixed `maf_effect_s` | Sampled `maf_effect_s` |
|---|---|---|---|---|
| `stblr_bed()` | `bayesc`, `bayesr`, `bayesrc` | global or probit-stick BayesRC | model-dependent | supported where the scalar kernel implements it |
| `stblr_csr()` | `sbayesc`, `sbayesr` | global | supported | supported by validated scalar routes |
| `stblr_csr_annot()` | `sbayesc`, `sbayesrc` | fixed-marker, group, learned-logistic, probit-stick | supported where validated | supported only where implemented |
| `stblr_block_eigen()` | `sbayesc`, `sbayesr`, `sbayesrc` | global or probit-stick | supported; retained low rank is default and keeps priors in marker space | supported where the shared scalar kernel implements it |
| `mtblr_bed()` | `bayesc` (corrected Cheng complete-pattern model) | sampled Dirichlet activity-pattern mass | not implemented | unsupported |
| MT CSR | unavailable pending a corrected MT summary likelihood | unavailable | unavailable | unavailable |
| MT block eigen | unavailable pending a corrected MT summary likelihood | unavailable | unavailable | unavailable |

At this historical checkpoint, unsupported combinations were intended to fail
before native execution and `.blr_model_capability_matrix()` was intended
to own dispatch truth. The current resolver is itself part of the documented
capability drift and must not be treated as the sole authority until the
Phase-1 capability consolidation is completed.
