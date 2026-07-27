# Unified BLR architecture

## Canonical dimensions

The package has two statistical families and four machine operator identifiers:

| dimension | identifiers |
|---|---|
| family | `stblr`, `mtblr` |
| model | `bayesc`, `sbayesc`, `bayesr`, `sbayesr`, `bayesrc`, `sbayesrc` |
| operator | `csr`, `block_eigen`, `packed_bed`, `dense_reference` |

STBLR fits traits independently; one logical task is `(trait, chain)`. MTBLR
fits a joint small-T model; one logical task is a complete joint chain. The
families share execution, convergence, output, memory, and warning principles,
not one numerical sampler.

## Canonical scientific-model and policy matrix

The six public models combine prior kernel and data level. Non-S names are
individual-level packed-BED models; S names are summary-statistics CSR or
block-eigen models. `selection_s` is a separate optional effect-scale policy,
so no summary model implicitly enables MAF scaling and no numerical kernel is
copied for the public prefix.

| family | model | CSR | block eigen | packed BED | probability policy | effect scale |
|---|---|---|---|---|---|---|
| stblr | `bayesc` | unsupported | unsupported | public_canonical | `global` | `unit` or `maf_s` when supported/requested |
| stblr | `sbayesc` | public_canonical | public_canonical | unsupported | `global` | `unit` or `maf_s` when supported/requested |
| stblr | `bayesr` | unsupported | unsupported | public_canonical | `global` | `component` or `component_maf_s` when supported/requested |
| stblr | `sbayesr` | public_canonical | public_canonical | unsupported | `global` | `component` or `component_maf_s` when supported/requested |
| stblr | `bayesrc` | unsupported | unsupported | public_canonical | `annotation_probit_stick` | `component` |
| stblr | `sbayesrc` | public_supported | public_canonical | unsupported | `annotation_probit_stick` | `component` or `component_maf_s` when supported/requested |
| mtblr | `bayesc` | unsupported | unsupported | public_canonical | `global` | `unit` |
| mtblr | `sbayesc` | public_canonical | public_canonical | unsupported | `global` | `unit` |
| mtblr | `bayesr` | unsupported | unsupported | public_canonical | `global` | `component` or `component_maf_s` when requested |
| mtblr | `sbayesr` | public_canonical | public_canonical | unsupported | `global` | `component` or `component_maf_s` when requested |
| mtblr | `bayesrc` | unsupported | unsupported | unsupported | `annotation_probit_stick` | `component` |
| mtblr | `sbayesrc` | unsupported | unsupported | unsupported | `annotation_probit_stick` | `component_maf_s` |

Annotation-aware BayesC remains a policy dimension of `stblr_csr_annot()`, not
three extra scientific models:

| model | probability policy | CSR | block eigen | packed BED |
|---|---|---|---|---|
| `sbayesc` | `fixed_marker` | public_supported | unsupported | unsupported |
| `sbayesc` | `learned_logistic` | public_supported | unsupported | unsupported |
| `sbayesc` | `group` | public_supported | unsupported | unsupported |

`bayesrc` and `sbayesrc` use `annotation_probit_stick`. Unsupported matrix
cells fail during public validation before numerical execution; they never
fall back to ordinary BayesC or a non-S model.

## Route inventory before Phase 18

| family | model | operator | public entry before Phase 18 | native entry | status before | topology | controls and traces | Phase 18 disposition |
|---|---|---|---|---|---|---|---|---|
| stblr | bayesc | csr | `stblr_csr` | `stblr_cpg_omp_csr` | public_canonical | trait x chain | multichain; compact traces optional | retain and align |
| stblr | bayesc | csr | `stblr_csr` scheduled | `stblr_cpg_omp_csr_scheduled` | public_supported | trait x chain | task-private post-burn `vbs`/`vgs`/`ves`/`vle`/`vld` traces | retain through the canonical CSR entry |
| stblr | bayesr | csr | `stblr_csr`, `stblr_csr_bayesr` | `stblr_cpg_omp_csr_bayesr` | public_supported | trait x chain | duplicate wrapper | retain only main entry |
| stblr | bayesc | csr | `stblr_csr_prior_annot` | `stblr_cpg_omp_csr_prior` | public_supported | trait x chain | fixed marker priors | internal dispatch via annotation entry |
| stblr | bayesc | csr | `stblr_csr_learn_annot` | `stblr_cpg_omp_csr_annot` | public_supported | trait x chain | learned annotations | internal dispatch via annotation entry |
| stblr | bayesc | csr | `stblr_csr_group_annot` | `stblr_cpg_omp_csr_group_annot` | public_supported | trait x chain | group model | internal dispatch via annotation entry |
| stblr | sbayesrc | csr | `stblr_csr_sbayesrc_generic` | `stblr_cpg_omp_csr_sbayesrc` | public_supported | trait x chain | annotation mixture | internal dispatch via annotation entry |
| stblr | bayesc | block_eigen | none | `stblr_cpg_omp_csr_block_eigen` | internal_canonical | trait x chain | validated block view | expose through `stblr_block_eigen` |
| stblr | bayesr | block_eigen | none | `stblr_cpg_omp_csr_bayesr_block_eigen` | internal_canonical | trait x chain | validated block view | expose through `stblr_block_eigen` |
| stblr | sbayesrc | block_eigen | none | `stblr_cpg_omp_csr_sbayesrc_block_eigen` | internal_canonical | trait x chain | validated block view | expose through `stblr_block_eigen` |
| stblr | bayesc | packed_bed | `stblr_bed` | scheduled-chain BED BayesC | public_canonical | trait x chain | multichain | retain and align |
| stblr | bayesr | packed_bed | `stblr_bed` | scheduled-chain BED BayesR | public_canonical | trait x chain | multichain | retain and align |
| stblr | bayesrc | packed_bed | `stblr_bed` | scheduled-chain BED BayesRC | public_canonical | trait x chain | multichain | retain and align |
| stblr | bayesc | packed_bed | `stblr_bed_marker` | sparse/scheduled low-level routes | experimental | backend-specific | overlapping public contract | unexport; retain internal reference |
| mtblr | bayesc | csr | `mtblr_csr` | `mtblr_csr_chains_raw_internal` | public_canonical | joint chain | shared CSR preparation; deterministic native multichain | aligned |
| mtblr | bayesc | block_eigen | `mtblr_block_eigen` | `mtblr_block_eigen_chains_raw_internal` | public_canonical | joint chain | shared block reconstruction; deterministic native multichain | aligned |
| mtblr | bayesc | packed_bed | `mtblr_bed` | `mtblr_bed_chains_internal` | public_canonical | joint chain | complete multichain and Tier 1 convergence | reference contract |
| mtblr | bayesc | dense_reference | `sblr` | `mtblr` | public_supported | one joint chain | legacy positional interface | unexport; retain internal_reference |
| mtblr | research | csr | none | `mtblr_cpg_omp_csr` | research_only | native research route | unsupported output contract | retain unexported |
| mtblr | research | dense_reference | none | `mtblr_eigen` | research_only | native research route | unsupported public contract | retain unexported |

All canonical raw adapters use schema version 1 and family-specific canonical
formatters. Numerical kernels, scientific fixtures, and reference routes remain
protected. Public status terms are exactly `public_canonical`,
`public_supported`, `internal_canonical`, `internal_reference`, `experimental`,
`research_only`, and `retired`.

For every `public_canonical` row, the public chain contract is `nchains >= 1`,
`ncores` requests logical-task workers, explicit `chain_seeds` preserve order,
`keep_chains` controls only compact records, convergence has
`auto`/`none`/`core`, raw schema remains version 1, and formatting ends in the
family canonical formatter plus `.blr_finalize_fit()`. Seeds are fit-local and
worker-independent. Memory warnings are pre-execution analytical bounds and
convergence warnings are aggregated main-thread advisories. Operator-specific
limitations are: CSR requires aligned standardized summaries and disk-backed
LD; block eigen requires BED provenance and records reconstruction/filter
policy; packed BED requires complete aligned phenotypes. Internal reference,
experimental, and research-only rows do not claim this public contract.

## Shared execution contract

Seeds attach to logical tasks before dispatch. Each task owns RNG and mutable
state. Static worker assignment cannot affect results. Workers create no R or
Rcpp objects, failures retain logical-task identity, and aggregation never
returns partial output after a task failure. `ncores` requests concurrent MCMC
tasks; `nthreads` is reserved for preparation, decoding, and construction.

## Shared output contract

Every supported fit has `family`, `model`, `operator`, `input`, `data`,
`diagnostics`, `convergence`, `convergence_traces`, `chains`, and
`memory_estimate`. Marker summaries and variance traces remain direct fields
when meaningful. Metadata has one owner: controls in `input`, data/provenance in
`data`, backend diagnostics in `diagnostics`, and analytical estimates in
`memory_estimate`.

## Shared convergence contract

One dependency-free scalar engine owns rank-normalized split/folded R-hat,
bulk/tail/mean ESS, posterior SD, and mean MCSE. Adapters supply unpooled,
post-burn, unthinned logical-chain traces. Compact-chain retention and
diagnostic-trace retention are independent.

The family-neutral core quantities are `vbs[trait]`, `vgs[trait]`,
`ves[trait]`, `vle[trait]`, and `vld[trait]`. Their public posterior traces use
iteration × trait orientation; convergence arrays use iteration × chain ×
quantity. For every recorded iteration, `vld = vgs - vle` within numerical
tolerance. MT `vbs`/`vgs`/`ves` are covariance diagonals, while the full
matrices remain separately owned by `cov_*_mean` and `cov_*_final`.

Statuses are `computed`, `computed_fewer_than_four_chains`,
`computed_partial`, `not_updated`, `not_applicable`, `structural_zero`,
`constant`, `constant_chain_mismatch`, `unavailable_single_chain`,
`insufficient_draws`, `nonfinite`, and `not_requested`.

## Warning and memory contract

Operator data are prepared as shared immutable state; task results are returned
in deterministic task/result order. Memory is an analytical upper-bound split
into operator storage, private sampler state per worker, result state per
logical chain, convergence capture/workspace, retained chains/traces, and
formatted output. It is never described as measured RSS. Memory, OpenMP
fallback, and convergence warnings are separate. Convergence uses one
main-thread aggregated advisory and no per-quantity warning flood.

Block-eigen runtime storage is float-packed reconstructed dense blocks. It is
not low-rank storage, even when filtering removes eigen-directions.

## Phase 19 MT pattern-by-component layer

All MT operators accept `bayesc`, `bayesr`, and `sbayesr`. The latter two use
one binding-neutral joint-state descriptor: one null state followed by each
supplied non-null trait pattern and its ascending positive mixture components.
CSR and block eigen share the summary core; packed BED retains the sample-space
core. The descriptor, probability normalization, scale removal for base-B
updates, chain dispatch, aggregation, convergence, and output contracts are
shared without a runtime-polymorphic operator hierarchy.

Each operator is prepared once per fit. Logical chains own component states,
joint probabilities, RNG, effects, residuals, accumulators, and core variance
traces. `keep_chains` controls compact records, while convergence trace capture
and retention remain independent. BayesR/SBayesR retain the common iteration ×
trait `vbs`, `vgs`, `ves`, `vle`, and `vld` traces.
## Phase 20 MT annotation-informed routes

The canonical MT matrix now includes `mtblr_bed(method="bayesrc")` and
`mtblr_csr()`/`mtblr_block_eigen(method="sbayesrc")`. All three routes share
`prior_kernel="bayesrc"`, `annotation_policy="annotation_probit_stick"`, one
Phase 19 joint-state descriptor, one base B, and the same chain/convergence
infrastructure. Annotation preprocessing and operator preparation each occur
once per fit; chain-private coefficient states never cross worker boundaries.
