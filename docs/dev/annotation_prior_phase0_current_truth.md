# Annotation-prior Phase 0 current implementation truth

**Status: CURRENT IMPLEMENTATION CHECKPOINT.** This document freezes observed
repository behavior at the baseline below. Sections labelled **CURRENT** report
source, dispatch, schema, formatter, and test evidence. Sections labelled
**PROPOSED FOR PHASE 1** are a gap list only; no proposed item is implemented by
this checkpoint.

## Baseline and authority

- Baseline commit: `fbd03b0891e66ca453fdcc32bfc98c4b9913e361`.
- Starting tracked worktree: clean.
- Raw schema versions: `stblr_raw` version 1 and `mtblr_raw` version 1.
- Formatted model-semantics version: 2,
  `s_prefix_means_summary_statistics`.

Implementation truth was resolved in this order:

1. current R/C++ source, executable dispatch, raw validators, and tests;
2. maintained developer contracts;
3. the approved annotation-prior audit, matrix, and implementation plan;
4. canonical methods theory;
5. qualification records.

Historical/research records and `docs/notes/` were not used as capability
authority. In particular, a registered native symbol was not counted as a
public capability unless a current public R route reaches it.

## Evidence inspected

The mandatory Phase 0/1 reading set from [the developer document authority
map](README.md) was read in full:

- `blr_development_guide.md`, `blr_architecture.md`,
  `blr_backend_inventory.md`, `blr_model_contracts.md`,
  `blr_output_schema.md`, `blr_convergence_contract.md`,
  `blr_test_ownership.md`, `blr_block_eigen_contract.md`,
  `stblr_low_rank_operator_design.md`, and
  `block_eigen_gctb_residual_contract.md`;
- `annotation_prior_architecture_audit.md`,
  `annotation_prior_architecture_matrix.md`, and
  `annotation_prior_architecture_implementation_plan.md`;
- `../methods/model_theory.qmd`, `../methods/annotation_priors.qmd`,
  `../methods/mt_bayesr_sbayesr.qmd`,
  `../methods/mt_bayesrc_sbayesrc.qmd`, and
  `../methods/sbayesrc_annotations.qmd`;
- `annotation_log_variance_implementation.md`.

Executable evidence included `NAMESPACE`, public fitter definitions, route-local
resolvers, annotation helpers and adapters, `stblr_raw`/`mtblr_raw` validators
and formatters, the BayesC/BayesR policy-aware cores, LV policies, retained
block execution, native registration, relevant model kernels, and the permanent
test files named below. Principal R sources were:

- `R/stblr-public.R`, `R/stblr-csr-annot.R`,
  `R/stblr-block-eigen.R`, `R/stblr-logvar-annotations.R`,
  `R/stblr-logvar-block-eigen.R`, `R/sparse_ld_bed_helper.R`;
- `R/annotation-helpers.R`, `R/stblr-csr-prior-annot.R`,
  `R/stblr-csr-learn-annot.R`, `R/stblr-csr-group-annot.R`,
  `R/stblr-csr-sbayesrc.R`, `R/stblr-bed-bayesrc-internal.R`,
  `R/sbayesrc-helpers.R`;
- `R/mtblr-csr.R`, `R/mtblr-bed.R`, `R/mtblr-block-eigen.R`,
  `R/mtblr-bayesr.R`, `R/mtblr-bayesrc.R`, and
  `R/mtblr-summary-chains.R`;
- `R/blr-unified.R`, `R/mtblr-convergence.R`, and `R/RcppExports.R`.

Principal native sources were the operator-aware BayesC/BayesR cores and
policies, the LV types/policies/kernel, the CSR and block LV bindings, the
ordinary CSR bindings, retained-block execution, and the MT joint-state and
BayesRC paths. These include `src/blr_csr_bayesc_core_impl.h`,
`src/blr_csr_bayesr_core_impl.h`, `src/blr_csr_bayesc_policy.h`,
`src/blr_csr_bayesr_policy.h`, `src/blr_csr_logvar_*`,
`src/st_logvar_annotation_prior.*`, `src/st_block_eigen_execution.h`,
`src/st_cpg_omp_csr*.cpp`, and the corresponding log-variance translation
units.

## Canonical public fitter inventory — CURRENT

`NAMESPACE` exports exactly the following seven canonical fitting interfaces.
Internal model adapters and native entry points are not exported fitters.

| Fitter | Family/operator | Accepted public models | Annotation/prior controls | Effect-scale controls | Dispatch and post-format behavior | Principal permanent tests |
|---|---|---|---|---|---|---|
| `stblr_csr()` | ST CSR summary statistics | `sbayesc`, `sbayesr` | none beyond global model controls | fixed or sampled `maf_effect_s`; aligned `effect_maf`; BayesC scheduled mode cannot sample S; BayesR is unscheduled | `.blr_resolve_st_model()` then `.stblr_csr_impl()` to ordinary BayesC/BayesR native routes; canonical ST finalizer and CSR contract attachment | `test-stblr-csr-interface.R`, `test-blr-operator-reductions.R`, reproducibility tests |
| `stblr_csr_annot()` | ST CSR summary statistics | `sbayesc` for fixed/group/learned; `sbayesrc` for probit stick; `sbayesc` or `sbayesr` for LV | `annotations`, `annotation_model`, policy-specific controls; LV has `theta_prior_sd`, `theta_init`, `updateTheta` | fixed/sampled S only with SBayesRC; rejected for fixed/group/learned/LV | route-local policy dispatch. Fixed/group/learned/SBayesRC use private R adapters. LV uses separate C++ translation units. All use the canonical ST finalizer, followed by route-specific metadata; LV identity and top-level logvar diagnostics are rewritten/restored after finalization | annotation backend/interface/chain/LD-swap tests; all `test-logvar-*.R` |
| `stblr_block_eigen()` | ST block-eigen summary statistics | `sbayesc`, `sbayesr`, `sbayesrc` | singular `annotation` for SBayesRC; `annotations` (or singular fallback) plus `annotation_model = "log_variance"` for LV | fixed or sampled S for ordinary C/R/RC; S rejected jointly with LV | route-local resolver; ordinary and LV BayesC/R share operator-aware engines and block construction; SBayesRC has its own route; canonical ST finalizer plus block metadata; LV identity/diagnostics are post-finalizer additions | `test-stblr-block-eigen.R`, `test-logvar-block-eigen.R`, residual-policy and operator-reduction tests |
| `stblr_bed()` | ST packed BED individual data | `bayesc`, `bayesr`, `bayesrc` | BayesRC annotations and probit-stick controls through `...` | MAF-S arguments are rejected as unsupported | route-local method check and packed-BED adapters; canonical ST finalizer | `test-stblr-bed-interface.R`, BED backend tests |
| `mtblr_csr()` | joint MT CSR summary statistics | `sbayesc`, `sbayesr`, `sbayesrc` | SBayesRC `annotations`, preprocessing and alpha controls | fixed S for BayesR/BayesRC; sampled S rejected; BayesC rejects S/component controls | `.mtblr_resolve_public_method()` then MT CSR native route; BayesRC raw enrichment and formatted `model_parameters$annotations` occur in route-specific R code | MT BayesR/BayesRC model, annotation, and operator tests |
| `mtblr_block_eigen()` | joint MT reconstructed block-eigen summary statistics | `sbayesc`, `sbayesr`, `sbayesrc` | same MT BayesRC controls | fixed S for BayesR/BayesRC; sampled S rejected | route-local resolver and MT block native route; BayesRC enrichment/formatting as above | MT BayesR/BayesRC operator tests |
| `mtblr_bed()` | joint MT packed BED individual data | `bayesc`, `bayesr`, `bayesrc` | same MT BayesRC controls | fixed S for BayesR/BayesRC; sampled S rejected | route-local resolver and MT BED native route; BayesRC enrichment/formatting as above | MT BED, BayesR, and BayesRC tests |

All seven expose the common chain/convergence controls. `stblr_csr_annot()` is
the only public route for fixed-marker, group, learned-logistic, and CSR LV.
There are no public top-level `stblr()` or `mtblr()` fitters.

## Executable capability truth — CURRENT

### Core family by route

Here “PUBLIC” means executable through one of the seven fitters and protected
by current route/model tests.

| Traits | Operator | BayesC | BayesR | BayesRC |
|---|---|---|---|---|
| ST | BED | PUBLIC as `bayesc` | PUBLIC as `bayesr` | PUBLIC as `bayesrc` |
| ST | CSR | PUBLIC as `sbayesc` | PUBLIC as `sbayesr` | PUBLIC as `sbayesrc` through `stblr_csr_annot()` |
| ST | block eigen | PUBLIC as `sbayesc` | PUBLIC as `sbayesr` | PUBLIC as `sbayesrc` |
| MT | BED | PUBLIC as `bayesc` | PUBLIC as `bayesr` | PUBLIC as `bayesrc` |
| MT | CSR | PUBLIC as `sbayesc` | PUBLIC as `sbayesr` | PUBLIC as `sbayesrc` |
| MT | block eigen | PUBLIC as `sbayesc` | PUBLIC as `sbayesr` | PUBLIC as `sbayesrc` |

The S prefix is a data-level name, not a different prior kernel and not an
implicit MAF-S switch.

### Annotation/prior mechanism by route

| Mechanism | ST BED | ST CSR | ST block eigen | MT BED | MT CSR | MT block eigen |
|---|---|---|---|---|---|---|
| fixed marker P/Q | NO_ROUTE | IMPLEMENTED_AND_PUBLIC, BayesC | NO_ROUTE | NO_ROUTE | NO_ROUTE | NO_ROUTE |
| group P/Q | NO_ROUTE | IMPLEMENTED_AND_PUBLIC, BayesC | NO_ROUTE | NO_ROUTE | NO_ROUTE | NO_ROUTE |
| learned logistic P/Q | NO_ROUTE | IMPLEMENTED_AND_PUBLIC, BayesC | NO_ROUTE | NO_ROUTE | NO_ROUTE | NO_ROUTE |
| probit-stick component P | IMPLEMENTED_AND_PUBLIC, BayesRC | IMPLEMENTED_AND_PUBLIC, SBayesRC | IMPLEMENTED_AND_PUBLIC, SBayesRC | IMPLEMENTED_AND_PUBLIC, BayesRC | IMPLEMENTED_AND_PUBLIC, SBayesRC | IMPLEMENTED_AND_PUBLIC, SBayesRC |
| BayesC-LV Q | EXPLICITLY_UNSUPPORTED | IMPLEMENTED_AND_PUBLIC | IMPLEMENTED_AND_PUBLIC | NO_ROUTE | NO_ROUTE | NO_ROUTE |
| BayesR-LV Q | EXPLICITLY_UNSUPPORTED | IMPLEMENTED_AND_PUBLIC | IMPLEMENTED_AND_PUBLIC | NO_ROUTE | NO_ROUTE | NO_ROUTE |

The public block fitter accepts both retained low-rank and historical dense
reconstructed representations for LV. Retained low rank is the canonical
scalable contract; both representations have current test coverage, while the
qualification emphasis is retained low rank.

No full model cell is merely internal-only: low-level R/native functions are
internal implementation details of the public cells above. The registered LV
math fixtures are internal testing hooks, not fitting routes.

### MAF-S

- ST CSR and ST block-eigen BayesC, BayesR, and SBayesRC expose fixed and
  sampled S in their ordinary routes. Scheduled CSR BayesC accepts fixed S but
  rejects sampled S; CSR BayesR is unscheduled.
- ST BED rejects MAF-S controls.
- fixed-marker, group, learned-logistic, and LV routes reject composition with
  S. ST CSR/block SBayesRC is the current annotation-aware route that accepts
  fixed or sampled S.
- MT BayesR/BayesRC accept fixed S on BED, CSR, and block eigen. MT BayesC
  rejects it. Sampled MT S is explicitly unsupported.

### Capability resolver discrepancy

There is no single canonical executable capability source. The effective
authority is route-local validation in each fitter. The table returned by
`.blr_model_capability_matrix()` is not consulted by public dispatch; its only
current callers are `tools/audit/blr_architecture_audit.R` and
`test-blr-unified-public-contract.R`.

That table has 46 rows, omits `log_variance`, and labels every MT
BayesRC/SBayesRC cell unsupported. Its permanent test explicitly expects those
stale values. This contradicts executable public MT BayesRC routes and their
passing operator tests. The discrepancy is a Phase 1 gap, not a Phase 0 fix.

## Current MT probability architecture

### Ordinary MT BayesR — CURRENT joint simplex

`.mtblr_bayesr_spec()` constructs one state table containing the null state and
every non-null trait-pattern by positive-component state. `joint_pi`, its
Dirichlet prior, and the native latent state all use this complete simplex.
Component and pattern probabilities are derived marginals. They are not
independent P and H factors. `test-mtblr-bayesr-model.R` fixes the state names,
component indices, multipliers, probabilities, and BayesC reduction; MT
operator tests exercise all three operators.

### MT BayesRC/SBayesRC — CURRENT P times H

`.mtblr_bayesrc_prior_probabilities()` computes marker-specific component
probabilities P from the probit sticks. `.mtblr_calibration_inputs()` forms each
non-null joint state probability as marker component probability times the
global conditional active-pattern probability `pattern_pi_init` (H); the null
uses the marker null probability. Native results retain marker component
probabilities and global pattern probabilities separately. MT BayesRC tests
verify row normalization, factorization, fixed-S reduction, annotation
alignment, and all three operators.

These current representations must not be conflated with the proposed future
provider vocabulary.

## Annotation preprocessing inventory — CURRENT

There is no common immutable annotation-design object. Current transforms are
provider- and route-specific.

| Path | Construction/alignment | Transform and validation | Stored provenance |
|---|---|---|---|
| fixed-marker BayesC, A + coefficients | `.stblr_prepare_annotation_matrix()`; reorders when row names contain all final marker names, otherwise accepts exact row count positionally; extra named rows can be dropped by subsetting | numeric finite; generated column names; optional intercept; continuous z-score by default; binary unchanged by default or z-scored when `center_binary_annotations = TRUE`; constants are not rejected by this generic helper; no duplicate/rank test | processed A and controls in `input`; wrapper adds top-level annotation/prior summaries |
| fixed direct marker P/Q | supplied trait lists are passed as marker-order vectors; native/R dimensional and numerical checks apply | P/Q bounds checked downstream; names are not a general marker-alignment contract | resolved vectors and use flags in `input`; raw `prior` namespace |
| group BayesC | `.stblr_prepare_group_index()` reorders a named group vector only when all marker names are present; otherwise exact-length positional | one disjoint factor assignment per marker; no overlapping A matrix; missing and empty groups rejected | factor, zero-based indices, names, sizes, priors, and initialization in `input`; raw `group` namespace |
| learned-logistic BayesC | same generic A helper as fixed A + coefficients | same optional intercept/scaling behavior; no duplicate/rank test; coefficient dimensions validated separately; induced P and Q are mean-centered then clamped | processed A and MH controls in `input`; raw `annotation`; formatter produces `eta_pi`/`eta_vb` and route code adds summaries |
| ST CSR SBayesRC | generic A helper | same generic transform profile; probit-stick initialization handled by `make_sbayesrc_alpha_init()` | processed A and alpha controls in `input`; raw `annotation`; formatter produces alpha/prior fields |
| ST BED BayesRC | `.stblr_align_bed_bayesrc_annotations()`; exact IDs when provided, otherwise exact selected-marker row count; factors expanded | rejects duplicate IDs/names, nonfinite and zero-variance non-intercepts, multiple intercepts; moves one intercept first; then generic standardization | explicit alignment and preprocessing records |
| ST block SBayesRC | SBayesRC block adapter using the current ST annotation helper profile | provider-specific current behavior, not the LV transform | current SBayesRC input/raw annotation metadata |
| MT BayesRC, all operators | `.mtblr_bayesrc_controls()` requires explicit marker IDs, calls the BED alignment helper, and rejects unused external rows | factor expansion and strict ID/name/finite/variance/intercept validation; binary unchanged unless requested; continuous standardized when requested | rich `metadata` and `model_parameters$annotations`, including alignment, transforms, names, and MAF-overlap flag |
| BayesC-LV/BayesR-LV, CSR and block | `.stblr_preprocess_logvar_annotations()` once in R; exact row count; non-default row names must exactly equal final marker order (no implicit reorder) | numeric finite; unique non-empty names; rejects all-ones intercept, constants, and rank deficiency/duplicates; binary center-only; continuous center and sample-SD scale; every X column is centered | `annotation_transform`, marker-alignment status in `input`, and top-level transform output |

The LV transform therefore has a materially stronger identification contract
than the shared historical helper. The strict MT BayesRC ID contract is also
different from ST CSR behavior. Current code does not mutate the caller's
annotation object intentionally, but it does not expose one common immutable
design/provenance structure.

## Provider and policy inventory — CURRENT

### Operator-neutral BayesC/BayesR policy surface

The ordinary BayesC and BayesR engines expose a narrow policy contract:

- optionally provide one marker `prior_scale` vector;
- run one `after_vb_update` hook with effect/state, current global marker
  variance, mixture multipliers for BayesR, and the fit-local RNG;
- capture, retain, and finish policy-owned summaries.

The ordinary no-op policies own no mutable state and do not draw from the RNG.
LV policies own theta, q, theta/q summaries and ESS diagnostics. They initialize
`q = exp(X theta)`, update theta by ESS in the post-vb hook when enabled, and
replace q after a successful update. The same engines are instantiated with CSR
and retained/dense block operators through binding-neutral block construction.

### Mechanisms not represented by that policy surface

- Fixed and sampled MAF-S are engine controls, not policy objects. The engine
  selects estimated S scale, fixed S scale, or policy scale; these scales are
  alternatives and are not composed.
- Fixed-marker, group, and learned-logistic behavior remains in separate
  BayesC CSR backends with mechanism-specific native state and updates.
- SBayesRC/BayesRC remains a separate component-probability kernel.
- MT marker scale is resolved as part of the MT mixture/specification path, not
  through the ST policy classes.

Thus only ordinary versus LV BayesC/R currently demonstrates a reusable
provider/policy seam across CSR and block operators. It is not yet a general
P/Q/H provider architecture.

## Raw-schema and formatter ownership — CURRENT

### Schema boundary

Active ST native backends return named `stblr_raw` version 1 objects. The ST
validator checks the schema tag, core dimensions, required pi state,
model-dependent component/annotation shapes, selected backend-specific fields,
and chain shapes. Active MT backends return named `mtblr_raw` version 1 objects;
its validator requires the core namespaces and has explicit BayesR/BayesRC,
block, BED, convergence, and annotation checks.

There is no active positional fallback, compatibility reader, or legacy schema
reinterpretation. Private historical positional formatter definitions remain
in `R/sparse_ld_bed_helper.R`, but no current public route calls them as a
fallback. Malformed active output stops. ST raw metadata does not itself carry
model-semantics version 2; the canonical finalizer records semantics in the
formatted `input` and `data`. MT raw producers record and the validator checks
version 2 when present.

### Field ownership map

| Quantity | Native/raw namespace | Family formatter | Route-specific work after family formatter/finalizer |
|---|---|---|---|
| marker means/final state/traces/variances/pi | ST `marker`, `trace`, `variance`, `pi`; analogous MT namespaces | `.as_stblr_fit()` / `.as_mtblr_fit()`, then `.blr_finalize_fit()` renames common pi/covariance/chain-summary fields | operator/data metadata only |
| posterior component allocation probabilities | ST `component$prob`; MT `marker$component_probabilities` | mapped to `component_probabilities`; ST formatter recomputes `dm = 1 - P(null)` | none for ordinary R; model summaries may be added later |
| ST SBayesRC alpha, sigma, marker prior | ST `annotation` | `.as_stblr_fit()` maps alpha/sigma and `marker_prior_final` to flat ST fields | `.standardize_stblr_annotation_fit()` and public wrapper add annotation effects, induced-prior summaries, and policy metadata |
| MT BayesRC alpha/sigma, P, H | MT `annotations`, enriched in R before validation/formatting | generic MT formatter does not own the final annotation layout | `.mtblr_bayesrc_format_fit()` copies them into `model_parameters$annotations`, derives prior-active probability, and installs explicit NULL common pi fields |
| fixed-marker | ST `prior` | ST formatter exposes `prior` | annotation matrix, `annotation_prior`, and summaries are route additions |
| group | ST `group` | ST formatter maps group means/counts | top-level annotation/group aliases and summaries are route additions |
| learned logistic | ST `annotation` | ST formatter maps `eta_pi` and `eta_vb` | annotation effects and summaries are route additions |
| LV theta/q | ST `annotation` contains theta, variance ratio, q and theta trace; `diagnostics$logvar` contains ESS diagnostics | `.as_stblr_fit()` has no LV mapping; `.stblr_attach_logvar_output()` runs after it and constructs all LV flat fields and theta diagnostics | canonical finalizer is called as ordinary `sbayesc`/`sbayesr`; public CSR/block wrappers then overwrite model identity and restore `diagnostics$logvar` at the top level |

LV diagnostics consequently occur both inside the preserved native diagnostic
payload (`diagnostics$native$logvar`) and as `diagnostics$logvar`. LV model
identity is not represented in `.blr_finalize_fit()`'s six-model enum. These are
ownership facts, not requests to change the current fit in Phase 0.

## Formatted output contract comparison — CURRENT

The following current contract statements are executable:

- canonical finalization owns `family`, `model`, `operator`, `input`, `data`,
  `diagnostics`, `convergence`, `convergence_traces`, `chains`, and
  `memory_estimate`;
- common formatted names include `pi_final`, `pi_mean`, `pi_trace`,
  `cov_b_mean`, `cov_g_mean`, `cov_e_mean`, `cov_b_final`, `cov_g_final`, and
  `cov_e_final` when their source fields exist;
- MT BayesRC common pi fields are present as NULL and P/H live under
  `model_parameters$annotations`;
- LV fields are the documented flat fields `theta`, `theta_summary`,
  `theta_chain_mean`, optional `theta_trace`, `annotation_variance_ratio`,
  `annotation_transform`, `marker_prior_scale`, and `diagnostics$logvar`.

Verified discrepancies with `blr_output_schema.md` are:

1. `model_parameters` is not made present by the canonical finalizer for every
   ordinary fit; it is model/route specific.
2. The implementation retains ST marker final-state names `b` and `d` (and
   `component` for ST mixture state), while the document describes common
   `b_final`, `d_final`, and `component_final`. MT explicitly constructs
   `component_final`; no current R formatter constructs `beta_final`.
3. `prior_component_probabilities` is not a common top-level field: MT BayesRC
   stores it under `model_parameters$annotations`; ST SBayesRC exposes its
   current flat/list annotation-prior representation.
4. LV scientific fields bypass family-formatter ownership and its model identity
   is installed after canonical finalization.

No proposed `fit$architecture` or `fit$model_spec$prior` structure exists.

## Numerical and RNG-neutrality baseline — CURRENT

No permanent Phase 0 test uses a digest/hash as its primary invariant. Current
protection is stronger where exact identity is expected: paths are compared
with `expect_identical()` or zero tolerance. Numeric fixture values reside in
the tracked current fixtures, not in historical reports.

| Invariant | Permanent owner and current expectation |
|---|---|
| ordinary BayesC/BayesR versus disabled LV policy | `test-logvar-bayesc.R`, `test-logvar-bayesr.R`, and `test-logvar-block-eigen.R`: exact same-seed trajectory identity at theta zero; fixed theta exactly matches fixed marker scale; q checked at about `1e-15` |
| LV math/ESS guards | `test-logvar-math.R`: eta/q about `1e-14`, likelihood about `1e-13`, fixed-q about `1e-15`, empty-active prior moments within `0.03`, 1D ESS moments within `0.04`, and explicit nonfinite/overflow failures |
| learned LV oracle | BayesC/R LV tests: MC-error-derived mean bounds, SD tolerance `0.12`/`0.15`, R-hat below `1.15`, induced-q correlation above `0.98`, plus pi/vb/ve tolerances |
| CSR/block LV and ordinary operator neutrality | `test-logvar-block-eigen.R` and `test-blr-operator-reductions.R`: exact theta/q and component paths where applicable; numeric operator tolerances range from exact/`1e-12` through `1e-7`, with explicitly approximate public comparisons up to `1e-4` |
| logical-task seed mapping and thread-count invariance | `test-blr-unified-reproducibility.R`: exact equality across core counts and chain-retention choices; explicit resolved task seeds |
| extended fixed fixtures | `test-blr-extended-reproducibility.R`: environment-gated current fixtures, generally `1e-12`; fixtures must not be regenerated during Phase 1 |
| raw-to-fit neutrality | `test-stblr-raw-schema.R`, `test-stblr-backend-field-inventory.R`, model semantics and public-contract tests |
| SBayesRC/MT architecture | SBayesRC backend/chain/LD-swap tests and MT BayesRC model/annotation/operator tests |

The no-op policy implementations are `noexcept`, own no state, and do not call
the supplied RNG. Exact zero-theta LV reductions are the permanent executable
proof that this remains true for the qualified BayesC/R engines.

## Permanent test ownership — CURRENT

| Phase 0 concern | Current owner | Ownership gap |
|---|---|---|
| public exports and model semantics | `test-blr-unified-public-contract.R`, `test-blr-model-semantics.R` | capability assertions are stale for MT RC and omit LV |
| effective route capability | route-specific ST/MT interface and operator tests | no single executable registry agrees with all routes |
| ST raw schema and formatter | `test-stblr-raw-schema.R`, `test-stblr-backend-field-inventory.R` | backend inventory does not include LV raw namespaces |
| MT raw schema and formatter | MT fitter/operator/BayesRC tests | no one cross-family schema inventory table |
| annotation preprocessing | annotation interface/backend tests, LV interface tests, MT BayesRC annotation tests | no shared cross-provider alignment/transform contract (because none exists) |
| LV scientific and numerical contract | all five `test-logvar-*.R` files | public learned-oracle tests do not make LV a canonical raw-schema namespace owner |
| operator neutrality | `test-blr-operator-reductions.R`, LV block tests | no unified provider-by-operator matrix owner |
| RNG neutrality | unified and extended reproducibility plus exact LV reductions | extended fixtures are opt-in; no generic provider no-op test interface |
| MT BayesR joint state | `test-mtblr-bayesr-model.R`, `test-mtblr-bayesr-operators.R` | current ownership is adequate |
| MT BayesRC P times H | `test-mtblr-bayesrc-model.R`, `test-mtblr-bayesrc-annotations.R`, `test-mtblr-bayesrc-operators.R` | current ownership is adequate |

## Focused baseline execution

All checks used existing code and fixtures; nothing was regenerated.

- `tools/audit/blr_architecture_audit.R`: 19/19 PASS.
- model semantics, public contract, ST raw schema: 171 PASS.
- LV math, public interface, block eigen: 128 PASS.
- LV BayesC/BayesR trajectory and oracle tests: 58 PASS.
- MT BayesR/BayesRC model and annotation tests: 69 PASS.
- ordinary operator reductions and reproducibility: 113 PASS.
- annotation interface and backend tests: 353 PASS.
- ST backend raw-field inventory: 369 PASS.
- MT BayesR/BayesRC operator tests: 203 PASS.

Focused test total: 1,464 PASS, 0 FAIL, 0 WARN, 0 SKIP. The only console
warning was that the installed `testthat` binary was built under R 4.4.3; it is
not a package-test warning. The full suite was not required or run for this
read-only checkpoint.

## Verified discrepancies

1. The standalone capability matrix and its test disagree with current public
   MT BayesRC routes and do not know LV.
2. Annotation alignment, intercept, binary/continuous scaling, constant/rank
   checks, and provenance are fragmented across at least the generic ST helper,
   strict BED/MT BayesRC helper, group helper, and LV helper.
3. The reusable BayesC/R policy seam covers LV only; historical fixed/group/
   learned/SBayesRC/MT scales have different owners.
4. MAF-S and a policy scale are exclusive branches, not composable providers.
5. Model-specific output ownership is split among native raw namespaces, family
   formatters, pre-finalizer annotation decorators, and post-finalizer public
   wrappers.
6. LV identity and flat fields are outside canonical formatter/model-enum
   ownership; top-level logvar diagnostics duplicate the native payload.
7. The maintained output document overstates universal `model_parameters` and
   uses final-state names not constructed by current formatters.
8. ST and MT raw validators differ substantially in strictness and namespace
   typing. This is current architecture, not evidence of invalid numerical
   behavior.
9. Dormant private positional formatter definitions remain in source, although
   active routes have no positional fallback.

## Exact Phase 1 gap list — PROPOSED FOR PHASE 1

Phase 1 must address architecture and ownership only after its own stop/go
review. The following list does not authorize a scientific change.

| Proposed item | Current problem and evidence | Intended owner | Behavior that must remain unchanged | Protecting tests/audits |
|---|---|---|---|---|
| one executable capability registry | route-local truth contradicts `.blr_model_capability_matrix()` | central R capability resolver consumed by all public dispatch and audits | every currently public/unsupported route above; early error wording may be deliberately migrated but support cannot drift | public contract plus all ST/MT interface/operator tests |
| immutable internal annotation design with transform profiles | four different alignment/transform contracts and incomplete provenance | one R preparation object; provider profile chooses scientific transform | each historical provider's transform until an explicitly approved model change; LV frozen transform exactly | annotation backend/interface, LV interface, MT BayesRC annotation tests |
| explicit current P/Q/H/provider descriptors | mechanisms are selected by model-specific wrappers and backend names | R model-spec resolution and binding-neutral descriptors, without combinatorial model enums | ordinary, fixed/group/learned, SBayesRC, LV, and MT state equations and defaults | model semantics, model-specific tests, trajectory fixtures |
| canonical model identity including LV | LV is finalized as ordinary then renamed | canonical resolver/finalizer identity mapping | public `sbayesc_logvar`/`sbayesr_logvar` identities and six existing names | LV interface/block tests and model-semantics tests |
| raw-schema ownership for all model-specific quantities | LV and MT BayesRC rely on route-specific decoration; ST extras have uneven validation | versioned named model namespaces and validators | existing numeric values, dimensions, optional trace behavior, and no positional fallback | raw-schema, backend inventory, LV/MT annotation tests |
| one canonical raw-to-fit ownership boundary | field attachment occurs before and after `.blr_finalize_fit()` | family formatter plus explicit typed model formatter stage | current formatted scientific values and present-NULL rules | raw schema, public contract, backend inventory, consistency checks |
| reconcile documented/current output names | maintained contract claims fields not constructed by code | output contract/schema owner, with an explicit later migration decision | no silent field rename in Phase 1; current `b`/`d`/ST `component` and model-specific layouts remain until approved | interface/output/consistency tests |
| normalize diagnostic ownership | LV diagnostics are duplicated at native and top-level paths | typed diagnostic namespace in canonical formatter | diagnostic values and availability; no new RNG or retained marker-by-iteration q history | LV diagnostics tests and convergence contracts |
| provider no-op and operator gates | only LV has a reusable policy seam and exact no-op reductions | common provider contracts where scientifically identical; operator remains orthogonal | zero extra RNG, exact ordinary trajectories, update order, LD swap, residual policies, MAF-S behavior | exact LV reductions, unified reproducibility, operator reductions, extended fixtures |
| test ownership repair | capability test freezes a false matrix and LV raw is absent from backend inventory | permanent capability/schema owner tests | no fixture regeneration and no weakened tolerance | `blr_test_ownership.md`, architecture audit, tests listed above |

The later informative-theta, external-q, MT-LV, BED-LV, composed P/Q, and new
public API phases are outside Phase 1 unless separately authorized by the
approved implementation plan.

## Invariants Phase 1 must preserve

1. Ordinary BayesC, BayesR, BayesRC/SBayesRC, fixed-marker, group,
   learned-logistic, MAF-S, LV, and MT joint-state mathematics and defaults.
2. Same logical-task seeds, RNG consumption, update order, marker order, LD
   swap, residual policy, and thread-count invariance.
3. Exact no-op provider behavior and current zero-theta/fixed-q reductions.
4. Current BED/CSR/retained-block numerical contracts and explicit dense
   reconstructed versus retained-factor distinction.
5. Current MT BayesR joint simplex and current MT BayesRC P times H
   factorization.
6. Named raw objects, schema validation, clear rejection of unsupported output,
   and no active positional fallback.
7. Current formatted scientific values and trace-retention/memory behavior,
   even if ownership is mechanically reorganized.
8. Frozen LV preprocessing, theta prior default 0.7, ESS, numerical guards, and
   empty-active prior draw.
9. Existing unsupported combinations, especially BED/MT LV, sampled MT S, and
   simultaneous learned probability plus LV variance.

Phase 0 ends with this checkpoint. It does not implement or authorize Phase 1.
