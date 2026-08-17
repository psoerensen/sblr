# BLR output and raw schemas

## Formatted fit structure

Every canonical fit contains the common structural fields below. Structural
fields are present even when their value is `NULL`; model-specific scientific
fields are present only where defined.

| Field | Meaning |
|---|---|
| `family` | `stblr` or `mtblr` |
| `model` | canonical public data-level/model name |
| `operator` | `packed_bed`, `csr`, or `block_eigen` |
| `input` | resolved controls, prior kernel, probability policy, effect-scale policy, seeds, and convergence request |
| `data` | marker/trait metadata, sample sizes, data level, scales, alignment, and MAF provenance |
| `model_parameters` | model-specific mixture, annotation, group, and selection-state descriptions |
| `diagnostics` | execution, numerical, scheduler, annotation, timing, and advisory diagnostics |
| `convergence` | scalar diagnostic summary, descriptors, overview, and thresholds |
| `convergence_traces` | retained iteration × chain × quantity bundle, or `NULL` |
| `chains` | optional compact logical-chain records, independent of formal traces |
| `memory_estimate` | analytical ownership-based memory estimate, never measured RSS |

Marker and trait ordering exactly follow `data$marker_metadata` and
`data$trait_metadata`.

Scalar fits also record an operator contract in `input$operator_contract`.
The contract distinguishes an individual-level reference, a complete-CSR
summary-statistics reference, an explicitly approximate hard-sparse CSR, a
retained projected block-eigen likelihood, and an unclassified supplied CSR.
Construction metadata are evidence about how an operator was formed, not a
source-fidelity measurement. When source probes are unavailable, diagnostics
say so rather than manufacturing a pass threshold.

Retained low-rank raw diagnostics include the configured residual-rebuild
interval, rebuild counts, and maximum absolute and relative drift by logical
task. SBayesRC invalid-scale diagnostics are failure-only and opt-in; they add
context to the error without changing the raw or formatted scientific fields.
LD-swap summaries live at `diagnostics$ld_swap` and
`diagnostics$ld_swap_chains`; these entries remain present with `NULL`
values when unavailable or inapplicable.

## Marker fields

| Field | Dimensions | Meaning and applicability |
|---|---|---|
| `bm` | marker × trait | posterior mean **effective** marker effect |
| `dm` | marker × trait | posterior probability of binary activity/non-null status; for MT patterns it is marginalized over pattern and component states |
| `b_final` | marker × trait | primary-chain final effective effect; exactly zero for inactive traits |
| `d_final` | marker × trait | primary-chain final binary activity, never a mixture component code |
| `beta_final` | marker × trait | final latent effect where the kernel exposes latent/effective separation |
| `component_final` | marker for MT, model-specific ST structure | final ordered mixture component; zero is null |
| `component_probabilities` | marker × component (or documented trait-specific ST collection) | posterior component allocation; ordinary mixtures include null and sum to one, whereas factorized Phase 7 MT-BayesR reports positive-scale sub-probabilities that sum to the marker non-null probability |

Posterior means and final mutable states are not interchangeable. Primary
chain 1 owns final state; posterior means pool retained draws across chains.
`dm` does not encode the component number.

The internal Phase 4a qualification route returns `blr_raw` version 2 directly
and does not create a promoted formatted MT fit. It preserves draw, chain,
marker, trait, activity-pattern, observation, and covariance axes. Collapsed
null-marker latent effects are unavailable values on the fixed latent array
axes, while realised effects are exactly zero. `activity_pattern_parameters`,
markerwise `activity_pattern_probabilities`, traitwise `pips`, and
`pleiotropic_probabilities` have separate meanings. Fixed $V_e$ is recorded as
input and final state; no residual-covariance draws or posterior mean are
manufactured. See
[the Phase 4a checkpoint](blr_phase4a_cheng_mt_bayesc_checkpoint.md).

The Phase 4b sampled policy adds explicit
`draws$residual_covariance` (`draw × chain × trait_row × trait_col`),
unthinned convergence covariance, `posterior$residual_covariance_mean`, and
`final$residual_covariance` (`chain × trait_row × trait_col`). All are finite,
symmetric, and positive definite; the posterior mean is the retained-draw mean,
and the final state equals the last completed convergence state when traces are
kept. Fixed mode continues to require sampled draws and the posterior mean to be
present as `NULL`. See
[the Phase 4b checkpoint](blr_phase4b_sampled_residual_covariance_checkpoint.md).

`posterior$pleiotropic_probabilities` has one finite value in $[0,1]$ per
declared marker, in global marker order. It is conditionally required when the
resolved model declares a joint activity-pattern probability policy with an
identifiable all-traits-active pattern. Its values must equal the corresponding
column of `posterior$activity_pattern_probabilities`, aligned by the declared
activity-pattern ID rather than by column position, using a raw-object
consistency tolerance of $10^{-12}$. Phase 4a declares $(1,1)$ with ID `1_1`.
The Phase 5A general-$T$ qualification route declares the complete binary
matrix in `input$prior$probability$activity_patterns`; the unique
$(1,\ldots,1)$ row defines the required column without assuming a fixed ID or
column position. Trait and activity-pattern axes grow dynamically while
size-one draw and chain axes remain present. See
[the Phase 5A checkpoint](blr_phase5a_general_t_cheng_mt_checkpoint.md).
Phase 5B preserves this raw-v2 contract unchanged as the scientific output of
public `mtblr_bed()`. The formatted fit is constructed from these fields, and
the validated raw object is retained as `attr(fit, "blr_raw")`. Fixed $V_e$
does not manufacture covariance draws or a posterior covariance mean. See
[the Phase 5B checkpoint](blr_phase5b_public_cheng_mt_checkpoint.md).
The mandatory markerwise activity-pattern posterior has marker by
activity-pattern axes and therefore requires $M2^T$ numeric values. The current
qualification object also retains the scientifically equivalent joint-state
view with its separately declared axis, and preflight counts both material R
matrices conservatively. Activity-transition diagnostics are not posterior
fields: the dense $2^T\times2^T$ matrix is absent, `transition_counts` is
present with value `NULL`, and compact per-chain occupancy vectors use the
declared activity-pattern IDs. Under `compact_occupancy_v1`, occupancy uses
exact chain and pattern axes and contains finite nonnegative integer counts;
`pattern_change_counts` contains one finite nonnegative integer per named
chain. Occupancy totals equal the number of marker updates across all completed
burn-in and sampling sweeps, and changes cannot exceed that total.
The field is not required for single-trait, independent-trait, or other models
for which pleiotropy is scientifically undefined.

The Phase 6A independent-summary route, public through `mtblr_csr()` and
`mtblr_block_eigen()` since Phase 6B, uses the same general-$T$
activity, effect, probability, covariance, task, and compact-diagnostic axes.
Its resolved input retains every provider, local-to-global map, fixed
$\phi_p$, and operator approximation descriptor. Because marginal independent
summary providers do not identify individual fitted values or a full residual
covariance, `derived$predictions`, `draws$residual_covariance`,
`posterior$residual_covariance_mean`, and `final$residual_covariance` remain
present with value `NULL`. No operator-relative quadratic is relabelled as SSE
or genomic covariance. See
[the Phase 6A checkpoint](blr_phase6a_summary_mt_checkpoint.md).
Public promotion does not change these fields; see
[the Phase 6B checkpoint](blr_phase6b_public_summary_mt_checkpoint.md).

Phase 7 factorized MT-BayesR preserves activity-pattern assignments and
probabilities and adds separate positive-scale component assignments,
marker-by-component sub-probabilities
$\Pr(\boldsymbol\delta_j\ne\mathbf 0,k_j=k\mid D)$, and sampled
$\boldsymbol\omega$ states. The unique null pattern uses component assignment
`-1`; positive scales use zero-based indices in declared scale order. No
marker-by-pattern-by-scale posterior array is required. Raw-v2 validation uses
the declared binary activity metadata to require each marker's component
sub-probabilities to sum to one minus its unique null-pattern probability. It
also validates compact pattern and scale occupancy/change diagnostics under the
Phase 7 policy, while the named dense transition diagnostic remains `NULL`. See
[the Phase 7 checkpoint](blr_phase7_pattern_scale_mt_checkpoint.md).

## Variance and covariance fields

`vbs`, `vgs`, `ves`, `vle`, and `vld` are iteration × trait traces of marker
effect variance, total genetic variance, residual variance, diagonal/LE
genetic contribution, and LD contribution. At every recorded iteration,
`vld = vgs - vle` within numerical tolerance.

MT fits additionally expose full trait × trait posterior-mean matrices
`cov_b_mean`, `cov_g_mean`, and `cov_e_mean`, and primary-chain final matrices
`cov_b_final`, `cov_g_final`, and `cov_e_final`. Their diagonals correspond to
the scientific quantities represented by the trait-level traces, but the
matrices are not trace arrays.

## Probability fields

`pi_final`, `pi_mean`, and `pi_trace` are respectively the final global/model
probability state, posterior mean state, and per-iteration global/model trace.
For factorized Phase 7 MT-BayesR, `pi_*` describes the activity-pattern
simplex and `omega_*` separately describes the positive-scale simplex
conditional on a non-null pattern. For BayesC, `pi_*` represents the applicable
global inclusion or pattern state. BayesRC/SBayesRC has marker-dependent
component priors, so these common
fields remain `NULL`; explicit `pattern_pi_*` and annotation parameters live
under `model_parameters$annotations`.

`prior_component_probabilities` are annotation-driven prior probabilities;
`component_probabilities` are posterior allocation probabilities. They must
not be conflated. Factorized MT-BayesR positive-scale entries are joint
sub-probabilities with non-null activity, rather than a null-inclusive simplex.

## Current log-variance annotation fields

Qualified ST CSR and retained block-eigen BayesC-LV/BayesR-LV fits currently
add `theta`, `theta_summary`, `theta_chain_mean`,
`annotation_variance_ratio`, `annotation_transform`, and
`marker_prior_scale`. `theta_trace` is retained when `keep_chains = TRUE`;
formal convergence-trace retention remains governed separately by the
convergence controls. ESS numerical diagnostics live under
`diagnostics$logvar`. These fields describe the existing flat formatted
contract.

The nested `fit$architecture` and `fit$model_spec$prior` layout discussed
in the annotation architecture audit is **Proposed** and is not a current raw
or formatted schema. See
[the annotation-prior implementation plan](annotation_prior_architecture_implementation_plan.md).

## Chain aggregation

Posterior marker/covariance/probability means pool retained samples across
logical chains. Public displayed traces are iterationwise chain means. Compact
chains retain deterministic task summaries only when `keep_chains = TRUE`.
Formal trace retention is controlled separately by
`convergence_control$keep_traces`. Chain-mean stability fields use explicit
`*_chain_mean_sd`, `*_chain_mean_min`, and `*_chain_mean_max` names.

## Raw schemas

Both `stblr_raw` and `mtblr_raw` use raw schema version 1 and named namespaces;
new fits carry `model_semantics_version = 2` and
`model_semantics = "s_prefix_means_summary_statistics"`. Raw producers are
validated before one family-specific formatter creates the canonical fit.
There is no positional fallback, compatibility reader, or automatic
reinterpretation of ambiguous older objects.

## `blr_raw` version 2 target and Phase 1 implementation

Native backends still return the current version-1 objects described above.
Phase 1 implements `blr_raw` schema version 2 in R, following Sections 26 and
28 of
[the unified framework design](blr_unified_framework_design.md), backed by
`tests/research/blr_framework_contract/`. Eligible one-chain maintained ST
fits now follow the explicit path `stblr_raw` v1 to validated `blr_raw` v2 to
formatted fit. The validated v2 object is retained as the `blr_raw` attribute
of the formatted fit. Native output and trajectories are unchanged.

The exact envelope is:

```text
schema, model, input, posterior, draws, final, derived, diagnostics, provenance
```

Raw v2 keeps `draw`, `chain`, `marker`, `trait`, state/component, region, and
covariance-row/covariance-column axes even when their size is one. Canonical
effect axes are `draw × chain × marker × trait`; covariance axes are
`draw × chain × trait_row × trait_col`; final objects omit only `draw`.
Independent-trait states and joint-state indices are separate fields, never
rank-varying versions of one field. Model-inapplicable contracted fields remain
present with `NULL`.

Every present scientific array is validated against IDs declared by the
resolved specification: counts, order, uniqueness, nonmissing dimnames, and
axis names must agree for draw, chain, marker, trait, state, joint-state,
component, activity-pattern, region, covariance-row, covariance-column, and
declared observation/provider axes. A compatibility identifier cannot bypass
axis, probability-simplex, covariance, or provenance validation.

The provenance namespace keeps `git_sha`, `dirty_build`, `compiler`, and
`timestamp` present with `NULL` when unavailable. A present SHA is a validated
Git identifier and `dirty_build` is a nonmissing logical scalar. Development
provenance is cached once per session; reliable build-time Git injection
remains future promotion work.

Raw v2 retires unqualified `pi`, `pis`, `pim`, `state_probabilities`, and
`pattern_probabilities`. It distinguishes sampled probability parameters,
markerwise posterior state/component/activity probabilities, and PIPs. A
formatted compatibility alias may view exactly one documented field; no
universal `fit$pis` is manufactured.

The native scientific state, stable raw-v2 R object, and formatted fit are
separate layers. Formatting may simplify axes explicitly and create tables or
read-only aliases. It must not replace sampled covariance with a descriptive
estimate or convert the current MT covariance hybrid into sampled schema-v2
$V_b$. Serialized v1 objects are converted only when semantics and dimensions
are established by schema/compatibility metadata; otherwise conversion fails
and a scientifically affected analysis must be rerun.

### Phase 1 conversion boundary

The initial production converter covers one-chain ST BayesC, BayesR, and
BayesRC-family raw semantics reached through maintained CSR, BED, block-eigen,
and supported annotation wrappers. Independent traits preserve their trait
axis. Multi-chain ST results remain explicitly on v1 because the native v1
object does not expose the final effect vector for every chain; Phase 1 does
not fabricate those states. Current MT v1 objects are rejected because their
hybrid covariance cannot be labelled as sampled schema-v2 $V_b$.

The temporary formatted aliases have exactly these sources:

| Alias | Schema-v2 source |
|---|---|
| `bm`, `dm` | `posterior$realised_effect_mean`, `posterior$pips` |
| `b`, `d` | `final$realised_effects`, `final$independent_trait_states` |
| `vbs`, `vgs`, `ves`, `vle`, `vld` | named fields under `derived$legacy_iteration_quantities` |
| `pi_trace` | `derived$legacy_iteration_quantities$non_null_probability_parameter` |
| `pi_mean`, `pi_final` | family-specific explicit probability-parameter mean/final fields |
| `component_probabilities` | `posterior$traitwise_component_assignment_probabilities` |
| `cov_b_*`, `cov_g_*`, `cov_e_*` | corresponding explicit scalar variance posterior/final fields |

These aliases are read-only views in meaning. In particular, `fit$pis` is not
manufactured as a universal schema-v2 probability field.

### Phase 3 execution metadata

Eligible newly created ST fits record `seed_contract_version = 1`,
`retention_contract_version = 1`, and `scheduler_version = 1` in the resolved
input and raw-v2 provenance. They also record canonical logical task IDs,
exact final uint32 task seeds, exact retained post-burn transition indices,
and unthinned convergence iteration indices. Native worker diagnostics keep
requested cores, configured workers, actual sampler-region team size, and one
zero-based worker ID per canonical task as distinct fields.

These fields describe execution and do not change posterior field meanings.
Scheduled CSR, log-variance, group, BayesRC/SBayesRC, multi-chain raw-v2 gates,
and current MT results remain explicitly version 0 or legacy as applicable.
Current MT covariance is never promoted to `draws$marker_covariance`.
