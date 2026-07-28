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

## Marker fields

| Field | Dimensions | Meaning and applicability |
|---|---|---|
| `bm` | marker × trait | posterior mean **effective** marker effect |
| `dm` | marker × trait | posterior probability of binary activity/non-null status; for MT patterns it is marginalized over pattern and component states |
| `b_final` | marker × trait | primary-chain final effective effect; exactly zero for inactive traits |
| `d_final` | marker × trait | primary-chain final binary activity, never a mixture component code |
| `beta_final` | marker × trait | final latent effect where the kernel exposes latent/effective separation |
| `component_final` | marker for MT, model-specific ST structure | final ordered mixture component; zero is null |
| `component_probabilities` | marker × component (or documented trait-specific ST collection) | posterior component-state probabilities, including null, with rows summing to one |

Posterior means and final mutable states are not interchangeable. Primary
chain 1 owns final state; posterior means pool retained draws across chains.
`dm` does not encode the component number.

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
For BayesR/SBayesR MT fits they are the named joint pattern × component
simplex. For BayesC they represent the applicable global inclusion or pattern
state. BayesRC/SBayesRC has marker-dependent component priors, so these common
fields remain `NULL`; explicit `pattern_pi_*` and annotation parameters live
under `model_parameters$annotations`.

`prior_component_probabilities` are annotation-driven prior probabilities;
`component_probabilities` are posterior allocation probabilities. They must
not be conflated.

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
