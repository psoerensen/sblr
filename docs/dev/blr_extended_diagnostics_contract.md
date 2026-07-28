# Extended BLR convergence diagnostics contract

## Purpose and authority

Phase 21 extends the shared scalar convergence layer without changing its
Phase 17U mathematics. Rank-normalized and folded split R-hat, bulk/tail/mean
ESS, posterior SD, MCSE of the mean, relative MCSE, threshold evaluation, and
warning aggregation have one R implementation. The generic bundle remains
`blr_convergence_trace_bundle`, schema version 1, with values ordered
iteration by chain by scalar quantity.

## Modes, tiers, and controls

`auto` is core-only and quiet with one chain; `none` captures nothing; `core`
explicitly requests `vbs`, `vgs`, `ves`, `vle`, and `vld`; `extended` adds the
applicable Tier 2 groups. Tier 3 requires an explicit, nonempty marker ID or
one-based-index vector. `NULL`, `"all"`, `"*"`, and logical shortcuts never
mean all markers.

The additive controls are `extended_groups`, `selected_markers`,
`selected_marker_quantities`, `full_probability_states`, `max_trace_gb`, and
`allow_large_traces`. The four group identifiers are `covariance`,
`probability`, `selection_s`, and `annotations`. Marker quantities are `b`,
`d`, and `component`.

## Completed-iteration checkpoint and ownership

Every formal trace is post-burn, unthinned, and read at the same completed
iteration checkpoint as posterior accumulation. Current state is repeated
between less-frequent parameter updates. All buffers are chain private; the
selection map and descriptor plan are immutable. Workers touch only plain C++
state. Capture is observational and cannot consume RNG, rebuild an operator,
rerun a chain, or depend on `nthin`, `keep_chains`, retained trace output, or
worker assignment.

## Inventory

| Parameter | Meaning/owner | State at checkpoint | Applicability | Physical storage |
|---|---|---|---|---|
| `vbs/vgs/ves/vle/vld` | existing ST/MT core | current iteration | scientifically defined core routes | double |
| strict-lower B | sampled base effect covariance | post-B update | MT; `not_updated` if fixed | double |
| strict-lower G | derived genetic covariance | current effective effects | MT | double |
| strict-lower E | residual covariance | post-E update | full-E MT; fixed is `not_updated` | double |
| diagonal-E off-diagonal | structural zero | no stochastic state | MT diagonal-E | no physical stochastic trace |
| BayesC active/pattern pi | global inclusion/pattern state | post-pi update | matching ST/MT models | double |
| BayesR component/pattern pi | global simplex or exact joint marginals | post-pi update | BayesR/SBayesR | double |
| full MT joint pi | global joint simplex | post-pi update | explicit BayesR opt-in | double |
| BayesRC pattern pi | conditional non-null pattern simplex | post-pi update | BayesRC/SBayesRC | double |
| selection S | sampled MAF-scale parameter | current MH state | sampled-S ST routes only | double |
| alpha | annotation coefficient | current scheduled-update state | sampled annotation routes | double |
| sigmaSqAlpha | annotation-effect variance | current scheduled-update state | sampled annotation routes | double |
| group parameters | group inclusion/variance state | current update state | sampled group routes | double |
| selected b | effective marker effect | current marker state | explicit selected markers | double |
| selected d | trait activity | current pattern/state | explicit selected markers | int32, numeric workspace |
| selected component | ordered mixture class | current component state | BayesR/BayesRC only | int32, numeric workspace |

Fixed marker-specific priors, `dm`, posterior component-probability matrices,
annotation-driven marker prior matrices, and processed annotations are not
default formal trace quantities.

## Covariance scalarization

Tier 1 owns diagonals. Tier 2 adds only strict-lower entries in R
column-major order: `(2,1), (3,1), (3,2), ...`. Names are
`cov_b[trait2,trait1]`, `cov_g[trait2,trait1]`, and
`cov_e[trait2,trait1]`; the descriptor stores the lower-numbered trait in
`trait_index` and higher-numbered trait in `trait2_index`. Raw covariances are
diagnosed, not correlations or Cholesky factors. Scalar rows make no claim
about joint positive-definite geometry.

## Probability semantics

Binary complements share one physical diagnostic key. ST BayesC stores active
mass only. Multicategory simplexes retain deterministic established order.
MT BayesR defaults to component and conditional-pattern marginals calculated
from the current joint probability state; the full joint simplex is explicit
opt-in. BayesRC diagnoses conditional pattern pi and its annotation parameters,
not marker-averaged annotation prior probabilities. These are marginal scalar
diagnostics and do not establish joint simplex convergence.

## Selection-S, annotation, and group semantics

Sampled S is diagnosed from its genuine task-private state. Fixed S is
`not_updated`; absent S is `not_applicable`; sampled MT S remains unsupported.
Coefficient matrices flatten annotation fastest within stick. Coefficient
names are `alpha[annotation,stick]`; variance names are
`sigmaSqAlpha[stick]`; intercept identity is explicit. Fixed alpha and fixed
variance states are `not_updated`. Group and learned policies diagnose only
their actual low-dimensional sampled states, never marker-expanded priors.

ST native ownership is route-specific but uses one adapter contract. CSR and
block-eigen SBayesR and packed-BED BayesR record the current global component
simplex in every trait-chain; binary null/active mixtures retain only the active
physical scalar. Group BayesC records group pi and variance multipliers,
learned-logistic BayesC records distinct inclusion and variance coefficients,
and SBayesRC/BayesRC records alpha and sigmaSqAlpha. Fixed-marker policies have
no low-dimensional sampled annotation parameter and therefore allocate none.
Every buffer is populated at the completed-iteration checkpoint, including
repeated current values between scheduled coefficient updates.

## Explicit selected-marker semantics

Character IDs resolve against final marker metadata; integers refer to the
final one-based marker order. Types cannot be mixed, duplicates and missing or
unknown values fail, and requested order is preserved. `b` is the effective
effect (zero while inactive), `d` is binary trait activity, and `component` is
the established ordered mixture code. No latent beta or redundant MT pattern
integer is captured. Rank diagnostics on component codes assess movement over
ordered variance classes, not the full categorical allocation distribution.

For ST, resolved selected indices are shared immutably across trait × chain
tasks. Each task copies only requested effective b, binary d, and (for
BayesR/BayesRC) component codes into task-private buffers. Fixed-marker,
ordinary/scheduled BayesC, BayesR, learned/group BayesC, and probit-stick routes
all use this direct indexed path where publicly applicable; component requests
for BayesC fail before execution.

## Descriptors, statuses, and output

Descriptors carry tier, group, parameter, trait pair, marker, component,
pattern, annotation, stick, intercept, update, derivation, structure, capture,
and stable diagnostic-key identity. One summary row represents one scalar.
Statuses are `computed`, `computed_fewer_than_four_chains`,
`computed_partial`, `not_updated`, `not_applicable`, `structural_zero`,
`constant`, `constant_chain_mismatch`, `unavailable_single_chain`,
`insufficient_draws`, `nonfinite`, and `not_requested`. Fixed or constant
states are never called perfectly converged.

`keep_traces` retains the generic bundle at `fit$convergence_traces`;
`keep_chains` independently retains compact logical-chain records. Group
overviews cover core, covariance, probability, selection S, annotations, and
selected markers. One main-thread convergence advisory summarizes all flags.

## Memory and safety

Before sampling the resolved plan accounts separately for numeric capture,
int32 state capture, numeric scalar workspace, retained numeric conversion,
descriptor metadata, and summary output. Workspace is O(chains × draws) because
the engine processes one scalar at a time. Genotype, phenotype, and operator
storage are excluded from diagnostic trace accounting. Above `max_trace_gb`,
execution fails without `allow_large_traces`; the override warns once and
never truncates, thins, drops, or reduces the request.

## Interpretation limits

Diagnostics are advisory chain-mixing summaries. Passing thresholds neither
proves absolute convergence nor model correctness. Selected-marker rows are
not automatic fine-mapping evidence. No ESS or MCSE for posterior SD, quantile
MCSE, median MCSE, matrix-valued diagnostic, or joint-simplex diagnostic is
introduced.
