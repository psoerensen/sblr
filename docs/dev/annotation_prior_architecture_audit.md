# Annotation-informed prior architecture audit

## Provenance

- Audit baseline branch: `master`
- Audit baseline commit: `2123699a9cc2e91059e7d81a745420b14eca7f6e`
- Package version: `0.2.0`
- Starting tracked worktree status: clean (`git status --short` produced no
  entries)
- Authoritative brief supplied for this audit:
  `local_reference/sblr_annotation_prior_architecture_audit_brief_v2.md`
- Audit date: 2026-08-10

The supplied v2 file is authoritative for this task. Its internal title and one
startup reference still say v1, while its appended section 27 contains the v2
refinements. This is a documentation provenance inconsistency, not a statistical
ambiguity.

This is a read-only scientific/code audit plus documentation reorganization.
No sampler, native binding, R function signature, default, output field, RNG
path, or generated file was changed.

## Scope and method

Claims were cross-checked across public R formals, dispatch, capability
resolution, preprocessing, prior calibration, native kernels, raw builders,
schema validation, canonical formatting, fit consistency checks, tests,
examples, methods documentation, developer records, and frozen local reference
material. The existence of a native symbol was not treated as evidence of a
supported production route.

### Principal files inspected

The audit included, in full or in the relevant implementation sections:

- governance: `AGENTS.md` and the v2 local audit brief;
- R public/dispatch: `R/blr-unified.R`, `R/stblr-csr-annot.R`,
  `R/stblr-csr-prior-annot.R`, `R/stblr-csr-group-annot.R`,
  `R/stblr-csr-learn-annot.R`, `R/stblr-logvar-block-eigen.R`,
  `R/stblr-csr-sbayesrc.R`, `R/stblr-bed-bayesrc-internal.R`,
  `R/sparse_ld_bed_helper.R`, `R/mtblr-csr.R`, `R/mtblr-bed.R`,
  `R/mtblr-block-eigen.R`, and `R/mtblr-bayesrc.R`;
- R preparation/calibration/output: `R/annotation-helpers.R`,
  `R/stblr-logvar-annotations.R`, `R/sbayesrc-helpers.R`,
  `R/stblr-maf-effect-s.R`, `R/prior-calibration.R`, raw-schema and
  raw-to-fit formatter code, consistency checks, and diagnostic helpers;
- native BayesC/R engines and policies: CSR core/type headers, ordinary and LV
  translation units, binding-neutral execution/policy headers, retained-block
  execution, LD-operator code, raw builders, and chain infrastructure;
- native annotation mechanisms: fixed-marker, group, learned-logistic,
  probit-stick BayesRC/SBayesRC, LV/ESS, MAF-S, MT joint-state and covariance-update
  code;
- methods: `docs/methods/annotation_informed_models.qmd`,
  `annotation_priors.qmd`, `sbayesrc_annotations.qmd`,
  `mt_bayesrc_sbayesrc.qmd`, and `mt_bayesr_sbayesr.qmd`;
- developer material: `blr_architecture.md`, `blr_backend_inventory.md`,
  `blr_output_schema.md`, `blr_model_capability_matrix.md`,
  `blr_prior_calibration_audit.md`, `blr_model_contracts.md`,
  `blr_convergence_contract.md`, LV implementation records, and SBayesRC
  development/qualification records;
- practical notes and examples under `docs/notes/`, including model overview,
  annotation, ST/MT, and MAF-S pages;
- frozen LV theory and executable reference scripts under `local_reference/`,
  plus relevant BayesRC/SBayesRC design/reference material;
- relevant tests under `tests/testthat/` covering annotation preparation,
  native math, fixed/group/logistic policies, BayesRC/SBayesRC, LV, MAF-S,
  prior calibration, capabilities, raw schemas, formatters, convergence,
  operator parity, trajectory identity, and RNG neutrality.

Generated `RcppExports` files were inspected as registration evidence only and
were not edited.

## Executive conclusions

1. **Probability, variance, sharing, and covariance are the right scientific
   decomposition.** The existing kernels already implement these ideas, but
   model-specific wrappers hide them.
2. **Operator and prior architecture should be orthogonal.** BED, CSR, and
   retained eigen are likelihood/operator implementations. `S`-prefixed names
   remain useful scientific display aliases, not different prior families.
3. **`annotation_model` should not survive as the long-term abstraction.** It
   forces mutually exclusive mechanisms even when the scientific quantities
   are orthogonal. Replace it with a typed prior specification after an
   additive migration.
4. **LV should become the canonical learned annotation-variance engine.** It
   has a frozen identifiable model, ESS, numerical diagnostics, ordinary-model
   trajectory gates, and CSR/eigen qualification.
5. **The variance half of `learned_logistic` should eventually be retired, not
   maintained independently.** Its probability half remains a useful logistic
   probability provider. Retirement is an intentional model change because
   the old q clipping, preprocessing, MH kernel, and update schedule differ.
6. **Fixed-marker q and future external q are the same scientific quantity.**
   They should share one variance provider and the public name
   `prior_variance_multiplier`.
7. **Group is a convenience specification plus a real hierarchical prior, not
   a separate core family.** Preserve its Beta and variance updates while
   expressing its outputs through probability/variance providers.
8. **The first MT-LV model should be MT BayesR/SBayesR with its existing
   global joint-state simplex and one shared q across traits:**
   \(\operatorname{Cov}(\beta_{j,p,k})=\gamma_kq_jV_{b,p}\). Do not start with
   BayesRC/SBayesRC, MT BayesC-LV, or trait-specific q.
9. **A common annotation-design object is needed, but a single silent transform
   is not.** Intercept and binary-centering policies are part of each model and
   historical paths currently differ.
10. **Raw-schema and formatter ownership lag behind sampler architecture.** LV
    and several annotation outputs are decorated after canonical formatting;
    ST and MT encode similar science differently. Typed raw namespaces should
    be the first structural cleanup.
11. **The capability resolver is materially stale.** It omits LV and reports MT
    BayesRC/SBayesRC unsupported despite public implementations and tests.
12. **SBayesRC plus LV must remain a research-gated composition.** Technical
    composability is not sufficient given alpha/occupancy mixing and likely
    P/Q confounding.

## Current mathematical architecture

### Core single-trait families

BayesC uses

\[
d_j\sim\operatorname{Bernoulli}(\pi_j),\qquad
b_j\mid d_j=1\sim N(0,v_bq_j).
\]

BayesR/BayesRC use

\[
c_j\sim\operatorname{Categorical}(\pi_{j0},\ldots,\pi_{j,K-1}),
\]

\[
b_j\mid c_j=k>0\sim N(0,v_b\gamma_kq_j),\qquad \gamma_0=0.
\]

The core mixture ladder \(\gamma\), global scale \(v_b\), marker probability
architecture \(P_j\), and marker variance architecture \(q_j\) are distinct.
Current model names frequently combine them.

The architecture uses lowercase v for scalar ST variances and uppercase V for
MT variance-covariance matrices:

| Quantity | ST | MT |
|---|---|---|
| marker effects | \(v_b\) | \(V_b\) |
| genomic | \(v_g\) | \(V_g\) |
| residual | \(v_e\) | \(V_e\) |

This is theory notation, not an R naming proposal. Existing code fields remain
`vb`/`covb`, `vg`/`covg`, and `ve`/`cove`; case-only R names are deliberately
avoided. Trait-specific scalar `v_le` and `v_ld` summaries remain lowercase.

### Current multi-trait effect model

For active trait pattern p and non-null component k, the implemented MT
mixture model is approximately

\[
\operatorname{Cov}(b_{j,p,k})=\gamma_kq_jV_{b,p}.
\]

Here \(V_{b,p}\) is the active-trait submatrix of \(V_b\), \(\gamma_k\) is the
BayesR relative component multiplier, q is the marker-specific relative
variance multiplier, and \(V_b\) carries the global active-effect covariance
scale.

Current MT BayesRC makes component probability annotation-dependent but retains
global conditional pattern probabilities. q may reflect supported fixed marker
scale such as MAF-S, but learned LV q is not implemented. \(V_b\) is global and
annotation-independent.

The current ordinary MT BayesR/SBayesR probability model is one `joint_pi`
simplex over the unique null state and all non-null
pattern-by-positive-component states. Its component and pattern probabilities
are marginals of that simplex; it is not currently parameterized as independent
\(P\) and \(H\) factors.

Current MT BayesRC/SBayesRC does use a factorization. Write
\(P(c_j=0\mid A_j)=P_{j0}\) and, for non-null component
\(k>0\) and non-null trait pattern p,

\[
P(a_j=p,c_j=k\mid A_j)=P_{jk}H_p.
\]

P controls null and non-null BayesR component probabilities; H controls the
conditional non-null pattern probability. Current MT BayesRC makes P
annotation-dependent, while H and \(V_b\) remain global. For MT BayesC, the single
non-null component reduces this to
\(P(d_j=0\mid A_j)=1-\pi_j\) and
\(P(a_j=p,d_j=1\mid A_j)=\pi_jH_p\).
That MT BayesC expression is an architectural special case, not a claim that
all current MT BayesC routes expose independent P/H providers.

### Orthogonal decomposition

The recommended internal model is the product of:

- a **core family**: BayesC or BayesR mixture ladder;
- a **probability provider P**: global, fixed marker, group, logistic, or
  probit-stick;
- a **variance provider Q**: unit, fixed external, group, LV, MAF-S, or a
  qualified composition;
- an MT **sharing provider H**: initially global pattern probabilities;
- an MT **covariance architecture \(V_b\)**: initially one global positive-definite
  covariance;
- a likelihood **operator**: BED, CSR, or retained block eigen.

This decomposition matches existing scientific update equations better than
the current `annotation_model` enumeration.

## Detailed current mechanism audit

### `fixed_marker`

#### Model and construction

The ST CSR BayesC route accepts fixed marker inclusion probabilities and/or
fixed relative variance multipliers. Either can be supplied directly or
constructed from A and fixed coefficients.

For probability, the historical implementation mean-centers the annotation
linear predictor, adds the global logit, applies `plogis`, then clamps to
configured limits (defaults approximately \(10^{-8}\) and 0.5). For variance,
it mean-centers the log-linear predictor, exponentiates, and clamps to configured
limits (defaults approximately \(10^{-3}\) and \(10^3\)). The pre-clamp q has
geometric mean one; clipping can change it.

Direct q must be positive. Marker updates, the \(v_b\) update
\(\sum_{j:d_j>0}b_j^2/q_j\), and scale-aware LD-swap density ratios use the
resolved multiplier. h2 calibration uses resolved marker probabilities and
variance multipliers.

#### Preprocessing

The common older annotation helper can reorder rows when complete marker row
names are supplied; otherwise it trusts row order if dimensions match. It
checks numeric finiteness and can add an intercept. Continuous columns are
standardized; binary columns are historically left uncentered by default.
Constant, duplicate, and rank-deficient columns are not rejected as strictly as
LV.

#### Output and tests

Current formatted output exposes `annotation_prior`, annotation metadata, and
summaries. Tests cover direct and A-derived values, h2 calibration, LD-swap
scale behavior, raw naming, multichain formatting, and RNG behavior. Support is
route-limited to ST CSR BayesC/SBayesC.

#### Conclusion

Keep the scientific fixed P/Q capability, but retire `fixed_marker` as a model
identity. Its q path should become the fixed form of the canonical variance
provider. Do not import its clipping defaults into LV.

### `group`

Markers belong to exactly one disjoint group. Each group may have a learned
inclusion probability with a Beta prior. Group variance multipliers can be
fixed or updated using a scaled-inverse-chi-square-style conditional based on
active marker effects and group-specific prior degrees/scale. Optional
normalization uses a marker-count-weighted arithmetic mean, not the LV
geometric-mean contract.

Outputs include group pi, variance multiplier, active counts, group sizes, and
extended traces/diagnostics. Tests cover fixed and sampled scales, update
controls, h2 calibration, chains, and raw fields. It is ST CSR BayesC only.

The hierarchy is scientifically distinct, but the sampler should eventually
consume ordinary group-valued P/Q providers. A group convenience constructor
can preserve the Beta and variance hyperpriors without a separate genomic
likelihood core.

### `learned_logistic`

#### Exact implemented model

For probability,

\[
\operatorname{logit}(\pi_j)=\operatorname{logit}(\pi)+
\operatorname{center}(A_j\eta_\pi).
\]

For variance,

\[
q_j=\exp\{\operatorname{clamp}[
\operatorname{center}(A_j\eta_{vb})]\},
\]

followed by configured q bounds. The default learns probability but not
variance. Coefficients have zero-centered isotropic Gaussian priors controlled
by `sigma_eta_pi` and `sigma_eta_vb`.

Each coefficient vector is proposed jointly by random-walk Metropolis, using
separate proposal SDs and an update interval (default every ten iterations).
The probability target uses the current Bernoulli marker states. The variance
target uses the active-effect Gaussian likelihood with current b, d, and
\(v_b\). When probability/q bounds bind, the clipped transform is part of the
implemented target.

The native update sequence is marker sweep, optional LD swap, annotation MH,
\(v_b\), \(v_e\), and global pi update. Acceptance counts/rates are reported;
coefficient traces are retained in extended modes.

#### Support and gaps

The route is ST CSR BayesC only. It has no MT, BED, eigen, BayesR mixture, or
probit-stick integration. It lacks LV-quality induced q/P stability summaries
and bulk/tail ESS/MCSE for coefficients.

#### Comparison with LV

The active-effect variance likelihood is the same underlying Gaussian
calculation as BayesC-LV where clipping is inactive. The models nevertheless
differ in preprocessing, bounds, coefficient prior parameterization, update
schedule/order, and RW-MH versus ESS. Exact draw equivalence is neither expected
nor an appropriate migration test.

#### Conclusion

LV should be canonical for learned Q. The learned-logistic probability half
should become a logistic probability provider. Freeze the old variance route,
qualify the replacement scientifically, then retire it rather than share two
implementations indefinitely.

### Probit-stick BayesRC/SBayesRC

#### Model

Each continuation step uses a probit regression based on \(A_j\alpha_k\); the
resulting stick probabilities define marker-specific prior component
probabilities. A small probability floor and row renormalization guard native
sampling. The null component remains the first component.

Non-intercept alpha coefficients have Gaussian priors with variance parameters
`sigmaSqAlpha`; one explicit intercept is allowed and has a separate proper
normal prior. `alpha_init` is initialization, not the prior mean. `updateAlpha`
and update-frequency controls govern latent-probit Gibbs updates;
`sigmaSqAlpha` has its own scaled-inverse-chi-square update.

#### Preprocessing

Current helpers align exact marker IDs, reject duplicate IDs and nonfinite
values, expand factors where relevant, standardize continuous columns, leave
binary columns uncentered under the historical profile, allow at most one
all-ones intercept, and reject nonintercept zero-variance columns. Exact
duplicate/rank-deficiency checks are not as uniform as LV, and ST legacy versus
BED/MT helpers are not completely consolidated.

#### ST, MT, and operator behavior

ST BayesRC/SBayesRC is implemented for BED, CSR, and eigen. MT BayesRC BED and
MT SBayesRC CSR/eigen are also implemented and tested. MT uses a common
component-probability alpha architecture and global conditional non-null trait
pattern probabilities. Neither sharing nor \(V_b\) is annotation dependent.

#### Output distinction

Some ST fields such as `annotation_pi` are plug-in probabilities derived from
posterior mean alpha and may be annotation-row diagnostics. MT raw output has a
clearer posterior-mean `prior_component_probabilities` field. In all routes,
posterior `component_probabilities` is a different object. Current ST/MT output
placement and naming are inconsistent.

#### Diagnostics and conclusion

Alpha/sigma traces, occupancy, eligible/continue counts, and convergence
diagnostics exist, with substantial development material on alpha mixing. Keep
probit-stick as the canonical component-probability provider. Do not expose
SBayesRC+LV until separate identifiability and mixing validation passes.

### Log variance: BayesC-LV and BayesR-LV

#### Preprocessing and identification

LV accepts an M by K annotation matrix aligned exactly to analysis marker order.
Binary columns are centered only. Continuous columns are centered and scaled to
sample SD one. It rejects nonfinite values, constants, all-one intercepts,
duplicates, and rank-deficient designs. Column names and transformation metadata
are retained; caller objects are not mutated.

With \(X\) centered,

\[
\eta_j=X_j\theta,\quad q_j=\exp(\eta_j),\quad
\overline\eta=0,\quad\operatorname{GM}(q)=1.
\]

#### Prior, sampler, and behavior

The default prior is \(\theta\sim N(0,0.7^2I)\); a different positive finite
scalar SD may be supplied. `theta_init` is separate. ESS samples theta, and an
empty active set triggers a direct prior draw. Log-scale guards prevent overflow
without silently clamping q.

BayesC-LV retains global pi. BayesR-LV retains global mixture probabilities and
uses one theta across non-null components. Scale-aware marker, global variance,
and LD-swap mathematics reuse the qualified ordinary engines.

#### Architecture and qualification

Separate scientific LV translation units attach policies to operator-neutral
BayesC/R engines. No-op ordinary policies consume no RNG; tests freeze ordinary
trajectories. LV is qualified for ST CSR and retained block eigen, including
theta-zero and fixed-q reductions, R-oracle learned-theta checks, CSR/eigen
concordance, convergence output, and numerical guards. BED and MT are absent.

#### Output

Outputs include theta, coefficient summaries (mean, SD, quantiles, R-hat,
bulk/tail ESS, MCSE), per-chain information/traces when requested, plug-in
`annotation_variance_ratio=exp(theta)`, posterior mean marker q, annotation
transform, and ESS evaluation/contraction/min-max-log-q diagnostics. Marker by
iteration q histories are not retained.

#### Conclusion

LV is the best existing template for a canonical learned variance provider and
for RNG-neutral policy design.

### MAF-S

MAF-S uses

\[
q_j(S)=h_j^{S+1}.
\]

It supports fixed S and sampled S on selected ST summary routes, including
scale-aware BayesC/R and SBayesRC configurations. Fixed MAF-derived scale is
used by established MT mixture routes; sampled MT S is explicitly unsupported.
BED and model-family coverage are not uniform despite similarly named public
arguments.

Sampled S uses bounded random-walk MH with explicit initialization, proposal,
and prior-range controls. Its historical q is not log-centered. Current
calibration includes the resolved q, but combining it with LV is unsupported.

Treat MAF-S as a variance provider. For new composition,

\[
\log q_j^{raw}=\log q_j^{external}+A_j\theta+(S+1)\log h_j.
\]

Adopt one declared identification convention for the new composed model, but
do not silently renormalize the validated historical MAF-S route. MAF columns
in A plus S require confounding warnings and separate simulation tests.

## Annotation preprocessing map

| Property | Older fixed/logistic helper | Group | Probit-stick BayesRC | LV |
|---|---|---|---|---|
| Marker-ID reordering | allowed when complete row names exist; otherwise row order trusted | group vector aligned separately | strict helper on BED/MT; ST legacy differs | exact aligned order required; no silent reorder |
| Binary columns | generally left uncentered | group labels, not A design | historically uncentered | centered, not SD-scaled |
| Continuous columns | centered/SD-scaled | not applicable | centered/SD-scaled | centered/SD-scaled |
| Intercept | optional helper behavior | implicit group baselines | one explicit intercept allowed and modeled | forbidden |
| Constant/all-one | not uniformly rejected | valid group labels checked | nonintercept zero variance rejected | rejected |
| Duplicate/rank deficient | not uniformly rejected | not applicable | not uniformly rank-checked | rejected |
| Nonfinite | rejected | invalid labels/parameters rejected | rejected | rejected |
| Transformation metadata | partial | group metadata | substantial | explicit and complete |
| Caller mutation | no intended mutation | no intended mutation | no intended mutation | tested absent |

### Recommended canonical design object

Use an immutable internal object with:

```text
X
marker_ids
original_names
processed_names
column_type
center
scale
intercept_policy
row_alignment
rank
duplicate/rank diagnostics
transform_profile
provenance
```

One constructor should own identity/alignment and metadata. Provider-specific
profiles own statistical transforms. Initially reproduce historical matrices
exactly; later harmonization is an intentional scientific migration.

## Identification, normalization, and composition

### Recommended new-provider contract

For fixed external q, LV correction, and future MAF-S composition, build raw
log q additively and subtract the unweighted marker mean:

\[
\log q_j=\log q_j^{raw}-M^{-1}\sum_l\log q_l^{raw}.
\]

Thus q has geometric mean one and redistributes variance relative to \(v_b\) or
\(V_b\). The aligned analysis-marker universe is part of the normalization
provenance. Binary variance ratios remain unchanged by the common shift.

### Existing conventions that must remain explicit

- LV already satisfies this contract through centered X.
- fixed-marker A-derived q is centered before clipping; clipping can alter GM.
- group optionally uses an arithmetic, marker-weighted normalization.
- historical MAF-S is uncentered.
- probability clipping and h2 calibration are different operations from q
  identification.

Do not force one convention into validated historical routes during a
mechanical reorganization. A new combined provider can adopt the coherent
contract while legacy adapters remain frozen until approved migration.

## Fixed, learned, and informative prior support

Current coverage is uneven:

| Quantity | Global | Fixed marker | Fixed A + coefficient | Learned coefficient | Informative coefficient prior | Fixed baseline + learned correction |
|---|---|---|---|---|---|---|
| BayesC inclusion | yes | ST CSR | ST CSR | ST CSR logistic | zero-centered isotropic only | no |
| BayesR component probability | yes | fixed alpha construction | probit-stick fixed alpha | probit-stick learned alpha | proper but limited alpha priors | no |
| Relative variance q | yes | ST CSR | ST CSR | LV ST CSR/eigen; older MH ST CSR | zero-mean scalar-SD LV only | no |
| MT sharing | yes, global | no | no | no | no | no |
| \(V_b\) | yes, global | not marker-specific | no | no annotation dependence | covariance priors/calibration are core-specific | no |

### Informative theta recommendation

Add independent named priors first:

\[
\theta_k\sim N(\mu_k,\sigma_k^2).
\]

- `theta_prior_mean`: named numeric vector; default zero.
- `theta_prior_sd`: positive finite scalar or named numeric vector; default 0.7.
- resolve names against processed annotation coefficient names;
- reject partial/extra/ambiguous names;
- record resolved values and preprocessing;
- keep `theta_init` scientifically separate from the prior mean;
- if `theta_init` is absent, resolve it to `theta_prior_mean`; otherwise use
  the explicit initialization;
- construct initial q and perform h2 calibration from that resolved
  initialization, while recording both values separately;
- add a later named `theta_prior_cov` only after positive-definite validation and
  a stable naming contract exist.

For shared-q MT-LV, one named vector is shared across traits. A matrix belongs
to the later trait-specific-q model, not to the first extension.

### External q recommendation

Expose the scientific object as `prior_variance_multiplier`, strictly positive
and finite, aligned by unique marker ID. Support external-only, learned-only,
and external+learned through one provider. Store both the supplied baseline and
normalization constant. Do not maintain both `fixed_vb_multiplier` and an LV
“offset” as independent implementations.

## Multi-trait design from the beginning

### Immediate model

Implement shared-q MT BayesR/SBayesR-LV while preserving the existing global
`joint_pi` state simplex:

\[
q_j=\exp(A_j\theta),\qquad
\operatorname{Cov}(b_{j,p,k})=\gamma_kq_jV_{b,p}.
\]

It preserves correlations implied by \(V_b\). The MT scale/\(V_b\) update and h2
calibration must divide/include q consistently. Component and pattern
probabilities remain marginals of the unchanged global joint simplex.
BayesRC/SBayesRC plus LV is excluded from this phase and belongs to the
separate simultaneous-P-and-Q research phase. MT BayesC-LV may remain a later
extension. ST and MT should share annotation design, log-q composition, theta
prior/ESS, diagnostics, and provider interfaces; they should not share effect
updates that have different sufficient statistics.

### Future-capable boundaries

- trait-specific q uses \(D_jV_bD_j\), requiring stronger scale constraints;
- annotation-dependent H requires a new pattern-probability provider;
- annotation-dependent \(V_b\) requires a positive-definite covariance-regression
  model.

These should remain possible at the type/interface level but absent from the
initial public API.

## Operator and RNG architecture

The LV extraction established a sound principle: operator-neutral sampler
engines consume narrow prior policies; ordinary policies are no-ops and consume
zero random numbers. Binding-neutral retained-block execution lets ordinary and
LV policies reuse one likelihood implementation.

This principle should extend to future providers and MT:

- probability/variance policies own prior state and their scheduled updates;
- BED/CSR/eigen operators own likelihood/residual mechanics;
- lifecycle hooks remain narrow (for example, the qualified post-vb theta
  update), not a general callback framework;
- disabled features add no RNG calls and negligible branching overhead;
- operator-specific code is unavoidable for residuals/data access, not prior
  mathematics.

Ordinary same-seed trajectories are permanent regression controls.

## Capability audit and drift

The detailed current/proposed table is in
`docs/dev/annotation_prior_architecture_matrix.md`.

Material findings:

- the executable capability resolver omits BayesC/R-LV entirely;
- it reports MT BayesRC/SBayesRC routes unsupported although public R/native
  routes and tests implement them;
- the older developer capability matrix also omits LV and contains stale public
  method claims;
- public signatures sometimes accept arguments for combinations that route
  validation later rejects;
- `stblr_block_eigen` exposes both singular `annotation` and plural
  `annotations` concepts;
- AGENTS' canonical public list omits the exported block-eigen and MT fitters.

Recommendation: one declarative capability registry, tested against public
dispatch, native registration, schema requirements, and generated docs. A
function argument is not capability evidence.

## Raw schema and formatter audit

### Intended path

```text
native named stblr_raw
  -> stblr_raw_schema validation
  -> canonical raw-to-fit formatter
  -> stable fit
  -> summaries and diagnostics
```

Legacy positional formatting is correctly not an active architecture and
should not return.

### Current strengths

- common ST and MT raw structures are named and versioned;
- canonical common fields and present-but-NULL diagnostics are formatted
  consistently;
- MT BayesRC validation checks its annotation component structures deeply;
- route tests cover named raw contracts and failure behavior.

### Drift and special cases

- LV native adapters decorate raw `annotation`/`diagnostics`, but the generic
  ST schema and formatter do not own the complete LV contract. R wrappers save,
  restore, and attach LV fields after canonical formatting.
- ST fixed/group/logistic metadata is standardized by route-specific
  post-format helpers rather than one typed prior namespace.
- ST BayesRC exposes several alpha/probability objects at top level; MT BayesRC
  places more of the same science under `model_parameters$annotations`.
- `annotation_pi`, `prior_component_probabilities`, and
  `component_probabilities` do not have uniform semantics/names across routes.
- fit consistency requirements enumerate several annotation backends but are
  not yet a complete LV capability/schema registry.
- LV-specific model identities need special finalizer handling because the
  generic model identity set predates them.

### Recommended raw namespaces

```text
raw$architecture$probability
raw$architecture$variance
raw$architecture$sharing
raw$architecture$covariance
raw$annotation_design
raw$diagnostics$prior_*
raw$model_spec$prior
```

Each namespace carries a `kind`, fixed/learned status, parameters, summaries,
traces when requested, and provenance. The schema validates dimensions and
semantics before formatting. Native code should return scientific quantities;
R should attach input transformation metadata without reconstructing sampler
state.

## Current and proposed output map

The long-term structure should separate posterior/fixed architecture results,
the supplied prior specification, provenance, and diagnostics:

```text
fit$architecture$probability
fit$architecture$variance
fit$architecture$sharing
fit$architecture$covariance
fit$model_spec$prior
fit$annotation_design
fit$diagnostics
```

`model_spec$prior` is what the user specified or the model assumed, including
hyperparameters, fixed inputs, normalization, and resolved initialization.
`architecture` is the scientific result: fixed quantities and posterior
inference about alpha, theta, induced probabilities, q, sharing, and covariance.
This avoids describing posterior inference about prior-governing parameters as
though it were itself the supplied prior.

| Current field/concept | Proposed location | Required semantic clarification |
|---|---|---|
| `alpha`, `alpha_summary` | `architecture$probability$coefficients` / `$summary` | continuation-probit coefficients |
| `sigmaSqAlpha` | `architecture$probability$coefficient_variance` | learned hyperparameter result, distinct from its specification |
| `eta_pi` | `architecture$probability$coefficients` | logistic inclusion coefficients |
| `eta_vb` | retire after migration; map historical output to `architecture$variance$coefficients` | historical clipped RW-MH model |
| `theta`, `theta_summary` | `architecture$variance$coefficients` / `$summary` | canonical log-variance coefficients |
| `annotation_variance_ratio` | `architecture$variance$coefficient_ratio` | plug-in `exp(E theta)`, not posterior mean q |
| `marker_prior_scale` | `architecture$variance$marker_multiplier_mean` | posterior mean `E[q_j]` for learned q |
| `fixed_vb_multiplier` | `architecture$variance$external_multiplier` | fixed result/provenance; supplied value also recorded in `model_spec$prior` |
| fixed marker pi | `architecture$probability$marker_probability` | fixed prior architecture, not PIP |
| `group_pi`, traces | `architecture$probability$groups` | learned group architecture probabilities |
| group variance multiplier | `architecture$variance$groups` | fixed/learned and normalization metadata |
| `pis`, `pi` | `architecture$probability$global` or core mixture summary | distinguish trace/mean and BayesC/R shapes |
| `annotation_pi` | replace with explicit `prior_component_probability_plugin` or posterior mean prior P | current route semantics differ |
| `prior_component_probabilities` | `architecture$probability$marker_component_mean` | prior architecture averaged over MCMC |
| `component_probabilities`, `comp_prob` | `posterior$component_probability` or stable existing posterior namespace | posterior allocation, not prior |
| `dm`, PIP | stable posterior marker field | posterior non-null probability |
| `maf_effect_s*` | `architecture$variance$maf_s` | fixed/sampled result; assumptions remain in `model_spec$prior` |
| annotation transform | `annotation_design` | preprocessing/provenance, not posterior result |
| prior/init controls in `input` | `model_spec$prior` and `$initialization` | keep prior and init separate |
| acceptance/ESS/q stability | `diagnostics$prior_probability` / `$prior_variance` | mechanism-specific diagnostics |

An additive transition can populate the new structure while retaining stable
fields. Because backward compatibility is not a long-term constraint, a later
intentional breaking release should remove redundant top-level fields rather
than maintain aliases forever.

## Public argument audit and proposed mapping

Current public interfaces mix singular/plural annotation names, mechanism
selectors, core models, and mechanism-specific parameters. `annotation_model`
is overloaded as both a scientific choice and a dispatch switch. ST/MT and
operator routes expose different alpha, eta, theta, group, fixed-marker, and
MAF-S spellings.

### Recommendation

Use:

```r
prior = prior_spec(
  probability = ...,
  variance = ...,
  sharing = ...,
  covariance = ...
)
```

The fitter chooses the operator; `model="bayesc"` or `model="bayesr"` chooses
the core family. Provider helpers accept `annotations` and scientifically named
controls. BayesRC and BayesR-LV may remain convenient aliases/display labels
that resolve to this specification.

| Current | Proposed | Decision |
|---|---|---|
| `annotation`, `annotations` | `annotations` | retain plural only |
| `annotation_model` | `prior` provider specification | retire selector |
| `method="sbayes*"` | `model="bayesc"` / `"bayesr"` | operator already implies summary route; retain S name as display alias |
| `alpha_*` | probability-provider `coefficient_*` with alpha shown in scientific output | consistent provider controls; preserve alpha notation in theory |
| `sigmaSqAlpha_*` | `coefficient_variance_*` or explicit alpha-prior object | remove ambiguous camel-case public names |
| `eta_pi_*` | logistic probability provider `coefficient_*` | eta is internal symbol |
| `eta_vb_*`, `learn_vb_annot` | LV `theta_*` then retire | one learned-Q engine |
| `sigma_eta_*` | provider `coefficient_prior_sd` | distinguish pi and q provider scopes |
| `fixed_pi_marker` | `marker_probability` | within fixed probability provider |
| `fixed_vb_multiplier` | `prior_variance_multiplier` | within fixed/LV variance provider |
| `group_*` | group provider arguments | preserve group prior hyperparameters with consistent names |
| `maf_effect_s_init` | `initialization` inside MAF-S provider | not a prior |
| `maf_effect_s_prior` | `prior_range` (until a proper density replaces it) | state exact bounded-prior meaning |
| `updateAlpha`, `updateTheta`, `updatePi`, `updateGroupVb` | provider `update=TRUE/FALSE` | scoped, consistent boolean |
| `*_init` | `*_init` or provider `initialization` | always starting state |
| `*_prior_mean`, `*_prior_sd` | retain these scientific suffixes | named and separate from init |

### Design-mockup use cases

```r
# Ordinary BayesR: no annotation preparation or RNG overhead
stblr_csr(stats, Glist, model = "bayesr")

# Binary, continuous, overlapping, or mixed A use the same provider
stblr_csr(
  stats, Glist, model = "bayesr",
  prior = prior_spec(
    variance = annotation_log_variance(annotations = A)
  )
)

# Informative coding variance ratio
stblr_csr(
  stats, Glist, model = "bayesr",
  prior = prior_spec(
    variance = annotation_log_variance(
      annotations = A,
      theta_prior_mean = c(coding = log(3)),
      theta_prior_sd = c(coding = 0.3)
    )
  )
)

# External q, optionally recalibrated by learned A theta
stblr_csr(
  stats, Glist, model = "bayesr",
  prior = prior_spec(
    variance = annotation_log_variance(
      annotations = A,
      prior_variance_multiplier = q_external
    )
  )
)

# SBayesRC as BayesR plus annotation-dependent component probability
stblr_csr(
  stats, Glist, model = "bayesr",
  prior = prior_spec(
    probability = annotation_probit_stick(annotations = A)
  )
)

# Proposed shared-q MT-LV
mtblr_csr(
  stats, Glist, model = "bayesr",
  prior = prior_spec(
    variance = annotation_log_variance(annotations = A,
                                       trait_scope = "shared")
  )
)
```

The exact constructor vocabulary should be prototyped in R before locking it,
but the decomposition should not be compromised to shorten one call.

The final fitter shape is deliberately unresolved. Before the intentional API
cleanup, compare retaining the six operator-specific fitters
(`stblr_csr`, `stblr_bed`, `stblr_block_eigen`, and their MT counterparts)
with adding thin `stblr(data, ...)` and `mtblr(data, ...)` wrappers whose typed
data object selects the operator. Operator-specific functions may remain
advanced public interfaces in either design. The internal operator/prior
architecture must support both choices.

## Model naming conclusions

- Preserve BayesC and BayesR as core scientific families.
- Treat SBayesC/SBayesR as summary-likelihood display names, not distinct prior
  implementations.
- Treat BayesRC/SBayesRC as useful aliases for BayesR plus a probit-stick
  component-probability provider.
- Treat BayesC-LV/SBayesC-LV and BayesR-LV/SBayesR-LV as useful aliases/display
  identities for core family plus LV variance provider.
- Avoid combinatorial native model enums for every P/Q/H/\(V_b\) combination.
- Do not expose arbitrary composition simply because the type system permits
  it; capability validation remains scientific.

## h2 and prior-calibration conclusions

Current calibration now resolves marker probabilities, component weights,
gamma, q, genotype scale, groups, and supported MAF-S when constructing initial
expected genetic variance. For ST the key weight is

\[
W=\sum_j\sum_kp^{init}_{jk}\gamma_kq_jv_j.
\]

For MT a trait-pair weight matrix additionally includes activity patterns. An
automatic diagonal \(V_b\) can target marginal trait variances; unrestricted \(V_b\) is
not identified by marginal h2 alone.

The interpretation must be explicit:

- with fixed architecture, h2 can describe the exact expectation under the
  declared fixed prior state;
- with learned alpha/theta/group/S, current calibration generally uses the
  initialized architecture and is an initialization approximation;
- for informative theta priors, missing `theta_init` resolves to
  `theta_prior_mean`, while an explicit `theta_init` is used as supplied;
- initial q and h2 calibration use that resolved initialization, and both the
  prior mean and initialization are recorded separately; no theta-prior
  integration is used in the first implementation;
- external q enters before calibration and after declared normalization;
- external q plus theta uses the combined initialized q;
- MT-LV must include q in every trait-pair calibration weight and in \(V_b\)-update
  sufficient statistics.

The existing developer calibration audit contains historical “before change”
tables alongside later corrected contracts. It should be clearly archived or
rewritten so old defects are not read as current behavior.

## Diagnostics audit and target

### Current

- learned logistic: MH acceptances and optional eta traces; limited canonical
  P/q stability summaries;
- group: optional group pi/q traces and update summaries;
- probit stick: alpha/sigma traces, occupancy, eligible/continue counts, and
  multichain diagnostics; documented mixing challenges remain;
- LV: coefficient R-hat, bulk/tail ESS, MCSE, per-chain summaries/traces,
  induced q stability, ESS likelihood evaluations/contractions, and min/max
  log q.

### Required layered interpretation

1. **Latent coefficients:** convergence and sampler efficiency of alpha, eta,
   theta, S, or group parameters.
2. **Induced prior architecture:** stability of marker probabilities, q,
   expected active counts, component occupancy, and MT pattern probabilities.
3. **Posterior marker inference:** stability of PIP, component allocation, and
   beta.
4. **Downstream inference:** genetic variance and prediction stability.

Weak coefficient identification with stable q can be scientifically acceptable;
stable coefficient means with unstable occupancy, PIP, or prediction are not.

## Documentation audit

### Canonical organization recommended

- `docs/methods/model_theory.qmd`: concise core ST/MT BLR and BayesC/BayesR
  theory;
- `docs/methods/summary_statistics.qmd`: cross-product and sufficient-statistic
  theory;
- `docs/methods/multitrait_overlap.qmd`: `SSY`, `Nmat`, sample-overlap,
  and residual-covariance assumptions;
- `docs/methods/sparse_ld_csr.qmd` and `block_eigen_operator.qmd`: operator
  mathematics, not prior theory;
- `docs/methods/annotation_priors.qmd`: canonical theory/model definitions,
  current capability, and clearly labelled proposals;
- `docs/methods/annotation_informed_models.qmd`: focused empirical-Bayes and
  posterior-to-prior derivations;
- `docs/methods/sbayesrc_annotations.qmd`: probit-stick details;
- MT BayesR and BayesRC methods pages: the current full-joint and factorized
  state architectures, respectively;
- `docs/notes/`: practical selection/workflows only;
- `docs/dev/`: architecture, schemas, audits, qualification, migration, and
  historical decisions;
- `local_reference/`: frozen local oracles and design material, not assumed to
  be rendered package documentation.

### Stale, duplicated, or ambiguous material found

- the old `annotation_priors.qmd` and `annotation_informed_models.qmd` overlapped
  probability/variance derivations while neither was a complete current LV/MT
  capability source;
- `annotation_informed_models.qmd` contains extensive conceptual pseudocode and
  research-stage direct MT extensions that could be mistaken for APIs without
  a stronger status banner;
- `blr_model_capability_matrix.md` and the executable resolver omit LV and
  disagree about MT BayesRC;
- annotation practical notes omit LV or imply broader composability than public
  dispatch permits;
- single-trait overview material does not consistently list LV;
- several SBayesRC developer files describe exploratory MCEM, mixing, or
  promotion variants; valuable history should remain, but an index/status label
  should distinguish production, qualified experiment, and abandoned proposal;
- the calibration audit mixes historical defects and current corrected rules;
- the v2 local audit brief has a v1 internal title/reference.

## Test audit

### Strong current coverage

- LV deterministic eta/q/log-likelihood fixtures, theta-zero/fixed-q
  reductions, empty-active prior behavior, 1-D ESS, numerical guards, R-oracle
  learned fits, CSR/eigen parity, block trajectory freezes, RNG neutrality, and
  output diagnostics;
- fixed/group/logistic raw fields, interface validation, h2 calibration,
  multichain behavior, LD-swap scale handling, and extended-diagnostic RNG
  tests;
- BayesRC exact annotation alignment, alpha updates and reductions,
  component/pattern probabilities, BED/CSR/eigen and ST/MT routes, retention and
  worker independence;
- MAF-S fixed/sampled formulas and supported route tests;
- named raw schemas, formatter stability, and fit consistency;
- prior calibration formulas for nonuniform P/Q and MT patterns.

### Tests required before future phases

- one common annotation-design fixture suite across every route/profile;
- capability registry versus public dispatch/native registration/schema/docs;
- informative theta nonzero-mean and named-vector/grid oracles;
- external q only and external+theta exact log-additive tests;
- combined q normalization and calibration tests;
- typed raw-schema validation for all P/Q/H/\(V_b\) providers;
- full current-to-future field mapping tests;
- shared-q MT-LV density, scale/\(V_b\)-update, calibration, ESS, reduction,
  operator-parity, recovery, and RNG tests;
- induced P/q cross-chain diagnostics for all learned providers;
- separate SBayesRC+LV identification/mixing simulation suite before enabling
  composition;
- performance tests ensuring annotation design is prepared once and disabled
  providers have negligible overhead.

## Explicit answers to the required architectural questions

1. **Is P/Q/sharing/covariance right?** Yes. Add core mixture ladder and data
   operator as orthogonal axes.
2. **Should `annotation_model` survive?** No, not long term. Replace it with a
   typed prior specification after an additive transition.
3. **Specify probability and variance independently?** Yes, subject to an
   explicit scientific capability registry.
4. **Should LV be canonical learned variance?** Yes.
5. **What happens to learned-logistic variance?** Freeze, compare, then retire;
   retain/reorganize its logistic probability provider.
6. **Should fixed-marker and external LV q share a provider?** Yes; they are the
   same scientific quantity.
7. **Should group remain distinct?** Preserve the hierarchy, express it as a
   convenience P/Q provider rather than a core model.
8. **Canonical A preprocessing?** One immutable alignment/metadata object with
   provider-declared transform/intercept profiles. LV's transform is canonical
   for log variance.
9. **How identify q?** New composable q uses mean log q zero over the aligned
   marker universe. Historical conventions migrate only explicitly.
10. **Informative theta interface?** Named `theta_prior_mean`, scalar or named
    `theta_prior_sd`; later named covariance. Init remains separate.
11. **External variance interface?** Positive aligned
    `prior_variance_multiplier`, not a generic `offset`.
12. **External q plus theta?** Add on log scale, then apply the declared combined
    normalization.
13. **MAF-S plus LV?** Compose on log scale eventually, warn for MAF-related A,
    and qualify a new normalization/calibration contract without changing the
    legacy route.
14. **First MT-LV?** MT BayesR/SBayesR with its global joint-state simplex and
    one shared q with covariance \(\gamma_kq_jV_{b,p}\). BayesRC/SBayesRC+LV is
    Phase 8; MT BayesC-LV is later.
15. **How share ST/MT implementation?** Share annotation design, q provider,
    theta prior/ESS, diagnostics, and lifecycle contract; keep ST/MT sufficient
    statistics and effect updates specialized.
16. **Future MT possibilities?** Trait-specific q, annotation-dependent H, and
    annotation-dependent \(V_b\) remain type-level possibilities and research models.
17. **Which model names survive?** BayesC/R core names; BayesRC and LV names as
    aliases/display summaries; S prefixes as likelihood-route display names.
18. **Which arguments survive?** `annotations`, `theta_init`, and explicit
    `*_prior_mean`/`*_prior_sd`; retire singular annotation, annotation_model,
    eta_vb, and fixed_vb naming after migration.
19. **Fit output?** Typed nested prior/result namespaces plus separate
    annotation design, model specification, and diagnostics.
20. **Redundant interfaces?** learned-logistic q, fixed-marker q versus future
    external q, singular/plural annotation, model-specific update names, and
    route-specific formatter decoration.
21. **Redundant/stale docs?** overlapping annotation theory pages, stale
    capability/calibration pages, and insufficiently labelled experimental
    SBayesRC material.
22. **Features to retire?** Retire redundant API/routes only after equivalence or
    intentional migration: learned-logistic variance, fixed_marker as a model
    selector, annotation_model, duplicate aliases, and special formatter paths.
23. **h2 rules?** Calibrate with the full resolved initial P/Q/H/gamma/genotype
    architecture; report whether exact fixed expectation or initialized learned
    approximation; include q in MT trait-pair weights.
24. **Raw schema changes?** Typed P/Q/H/\(V_b\) and annotation-design namespaces,
    validated before formatting; no post-formatter reconstruction.
25. **Tests before each phase?** Mathematical oracle, ordinary/fixed reductions,
    same-seed RNG trajectory, route/operator parity, raw/fit schema, calibration,
    diagnostics, performance, and full regression gates.
26. **Lowest-risk sequence?** Freeze; common design/capability/schema; canonical
    LV; informative theta; external q/composition; retire redundancy; shared-q
    MT-LV; operator completion; separately validate P+Q composition.

## Open scientific decisions requiring approval

1. Confirm unweighted marker mean log q as the normalization for **new composed
   providers**, while leaving historical MAF-S/group/fixed clipping behavior
   unchanged until explicit migration.
2. Decide whether old learned-logistic q should remain available under a
   “historical clipped MH” name or be removed once LV replacement is qualified.
3. Decide whether probability-provider coefficient names should remain
   scientific (`alpha`, `eta`) in public controls or use generic provider names
   while retaining alpha/theta in outputs/theory.
4. Decide whether group variance's arithmetic normalization is retained as a
   distinct hierarchical convention or migrated to geometric normalization in
   a new model version.
5. Decide whether BayesRC/SBayesRC remains a first-class public alias after the
   structured prior API is stable. Keeping it as a convenience alias is
   recommended.
6. Decide whether the final easy-to-use API remains operator-specific or gains
   thin `stblr()`/`mtblr()` wrappers over typed data objects. Do not remove the
   operator-specific advanced interfaces merely by adding wrappers.

## Recommended implementation sequence

The gated, file-level plan is in
`docs/dev/annotation_prior_architecture_implementation_plan.md`. In summary:

1. freeze current behavior and repair capability truth;
2. add the common annotation design and typed raw namespaces without changing
   matrices or samplers;
3. consolidate learned Q around LV;
4. add informative theta priors;
5. add external q and log-additive correction;
6. intentionally retire redundant APIs/outputs;
7. implement shared-q MT BayesR/SBayesR-LV with the existing global joint-state
   simplex unchanged;
8. complete BED/other operators only through clean adapter reuse;
9. validate SBayesRC+LV as a separate research model before exposure.

This sequence maximizes reuse of qualified provider/policy machinery, keeps
ordinary trajectories as controls, and separates mechanical reorganization
from changes to the statistical model.
