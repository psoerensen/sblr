# Annotation-prior architecture implementation and reorganization plan

## Purpose

This plan turns the repository audit into gated work. It is deliberately not a
promise to expose every composable prior. Each phase has a scientific target,
an implementation boundary, and stop/go criteria. Existing validated samplers
remain regression controls until an explicitly approved phase changes a model.

The target decomposition is:

\[
\text{data/operator} \;\perp\;
\{\text{core family}, P, Q, H, V_b\},
\]

where \(P\) is marker state/component probability, \(Q\) is relative marker
effect variance, \(H\) is MT trait-pattern sharing probability, and \(V_b\) is
cross-trait effect covariance. BayesR's fixed component ladder \(\gamma\)
belongs to the core family, not to the annotation provider.

Theory uses scalar \(v_b,v_g,v_e\) for ST variance quantities and matrix-valued
\(V_b,V_g,V_e\) for their MT covariance counterparts. This notation does not
rename package fields: code retains `vb`/`covb`, `vg`/`covg`, and `ve`/`cove`,
and future R APIs must avoid case-only names. Trait-specific scalar `v_le` and
`v_ld` summaries remain lowercase.

## Non-negotiable gates used in every phase

- Ordinary disabled-feature paths consume no additional random numbers and
  preserve same-seed trajectories.
- BED decoding, CSR alignment, retained-eigen residual policy, marker update
  order, OpenMP scheduling, and LD-swap mathematics change only in a separately
  authorized and qualified phase.
- R preprocessing happens once. Native code validates dimensions and finite
  state but does not silently transform annotations again.
- Named `stblr_raw` fields are validated before canonical formatting. No
  positional fallback is introduced.
- Prior component probability, posterior component allocation probability,
  PIP, beta, and prediction diagnostics are never conflated.
- A phase stops if its independent mathematical oracle, reduction gate,
  trajectory gate, or full regression suite fails.

## Phase 0 — Freeze the current scientific checkpoint

### Scientific objective

Turn the current repository behavior into an explicit baseline before moving
ownership or changing public abstractions.

### Mathematical target

No model change. Freeze global BayesC/R, fixed-marker, group,
learned-logistic, BayesRC/SBayesRC, LV, and MAF-S behavior exactly as currently
implemented, including bounds, normalization, update order, and RNG streams.

### Current implementation affected

Only fixtures, capability metadata, and audit documentation. Do not alter
sampler kernels.

### R and native files

- R: capability resolver, raw-schema tests, public-dispatch tests.
- Native: no edits expected; record registered entry points and raw namespaces.

### Raw schema and formatter

Freeze representative named raw objects and formatted fits for every supported
route. Add explicit assertions for currently present annotation namespaces and
for present-but-`NULL` stable fields.

### Public API and output

No change.

### Tests

- same-seed ordinary C/R trajectories for BED, CSR, and eigen;
- fixed-marker/group/logistic trajectory fixtures;
- BayesRC alpha and allocation fixtures;
- LV zero-theta, fixed-theta, ESS, CSR/eigen, and RNG-neutrality fixtures;
- MAF-S fixed and sampled route fixtures;
- raw-schema and formatter snapshots;
- executable capability registry versus exported dispatch and tests.

### Stop/go

Go only after every claimed current cell is backed by a route test or is
downgraded to “not fully qualified.” Stop if documentation claims cannot be
reconciled with executable dispatch.

## Phase 1 — Common terminology, annotation design, and capability registry

### Scientific objective

Create one lossless annotation preparation contract without changing the
statistical transforms used by existing routes.

### Mathematical target

Introduce an immutable internal `AnnotationDesign` concept containing aligned
raw values, processed \(X\), marker IDs, column names/types, centers, scales,
intercept policy, rank information, and provenance. The provider—not a generic
helper—declares whether an intercept is part of its model.

### Current implementation affected

Consolidate duplicated logic in `annotation-helpers.R`,
`stblr-logvar-annotations.R`, and BayesRC/MT alignment helpers behind explicit
profiles. Initially, each profile must reproduce its historical transformation
exactly. Do not silently make learned-logistic use the LV transform in this
phase.

### R files/functions

- add an internal annotation-design constructor and validators;
- make public wrappers consume a prepared object;
- replace duplicated capability tables with one declarative registry;
- test the registry against public model resolution.

### Native files/functions

No statistical changes. Native boundaries receive the same numeric matrices as
before plus typed metadata only if necessary.

### Raw schema and formatter

Add a common raw annotation-design/provenance namespace. Format it as
`fit$annotation_design`; retain existing top-level fields during the additive
migration.

### Public API and output

No breaking change yet. Internally standardize on `annotations`; deprecate the
singular `annotation` alias only after parity tests pass.

### Documentation

Make `docs/methods/annotation_priors.qmd` canonical. Recast
`annotation_informed_models.qmd` as a focused empirical-Bayes construction
note. Generate capability documentation from the registry where practical.

### Tests

- row alignment/reordering and strict failure fixtures;
- duplicate row IDs, duplicate columns, rank deficiency, constants, nonfinite
  values, factor expansion, and caller-object immutability;
- exact matrix parity for every historical profile;
- cross-operator prepared-object parity;
- capability registry versus wrapper/native/schema/test inventory.

### RNG/performance gates

Prepared numeric matrices passed to native code must be byte-identical to the
historical path. Measure preparation cost on large \(M\); construction must be
once per fit, not once per chain or iteration.

### Stop/go

Stop if one common object cannot represent all required metadata without
changing a historical transform. The object may retain multiple named profiles;
the statistical harmonization happens later and explicitly.

## Phase 2 — Canonical effect-variance provider and LV consolidation

### Scientific objective

Make LV the single learned annotation-variance engine while retaining an
independent probability provider.

### Mathematical target

The canonical learned variance model remains

\[
\log q_j=A_j\theta,\qquad
\theta\sim N(0,0.7^2I),
\]

with LV preprocessing, geometric-mean identification, ESS, and no clipping.
The existing learned-logistic probability half remains
\(\operatorname{logit}(\pi_j)=\operatorname{logit}(\pi)+A_j\eta_\pi\)
until separately redesigned.

### Current implementation affected

- LV policy/provider and `st_logvar_annotation_prior` become canonical for
  learned \(Q\).
- Freeze the old learned-logistic variance route as a historical MH/clipped
  model; do not call it equivalent to LV.
- After scientific side-by-side qualification, retire `learn_vb_annot` rather
  than maintain two learned-\(Q\) implementations.

### R files/functions

Model specifications resolve probability and variance independently. Split
learned-logistic argument preparation and formatting into probability and
variance namespaces.

### Native files/functions

Reuse the qualified LV policy with BayesC/R engines and all qualified summary
operators. The probability policy is independent. No generic callback system;
retain narrow lifecycle hooks and no-op policies.

### Raw schema and formatter

Add typed architecture-result namespaces for probability and variance, plus a
separate `model_spec$prior` namespace for the supplied assumptions. Move LV
attachment out of post-formatter decoration and validate theta, q, traces, and
ESS diagnostics in the raw schema.

### Public API and output

Introduce an experimental structured `prior` input specification while keeping
current wrappers as adapters. New scientific output is additive under
`fit$architecture`; the resolved input specification is recorded separately
under `fit$model_spec$prior`.

### Tests

- all existing LV oracles and trajectory gates;
- historical learned-logistic variance fixtures until retirement;
- independent probability-only and variance-only reductions;
- no-RNG ordinary and disabled-provider tests;
- schema failure tests for malformed provider output.

### Performance gates

No material cost for global \(Q=1\); LV cost scales with active markers times
annotation dimension plus ESS evaluations. Record ESS contraction/evaluation
distributions.

### Stop/go

Retire the old variance route only after users can reproduce the intended
scientific use through LV and after the change is explicitly accepted as a
model migration, not a refactor.

## Phase 3 — Informative theta priors

### Scientific objective

Allow external scientific information about annotation variance ratios to be
updated by current data.

### Mathematical target

First production form:

\[
\theta_k\sim N(\mu_k,\sigma_k^2),
\]

with named, independently distributed coefficients. Later covariance support
may use \(N(\mu,\Sigma)\), but it is not required here.

### Current implementation affected

Generalize the LV ESS prior ellipse from zero/isotropic to mean-shifted,
diagonal scale. The likelihood and q identification do not change.

### R files/functions

Add `theta_prior_mean` and `theta_prior_sd` to the variance-provider
specification. Accept a scalar SD or named vector; accept a named mean vector.
Resolve names after annotation preprocessing. Keep `theta_init` scientifically
separate. If it is absent, resolve it to `theta_prior_mean`; if supplied, use
the explicit value. Construct initial q from the resolved `theta_init` and use
that q for h2 calibration.

### Native files/functions

The ESS kernel draws an auxiliary vector from the declared Gaussian prior and
operates on \(\theta-\mu\). Empty-active behavior draws from that same prior.

### Raw schema and formatter

Store resolved prior mean/SD, coefficient names, and initialization in model
specification/provenance as distinct fields; posterior summaries remain
scientific architecture results.

### Public API and output

Example target: coding annotation prior
`theta_prior_mean = c(coding = log(3))`,
`theta_prior_sd = c(coding = 0.3)`. A scalar SD is recycled only after names and
dimensions are resolved. With no explicit `theta_init`, the resolved
initialization is therefore `c(coding = log(3))`, not zero.

### Documentation, migration, and performance

Document coefficient-name resolution, the prior variance-ratio interpretation,
and the difference between prior mean and initialization. Existing zero-mean
scalar-SD calls must pass through an identity adapter. ESS cost should remain
the same order as current LV; benchmark diagonal-scale draws and require no
material regression for the default isotropic case.

### Tests

- analytic 1-D and deterministic-grid 2-D posterior fixtures;
- nonzero-mean empty-active prior draws;
- named reordering and missing/extra name failures;
- scalar/vector equivalence;
- init/prior separation;
- missing-init resolution to the prior mean and explicit-init override;
- initial-q and h2-calibration agreement with the resolved initialization;
- separate recording of prior mean and resolved initialization;
- old zero-mean 0.7 trajectory identity.

### Stop/go

Stop on any ambiguity between raw and processed annotation names. Do not add a
full covariance parameter until a stable named-matrix interface and Cholesky
validation are designed.

## Phase 4 — External marker variance and learned recalibration

### Scientific objective

Unify fixed marker variance information and LV correction.

### Mathematical target

\[
\log q_j^{\rm raw}=\log q_j^{\rm external}+A_j\theta,
\qquad q_j^{\rm external}>0.
\]

For the new provider, subtract the unweighted mean log scale over the aligned
analysis-marker universe so \(\operatorname{GM}(q)=1\). This preserves relative
variance and leaves the global scale to \(v_b\) or \(V_b\).

### Current implementation affected

Extract the fixed `vb_multiplier` implementation behind the same variance
provider used by LV. Existing `fixed_marker` remains a regression adapter until
the new route is qualified.

### R files/functions

Public scientific name: `prior_variance_multiplier`, not `offset`. Resolve by
marker ID, require unique complete alignment and positive finite values, and
record pre/post-normalization metadata.

### Native files/functions

Provider supplies current q and the sufficient denominators required by
marker, \(v_b\)/\(V_b\), and LD-swap updates. No second external-q mechanism.

### Raw schema and formatter

Record external baseline, normalization constant, learned correction, combined
posterior mean q, and whether q is fixed or learned. Do not retain marker by
iteration q histories by default.

The formatter maps these fields into one variance result and does not
reconstruct the combined q from rounded summaries.

### Tests

- external-only equivalence to qualified fixed-marker scale;
- external vector versus fixed \(A\)+coefficient construction;
- external+theta log-additivity;
- normalization and calibration;
- LD swap and operator parity;
- invalid zero, negative, nonfinite, duplicate, missing, or reordered inputs;
- theta-zero and no-provider ordinary reductions.

### Documentation, migration, and performance

Document the marker universe used for normalization and show external-only and
external+learned examples. Existing `fixed_vb_multiplier` calls become exact
adapters before that name is deprecated. Construction is O(M) and performed
once; each learned update adds only the existing O(MK) LV calculation. Record
memory for the fixed baseline and posterior mean q only.

### Stop/go

Stop if a route cannot thread q consistently through marker density, global
scale update, and LD-swap ratio. Do not expose partial variance scaling.

## Pre-Phase-5 public fitter design checkpoint — decision deferred

Before intentional breaking cleanup, compare two compatible public shapes:

1. retain the six operator-specific fitters:
   `stblr_csr`, `stblr_bed`, `stblr_block_eigen`, `mtblr_csr`, `mtblr_bed`, and
   `mtblr_block_eigen`;
2. add thin `stblr(data, model=..., prior=...)` and
   `mtblr(data, model=..., prior=...)` wrappers, where a typed data object
   selects BED, CSR, or retained block eigen.

Do not decide or implement this during the architecture phases. The
operator-specific functions may remain advanced public interfaces even if
high-level wrappers are later added. Annotation design, prior providers, raw
schemas, and formatters must remain operator-neutral enough to support either
choice. Record the decision before Phase 5 removes or renames public arguments.

## Phase 5 — Retire and reorganize redundant historical routes

### Scientific objective

Remove redundant interfaces only after the common providers reproduce or
deliberately supersede their scientific roles.

### Target changes

- Replace `annotation_model` with an explicit prior specification.
- Retire `learn_vb_annot`, `eta_vb_*`, and duplicated fixed variance arguments.
- Make `group` a convenience constructor for group probability and/or variance
  providers while retaining its hierarchical priors.
- Retire singular `annotation`, `fixed_vb_multiplier`, and ambiguous
  `sigma_eta_*` names.
- Keep BayesRC and LV as convenient display aliases, not separate combinatorial
  internal core families.

### Current implementation, R, and native ownership

Replace route-specific R argument parsing and post-format decorators with
provider constructors and one formatter. Ordinary C/R engines, BayesRC kernels,
and LV kernels remain unchanged; delete native duplicate paths only where a
qualified provider has already assumed ownership. Generated bindings change
only if a genuinely retired native entry point is removed in the approved
breaking release.

### Raw schema, formatter, and output

Make nested `fit$architecture` results canonical and record the supplied prior
under `fit$model_spec$prior`. Remove redundant top-level fields only in one
intentional breaking release. Provide a machine-readable field mapping in the
release notes; no silent aliases.

### Public API and documentation

Adopt `annotations`, `model`, and structured `prior`; remove
`annotation_model`, singular `annotation`, and redundant model-specific q
controls. Update all methods, notes, examples, capability output, and reference
indexes together. Preserve old scientific names in migration documentation,
not as indefinite runtime aliases.

### Tests

- golden conversions from every old specification to its new provider;
- explicit errors for removed ambiguous arguments;
- complete raw and formatted schema tests;
- documentation examples for all common analyses.

### RNG, performance, and migration gates

Calls translated from old syntax must reach the same provider state before
native execution and preserve same-seed trajectories. The structured API may
not add per-iteration R dispatch. Publish an explicit removed/renamed field and
argument table, with a one-release conversion helper if useful to the sole
developer; do not silently reinterpret serialized fits.

### Stop/go

Each retirement requires an explicit decision that the old scientific model is
duplicated or no longer desired. A route is not removed merely because its API
is awkward.

## Phase 6 — Shared-q MT BayesR/SBayesR-LV

### Scientific objective

Add shared-q annotation-informed effect magnitude to MT BayesR/SBayesR while
preserving its existing global `joint_pi` over complete MT states, without
changing cross-trait correlation or state probabilities.

### Mathematical target

The only implementation target in this phase is MT BayesR/SBayesR-LV with
the ordinary model's full joint-state simplex unchanged:

\[
q_j=\exp(A_j\theta),\qquad
\operatorname{Cov}(\beta_{j,p,k})=\gamma_k q_j V_{b,p}.
\]

One theta and one q are shared across traits. Complete joint-state
probabilities and \(V_b\) remain global and annotation independent. Component
and pattern probabilities remain marginals of `joint_pi`; Phase 6 does not
refactor ordinary MT BayesR into independent \(P\times H\) factors.
\(V_{b,p}\) is the active-trait submatrix of \(V_b\); \(\gamma_k\) and q are
relative component and marker multipliers, while \(V_b\) carries the global
active-effect covariance scale.

### Current implementation affected

Only MT BayesR/SBayesR is affected. BayesRC/SBayesRC is not part of Phase 6,
even in a restricted variant; all simultaneous learned annotation-dependent P
plus Q work belongs to Phase 8. MT BayesC-LV may be considered in a later
extension. The MT global-scale update must use q-aware sufficient statistics,
and automatic \(V_b\) calibration must use q in its weight matrix.

### R files/functions

All MT wrappers consume the common annotation design and variance provider.
Reject trait-specific theta specifications in this phase.

The public structured prior accepts `trait_scope="shared"`; no new combinatorial
MT method string is required. Existing MT calls remain ordinary no-op paths.

### Native files/functions

Share annotation math and provider lifecycle with ST, but use MT-specific
effect and \(V_b\) update adapters. For an active pattern \(p\), q multiplies the
whole \(V_{b,p}\); it must not scale traits separately or change the correlation
structure represented by \(V_b\).

### Raw schema and formatter

The same `architecture$variance` result namespace is used by ST and MT. Store
shared theta once, not once per trait. Report q stability alongside \(V_b\) and
pattern diagnostics; record the global probability specification separately in
`model_spec$prior`.

### Documentation and migration

Update MT BayesR theory, the capability registry, calibration documentation,
and worked examples. Label BayesRC/SBayesRC+LV, MT BayesC-LV, and trait-specific
q as outside this phase.
There is no old MT-LV route to migrate, but ordinary MT serialized output must
remain schema-compatible during the additive phase.

### Tests

- deterministic covariance-density and scale-update fixtures;
- theta-zero exact ordinary MT BayesR reduction;
- fixed-theta equivalence to fixed marker q;
- isolated-theta grid/ESS oracle;
- BED/CSR/eigen exact or MC-aware parity;
- \(V_b\) and q recovery, R-hat/ESS, pattern occupancy, beta/PIP/prediction stability;
- h2 calibration with nonuniform q;
- ordinary MT trajectory/RNG neutrality.

### Performance gates

No marker-by-trait-by-iteration q storage. Shared q is computed once per theta
update. Benchmark added annotation cost separately from MT state enumeration.

### Stop/go

Stop unless global-joint-state MT BayesR/SBayesR-LV and shared-q covariance
behavior qualify independently. Do not add BayesRC/SBayesRC code to this phase.

## Phase 7 — Operator completion, including BED LV

### Scientific objective

Attach the already qualified variance policy to remaining operators without
cloning scientific samplers.

### Mathematical target

No change from the qualified ST or MT LV model.

### Current implementation affected

Binding-neutral BED execution adapters and operator construction only. Genotype
decoding and residual policy remain untouched.

### R, native, schema, formatter, API, and output

R dispatch attaches the existing variance provider to the BED execution input;
native code reuses the qualified engine and BED operator. Extend the same raw
variance namespace and canonical formatter—do not add a BED-specific LV output
path. The structured public prior call is operator-independent and the fit
exposes the same theta/q/diagnostic contract as CSR/eigen.

### Tests and gates

- ordinary BED same-seed trajectories before/after extraction;
- theta-zero and fixed-q equivalence;
- individual-level R oracle on tiny fixtures;
- BED/summary agreement where the data construction makes that meaningful;
- memory and OpenMP safety checks.

### Documentation, performance, and migration

Update the generated capability table and BED method page only after
qualification. Benchmark decode/residual time separately from theta updates;
ordinary BED overhead must be negligible. No historical route is removed in
this phase.

### Stop/go

Stop if reuse requires redesigning BED decoding or broad native refactoring.
BED support is not justification for a cloned sampler.

## Phase 8 — Separately validate probability plus variance annotation models

### Scientific objective

Determine whether simultaneous annotation-dependent \(P_j\) and \(Q_j\) is
identifiable and computationally useful—not merely implementable.

### Candidate mathematical target

\[
P_j=\operatorname{probit\mbox{-}stick}(A_j\alpha),\qquad
q_j=\exp(A_j\theta).
\]

### Current implementation, R, and native ownership

Compose already-qualified probability and variance providers in an experimental
internal capability only. R requires explicit separate specifications for P and
Q; native engines provide the two narrow update hooks without merging alpha and
theta state. No public alias is added during the research phase.

### Raw schema, formatter, output, and documentation

Use the existing typed `architecture$probability` and
`architecture$variance` result namespaces, plus joint diagnostics for
alpha/theta dependence. The supplied composition remains in
`model_spec$prior`. Experimental outputs must state that the model is
unqualified. Theory documentation may derive the candidate model; usage
documentation must not advertise it before promotion.

### Required research before public exposure

- prior predictive identification under correlated and overlapping A;
- signal/null calibration with separate and shared annotations;
- alpha/theta posterior correlation and occupancy mixing;
- induced P and q stability across chains;
- PIP, beta, and prediction robustness;
- sensitivity to informative priors and component ladders;
- clear conditions or warnings for weak identification.

### RNG, performance, and migration gates

Disabling either provider reproduces its independently qualified model exactly.
Measure the combined update cost and occupancy degradation; do not optimize
away diagnostics. There is no public migration until the model is promoted.

### Stop/go

The architecture may permit composition internally, but public dispatch remains
disabled until this phase passes a pre-registered simulation suite. Existing
SBayesRC alpha/occupancy mixing concerns make this a hard scientific gate.

## Later research, not implementation backlog

### Trait-specific variance

\[
q_{jt}=\exp(A_j\theta_t),\qquad
D_j=\operatorname{diag}(\sqrt{q_{j1}},\ldots,\sqrt{q_{jT}}),\qquad
\operatorname{Cov}(\beta_j)=\gamma_k D_jV_bD_j.
\]

This requires new identification constraints between \(D_j\) and \(V_b\), a
trait-structured coefficient prior, and stronger calibration.

### Annotation-dependent sharing

Treat annotations in pattern probabilities as a new sharing provider. It needs
its own sparsity, label, and occupancy theory and must not be implied by current
MT BayesRC support.

### Annotation-dependent covariance

Changing \(V_b\) by marker annotation is a substantially different covariance
regression model. It remains research-stage and unsupported by design until
positive-definite parameterization and identification are established.

## Intended final cleanup order

1. Freeze tests and correct capability truth.
2. Add common design and typed raw namespaces without model changes.
3. Consolidate learned Q around LV.
4. Add informative theta priors.
5. Add normalized external q and composition.
6. Remove redundant arguments/routes in one intentional API break.
7. Implement and qualify shared-q MT BayesR/SBayesR-LV with its existing
   global joint-state simplex unchanged.
8. Complete operators where adapter reuse is clean.
9. Consider simultaneous probability+variance annotation only after dedicated
   scientific validation.

This order separates mechanical ownership changes from statistical changes,
keeps ordinary trajectories as controls, and prevents the public API from
promising combinations before their posterior geometry is understood.
