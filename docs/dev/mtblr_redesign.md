# MTBLR Statistical and Software Redesign

## Status

**Document status:** Design draft 0.3; covariance, Cheng-pattern, regional
reference, and heterogeneous-provider phases verified by standalone evidence

**Package:** `sblr`

**Scope:** Major redesign of the multi-trait Bayesian linear regression framework

**Research-design checkpoint baseline:** [`psoerensen/sblr@52d5021`](https://github.com/psoerensen/sblr/commit/52d50219ead086245dc38cb5e5958a1403c78d1d)

**Earlier Phase 3A package checkpoint:** [`psoerensen/sblr@e11abdf9`](https://github.com/psoerensen/sblr/commit/e11abdf9f2f49b03169449c680ae033df4992acf)

**External design reference:** [`reworkhow/JWAS.jl@038c6c8`](https://github.com/reworkhow/JWAS.jl/commit/038c6c8a497ec6e17e957657a5c13ae698377132)

**Intended repository location:** `docs/dev/mtblr_redesign.md`

This document defines the proposed statistical and software architecture for a
major overhaul of the MTBLR framework in `sblr`. It is a design document, not a
description of current package functionality. Current behavior remains defined
by executable code, public dispatch, maintained tests, developer contracts, and
the canonical Methods pages at the pinned package baseline.

The document deliberately separates:

- confirmed deficiencies in the current implementation;
- design principles that should govern the replacement;
- recommended phase-1 decisions;
- alternatives requiring an explicit scientific decision;
- later extensions that should not enlarge the first implementation phase.

No production implementation should begin until the remaining phase-1 decisions
in the decision register have been approved. The Cheng state geometry,
inverse-Wishart reference parameterization, null-collapsed conditional
completion and joint state transition are resolved by
`mtblr_covariance_design.md`; this does not promote the current hybrid source
transition.

## 1. Motivation

The current MTBLR framework already contains valuable infrastructure:

- packed-BED individual-level analysis;
- CSR and retained block-eigen summary-statistic operators;
- BayesC, BayesR, and BayesRC marker-state models;
- strict marker and allele alignment;
- multi-chain execution with chain-local random-number generators;
- convergence and memory diagnostics;
- named raw backend output and canonical formatted-fit processing.

However, the scientific MCMC layer grew across several implementation paths and
now contains duplicated update schedules, ambiguous covariance semantics, and
legacy covariance routines. The source audit identified a confirmed mismatch
between the marker-effect covariance used in marker transitions and the
covariance subsequently traced and returned. Further inspection shows that the
problem is structural rather than a one-line reporting error.

The redesign therefore has two objectives:

1. define one coherent family of multi-trait statistical models; and
2. make every supported data representation evaluate those models through a
   shared, testable sampler architecture.

## 2. Authority and evidence

The redesign follows the repository authority order in `AGENTS.md` and
`docs/dev/README.md`:

1. current executable code, public dispatch, tests, and schemas;
2. maintained developer contracts under `docs/dev/`;
3. approved architecture and implementation plans;
4. canonical theory under `docs/methods/`;
5. qualification records;
6. historical and research material;
7. practical notes under `docs/notes/`.

JWAS is used as an external implementation and design reference. It does not
override the `sblr` Methods theory. Where JWAS and `sblr` use different model
definitions, the distinction must be explicit.

## 3. Goals

The target MTBLR framework should:

1. provide one authoritative statistical definition for each model;
2. maintain exactly one authoritative sampled covariance state for each
   declared statistical covariance unit at each MCMC iteration: one $V_b$ in
   the global model, or one persistent $V_{b,r}$ per region in a regional
   model;
3. update global and regional statistical parameters according to their model
   definitions, independently of marker-set, provider, LD-block, or eigenblock
   partitioning;
4. distinguish latent base effects, realized marker effects, marker states,
   covariance parameters, and derived covariance quantities;
5. separate individual-level, no-overlap summary, and overlap-aware summary
   likelihoods;
6. share marker-state and covariance logic across BED, CSR, and block-eigen
   representations where their likelihood contracts permit it;
7. support exact small-state joint transitions and scalable coordinate or
   structured transitions;
8. retain deterministic, chain-local RNG behavior and static chain-level
   parallelism;
9. return a validated, named raw object with unambiguous posterior summaries,
   final states, traces, and identification metadata;
10. allow intentional API and schema improvement without preserving obsolete
    behavior solely for compatibility;
11. be testable against analytical conditionals, exact tiny-state oracles,
    operator reductions, and known-truth simulations.

## 4. Non-goals for the first implementation

The first implementation should not attempt to deliver all planned MT methods.
In particular, phase 1 should not include:

- arbitrary summary-statistic sample overlap;
- annotation-dependent $V_b$;
- annotation-dependent trait-sharing architecture $H$;
- learned MT MAF-$S$;
- MT log-variance architecture;
- automatic selection between multiple covariance-prior families;
- large-trait approximate state-space algorithms;
- single-step pedigree integration;
- missing phenotypes;
- categorical traits;
- every current MT model/operator combination.

The architecture must leave clear extension points for these features without
pretending that they are phase-1 functionality.

## 5. Current framework assessment

### 5.1 Components worth preserving

The following current components are conceptually sound and should be reused or
adapted:

| Component | Design disposition |
|---|---|
| Marker identity, order, allele, and genotype-scale validation | Keep |
| Summary-statistic and LD provenance | Keep |
| Packed-BED decoding and phenotype preparation | Keep |
| CSR sparse-LD view | Keep |
| Retained block-eigen view | Keep |
| Representation-neutral residual-score operations | Keep and formalize |
| One complete joint MT model per logical chain | Keep |
| Chain-local `std::mt19937` state | Keep |
| Static OpenMP dispatch across chains | Keep |
| BayesR joint-state descriptor | Keep as one state-policy implementation |
| BayesRC annotation preprocessing and alignment | Keep initially |
| Convergence calculations | Keep, with revised trace inputs |
| Memory estimation and hard trace guards | Keep |
| Named raw backend objects and validation | Keep, with schema migration |

### 5.2 Components requiring redesign

| Component | Reason |
|---|---|
| Global $V_b$ transition | Current transition and returned value are not one authoritative state |
| BED and summary sampler loops | Scientific update order is duplicated |
| Covariance-prior parameterization | The meaning of prior scale inputs is ambiguous |
| Global update frequency | Current $V_b$ updates depend on marker-set partitioning |
| Retention schedule | Marker summaries and covariance summaries do not consistently use one retained-draw definition |
| Genomic covariance output | Some summary cross-trait values are descriptive bilinear forms rather than identified covariances |
| Raw schema | Posterior means, traces, final states, and derived values are insufficiently separated |
| Model capability reporting | Must be derived from current dispatch/contracts, not a stale standalone matrix |

### 5.3 Components to retire or quarantine

Subject to a final reference search, the redesign should retire:

- the internal positional `R/interface_mtblr.R::sblr()` route;
- heuristic `sampleB()` and `sampleBset()` covariance routines;
- superseded `sampleB_old()` and `sampleBset_old()` implementations;
- commented alternative MT covariance kernels in `src/mtblr.cpp`;
- positional raw-output assumptions;
- legacy fields whose accumulators have no supported statistical meaning;
- any capability claim derived only from `.blr_model_capability_matrix()`.

Removal should happen only after the new vertical slice is validated and active
references have been checked.

## 6. Statistical object model

Let $T$ be the number of traits, $M$ the number of markers, and $N$ the number
of individuals in a common-sample individual-level analysis.

For individual-level data, write

$$
Y=XB+E_Y,
$$

where:

- $Y\in\mathbb R^{N\times T}$ is the centered phenotype matrix;
- $X\in\mathbb R^{N\times M}$ is the centered or standardized genotype
  matrix;
- $B\in\mathbb R^{M\times T}$ contains realized marker effects;
- $E_Y\in\mathbb R^{N\times T}$ contains residuals.

For marker $j$, distinguish four objects:

1. $\boldsymbol{\delta}_j$: the discrete trait-activity or mixture state;
2. $\boldsymbol{\beta}_j\in\mathbb R^T$: a complete latent base-effect vector;
3. $D_j=\operatorname{diag}(\boldsymbol{\delta}_j)$: the activity mask;
4. $\boldsymbol{\alpha}_j=D_j\boldsymbol{\beta}_j$: the realised BayesC effect
   entering the phenotype model.

The rows of $B$ are $\boldsymbol{\alpha}_j^\top$. This distinction is central.
A null or trait-inactive component of $\boldsymbol{\alpha}_j$ is exactly zero,
while the corresponding component of $\boldsymbol{\beta}_j$ may exist as an
augmentation variable. Component or marker scaling in BayesR is represented by
a separate scale matrix, not folded into $D_j$. Augmentation variables must not
be reported as realised marker effects.

The shared base marker-effect covariance is

$$
V_b\in\mathbb R^{T\times T}.
$$

The observation-level residual covariance is

$$
V_e\in\mathbb R^{T\times T}.
$$

The genomic covariance $V_g$ is derived from realized effects and a declared
genotype covariance or cross-product. It is not the same object as $V_b$.

## 7. Marker-state policies

The engine should not hard-code one state geometry. A state policy should
provide:

- a state descriptor;
- prior state probabilities;
- the diagonal transformation $D_j$;
- marker-transition calculations;
- probability-parameter sufficient statistics;
- state summaries and names;
- memory and state-count guards.

### 7.1 MT BayesC

General trait-pattern MT-BayesC$\Pi$, as in Cheng et al. (2018), is an existing
intended and publicly selectable model geometry. For a binary activity pattern

$$
a_p\in\{0,1\}^T,
$$

define

$$
D_j=\operatorname{diag}(a_p)
$$

and

$$
\boldsymbol{\alpha}_j=D_j\boldsymbol{\beta}_j,
$$

where $\boldsymbol{\beta}_j\sim N_T(0,V_b)$ is the complete latent slab vector
and $\boldsymbol{\alpha}_j$ is the realised effect. The state space contains one
null pattern and a declared set of non-null trait patterns. For $T=2$, the
complete default state space is

$$
\{(0,0),(1,0),(0,1),(1,1)\}.
$$

The restricted shared-state model containing only $(0,0)$ and $(1,1)$ is a
validation special case, not the general phase-1 target.

### 7.2 Shared-component MT BayesR

This is the current `sblr` statistical model. A non-null state consists of a
trait pattern $p$ and one positive component $k$ shared across active traits.

With component multiplier $\gamma_k$ and marker multiplier $q_j$, the
standalone shared-component prototype stores
$\boldsymbol{\theta}_j\sim N_T(0,\gamma_kq_jV_b)$ and, for an active
all-traits component, sets $\boldsymbol{\alpha}_j=\boldsymbol{\theta}_j$.
Equivalently, define a base latent vector
$\boldsymbol{\beta}_j\sim N_T(0,V_b)$ so that

$$
D_j
=
\operatorname{diag}(a_p),
$$

so that

$$
\boldsymbol{\alpha}_j
=
\sqrt{\gamma_kq_j}\,D_j
\boldsymbol{\beta}_j.
$$

The active-subspace covariance is

$$
\operatorname{Var}(\boldsymbol{\alpha}_{j,p}\mid V_b,z_j)
=
\gamma_kq_jV_{b,p}.
$$

Under the prototype's scaled-effect storage, the covariance sufficient
statistic divides
$\boldsymbol{\theta}_j\boldsymbol{\theta}_j^\top$ by $\gamma_kq_j$ exactly
once. An implementation storing completed base vectors instead uses
$\boldsymbol{\beta}_j\boldsymbol{\beta}_j^\top$ directly and must not remove
the scale again. The $q_j=1$ case is only a reduction of this contract.

There is one null state and $LK$ non-null states when $L$ non-null activity
patterns and $K$ positive components are supplied.

### 7.3 Trait-specific-component MT BayesR

JWAS motivates a more general state policy in which every trait has its own
component class:

$$
c_j=(c_{j1},\ldots,c_{jT}),
\qquad
c_{jt}\in\{0,1,\ldots,K\}.
$$

Define the activity mask and scale matrix separately:

$$
D_j
=
\operatorname{diag}
\left(
\mathbf 1(c_{j1}>0),
\ldots,
\mathbf 1(c_{jT}>0)
\right),
\qquad
G_j
=
\operatorname{diag}
\left(
\sqrt{q_j\gamma_{c_{j1}}},
\ldots,
\sqrt{q_j\gamma_{c_{jT}}}
\right),
$$

where $\gamma_0=0$. Then

$$
\boldsymbol{\alpha}_j=D_jG_j\boldsymbol{\beta}_j,
\qquad
\operatorname{Var}(\boldsymbol{\alpha}_j\mid V_b,c_j)
=
D_jG_jV_bG_jD_j.
$$

Because $\gamma_0=0$, $G_j$ already has zero entries for inactive traits, so
$D_jG_j=G_j$. The mask is retained in the notation to distinguish activity
from scale and to support implementations in which null-state masking is
represented explicitly.

This permits a marker to occupy a large-effect class for one trait and a
smaller-effect class for another. The unrestricted state count is

$$
(K+1)^T,
$$

so full joint enumeration is suitable only for small $T$ and $K$. Coordinate or
structured transitions are required for larger problems.

This policy is recommended as a planned extension, not as the first vertical
slice. Its marker transition, conditional completion, and covariance
sufficient statistic require a separate derivation because $G_j$ is generally
not a scalar multiple of the identity. The validated shared-component BayesR
derivation must not be applied to this trait-specific scale geometry without
that derivation.

### 7.4 Transition strategies

The engine should support at least two transition strategies:

- **joint transition:** sample the complete state $z_j$ in one categorical
  draw;
- **coordinate transition:** update one trait indicator or component at a time
  while conditioning on the remaining state.

Coordinate updating additionally requires the declared state graph to be
connected under one-coordinate moves. Disconnected restricted spaces, such as
only $(0,0)$ and $(1,1)$, require the joint transition or another separately
derived irreducible kernel.

When both strategies represent the same state model, they must be shown to
leave the same posterior invariant. One-marker exact-probability tests should
validate both.

Joint transitions are preferred for the phase-1 BayesC reference because they
are easier to verify. Coordinate transitions are an important scalability path,
not merely an optimization.

## 8. Marker-effect covariance $V_b$

### 8.1 Required invariant

For the global model, at iteration $s$ there must be exactly one authoritative
covariance state

$$
V_b^{(s)}.
$$

The covariance used by marker transitions, stored in retained draws, included
in posterior summaries, and returned as the final chain state must refer to the
same MCMC state. A descriptive estimator must never overwrite the sampled
covariance parameter.

The corresponding invariant for a regional model is one persistent,
authoritative $V_{b,r}^{(s)}$ for every declared statistical region $r$. Each
regional state must be used consistently by its marker transitions, retained
draws, summaries, and final chain state. Computational blocks and providers do
not create additional covariance states unless the statistical model
explicitly declares them.

### 8.2 Phase-1 prior family

The recommended phase-1 reference prior is

$$
V_b\sim\operatorname{IW}(\nu_b,\Psi_b),
$$

using the convention

$$
\mathbb E(V_b)
=
\frac{\Psi_b}{\nu_b-T-1},
\qquad
\nu_b>T+1.
$$

If the user supplies a desired prior mean $V_{b0}$, the resolved scale is

$$
\Psi_b
=
(\nu_b-T-1)V_{b0}.
$$

The public and native contracts must state whether an input is a prior mean,
scale, or mode. A generic name such as `ssb_prior` is insufficient for the new
API unless accompanied by explicit parameterization metadata.

### 8.3 Cheng full-latent reference augmentation

The Cheng/JWAS reference formulation samples a complete latent
$\boldsymbol{\beta}_j$ for every marker, including fully null
markers. Under

$$
\boldsymbol{\beta}_j\mid V_b\sim N_T(0,V_b),
$$

the conjugate update is

$$
V_b\mid \boldsymbol{\beta}_1,\ldots,\boldsymbol{\beta}_M
\sim
\operatorname{IW}
\left(
\nu_b+M,
\Psi_b+
\sum_{j=1}^M
\boldsymbol{\beta}_j\boldsymbol{\beta}_j^\top
\right).
$$

This is a coherent data-augmentation scheme when inactive latent components are
actually sampled from their correct conditional distributions. Setting those
components to zero while retaining the $\nu_b+M$ update is not equivalent.

### 8.4 Production partially collapsed and conditionally completed augmentation

Fully null markers have no likelihood contribution. Their latent slab effects
can therefore be integrated out. For a non-null marker with only a subset of
traits active, inactive elements of $\boldsymbol{\beta}_j$ must still be
augmented from their
conditional normal distribution if an inverse-Wishart update based on full
outer products is used.

Let $\mathcal A$ denote markers with at least one active trait, and let

$$
M_+=|\mathcal A|.
$$

After completing the latent vector for each $j\in\mathcal A$, the update is

$$
V_b\mid\{\boldsymbol{\beta}_j:j\in\mathcal A\}
\sim
\operatorname{IW}
\left(
\nu_b+M_+,
\Psi_b+
\sum_{j\in\mathcal A}
\boldsymbol{\beta}_j\boldsymbol{\beta}_j^\top
\right).
$$

The inactive conditional is

$$
\boldsymbol{\beta}_{j,I}
\mid
\boldsymbol{\beta}_{j,A},V_b
\sim
N
\left(
V_{b,IA}V_{b,AA}^{-1}\boldsymbol{\beta}_{j,A},
V_{b,II}-V_{b,IA}V_{b,AA}^{-1}V_{b,AI}
\right).
$$

The standalone covariance prototype has derived and validated this partially
collapsed Gibbs step against full augmentation and independent enumeration. It
is not valid to embed observed active subvectors in a full vector of zeros and
treat the result as a complete draw from $N_T(0,V_b)$.

### 8.5 Recommended augmentation decision

The partially collapsed and conditionally completed scheme is the approved
production recommendation because:

- null markers need not introduce auxiliary covariance information;
- sparse models avoid sampling large numbers of prior-only latent vectors;
- the degrees of freedom reflect the number of augmented non-null vectors;
- it should reduce dependence between $V_b$ and irrelevant null-marker
  augmentations.

The required derivation, exact-reference comparison, stationary-moment checks
and sparse-inclusion mixing comparison are recorded in
`mtblr_covariance_design.md`. Cheng full augmentation remains the principal
reference sampler.

### 8.6 Later covariance-prior providers

Inverse-Wishart is recommended as the exact reference implementation, not as
the permanent only option. The covariance-policy interface should permit later
providers such as:

- separate standard deviations and a correlation matrix;
- LKJ-style correlation priors;
- Huang-Wand weakly informative covariance priors;
- structured or factor-analytic covariance;
- annotation- or component-dependent covariance models after separate
  scientific approval.

These alternatives must not be hidden behind the same arguments if they imply
different priors or samplers.

### 8.7 Regional covariance and pattern probabilities

Region-specific activity and covariance are separate policies. For a declared
region map $r(j)$, the first coherent regional extension is

$$
\boldsymbol{\delta}_j\mid r(j)=r
\sim
\operatorname{Categorical}(\boldsymbol{\Pi}_r),
$$

$$
\boldsymbol{\beta}_j\mid r(j)=r,V_{b,r}
\sim
N_T(0,V_{b,r}).
$$

Each $V_{b,r}$ is a persistent chain state updated exactly once per complete
iteration from markers in region $r$. Region-specific
$\boldsymbol{\Pi}_r$, latent-effect covariance $V_{b,r}$, regional genomic
covariance and descriptive regional summaries are not interchangeable.

The current native `sets` loop is not a coherent implementation of this model.
It owns one global covariance, immediately overwrites `sampleBset()` with an
all-marker latent draw, repeats that global draw once per set, updates one global
pattern simplex and retains only a final global heuristic covariance. Until a
regional policy is implemented, sets must be treated only as traversal
partitions and must not alter global parameter update frequency.

The first regional **reference** should provide independent proper
inverse-Wishart and Dirichlet priors, plus fixed-covariance and shared-global
reductions for testing. Independent unrestricted matrices should not become the
production default for low-information regions. Production should select a
valid partial-pooling prior or a shared fixed/pre-estimated covariance-template
library with region-specific weights after the reference passes exact
small-region tests.

Regional pattern probabilities use separately declared priors
$\mathbf a_{\Pi,r}$. A shared-global pattern state instead uses one declared
$\mathbf a_{\Pi,\mathrm{global}}$; it must not be constructed by summing
replicated regional priors.

## 9. Residual covariance $V_e$

Residual covariance is determined jointly by the statistical model and the
available data, not merely by a software option.

### 9.1 Common-sample individual-level likelihood

For complete phenotypes on one common sample,

$$
e_i\mid V_e\sim N_T(0,V_e).
$$

With

$$
V_e\sim\operatorname{IW}(\nu_e,\Psi_e),
$$

the full conditional is

$$
V_e\mid Y,B
\sim
\operatorname{IW}
\left(
\nu_e+N,
\Psi_e+(Y-XB)^\top(Y-XB)
\right).
$$

This is the phase-1 full-residual-covariance reference.

### 9.2 Diagonal residual model

Diagonal mode should be defined as independent traitwise residual-variance
models, for example

$$
v_{e,t}\sim\operatorname{ScaleInv}\chi^2
(\nu_{e,t},s_{e,t}^2).
$$

It should not be described as inference for an unrestricted covariance matrix
whose off-diagonal entries are subsequently forced to zero.

### 9.3 Summary statistics without modeled overlap

When only marginal trait summaries are supplied and cross-trait phenotype
products are unavailable, the likelihood supports traitwise residual
variances. It does not identify arbitrary residual covariance.

The public contract should therefore specify:

- `sample_overlap = "not_modeled"`;
- marginal phenotype cross-products only;
- diagonal or fixed residual variance;
- no inferred off-diagonal $V_e$.

### 9.4 Future overlap-aware summary likelihood

An overlap-aware model requires a separate likelihood contract including the
necessary cross-trait information, such as an appropriately defined phenotype
cross-product matrix $S_Y$, overlap counts or their equivalent, compatible
marker ordering and scaling, and operator provenance.

This must be a later model rather than an extra switch inside the current
no-overlap likelihood.

## 10. Likelihood regimes

The architecture should represent three regimes explicitly.

### 10.1 Regime A: common-sample individual data

Inputs include:

- one aligned genotype representation for the common individuals;
- complete centered phenotype matrix $Y$;
- one sample size $N$;
- optional pre-adjustment/covariate provenance;
- full or diagonal $V_e$ policy.

This regime identifies common-sample residual and genomic covariance.

### 10.2 Regime B: marginal summary data with no modeled overlap

Inputs include one or more declared providers $d$ for each trait $t$:

- $\mathbf s_{dt}=X_{dt}^\top\mathbf y_{dt}$;
- $\mathbf y_{dt}^\top\mathbf y_{dt}$ when required by the residual policy;
- sample size $N_{dt}$;
- provider-local cross-product operator $C_{dt}=X_{dt}^\top X_{dt}$;
- a local-to-global marker map $m_{dt}$;
- marker, allele, coding, centering, standardization, effect-scale, population,
  and summary-error provenance.

This regime supports a joint marker prior but trait-marginal likelihood
contributions. Providers may have different sample sizes, marker coverage, LD
populations, retained eigenspaces, ranks, and blocks. Cross-trait borrowing is
through the marker prior; equal $C_{dt}$ is not required. Off-diagonal residual
covariance is not identified.

### 10.3 Regime C: overlap-aware summary data

This future regime must define the joint summary likelihood from first
principles. It should not infer overlap from marginal inputs or silently reuse
the Regime B likelihood.

## 11. Likelihood-operator interface

The shared sampler should depend on a small statistical operator contract, not
on BED, CSR, or eigen implementation details.

At minimum, a likelihood operator should provide:

1. trait and marker counts;
2. marker diagonal or local precision contribution;
3. current marker score with the marker's existing effect restored;
4. application of an effect difference to the residual state;
5. reconstruction or validation of residual scores;
6. identified quadratic forms required by the likelihood;
7. declared genomic-covariance capabilities;
8. provenance and approximation status.

Conceptually:

```text
trait_count()
marker_count()
marker_score(j, state)
marker_diagonal(j)
apply_effect_delta(j, delta, state)
rebuild_residual_state(effects)
likelihood_quadratic(state)
genomic_quantity_capability()
```

The exact C++ interface may use templates for zero-overhead dispatch. The
statistical contract should nevertheless be representation-neutral and tested
independently of that implementation choice.

### 11.1 BED operator

The BED adapter owns genotype decoding and individual residuals. It can provide
exact common-sample cross-products and full residual cross-products.

### 11.2 CSR operator

The CSR adapter owns disk-backed sparse LD traversal and trait-specific or
shared operators. It provides marginal score updates and operator-relative
quadratics.

### 11.3 Retained block-eigen operator

The block-eigen adapter owns retained eigenspaces, filters, and reconstruction
rules. Approximation status must remain visible in provenance and output.

Filtered block-eigen output should not be described as exact full-data output.

Every provider owns its own block representation,

$$
C_{dt,b}
\approx
Q_{dt,b}\Lambda_{dt,b}Q_{dt,b}^\top,
$$

including its local marker map, block structure and retained rank $r_{dt,b}$.
No common eigenbasis across traits or providers is required. A full-rank
representation must reduce to the corresponding provider cross-product;
retained rank defines that provider's explicit approximate likelihood.

### 11.4 Provider collection and overlap contract

For providers estimating one shared trait effect, independent likelihood terms
are accumulated through the operator interface rather than by requiring a
materialized $\sum_d C_{dt}$. Combination is valid only after allele
orientation, genotype coding, centering, standardization, phenotype/effect
scale, marker maps, and operator scale are compatible. A marker missing from a
provider supplies no likelihood term; it is not a zero-effect observation.
Provider-local order is arbitrary when marker identifiers, scores, effect
alleles, operator rows and columns, and the local-to-global map are permuted
consistently. Reordering only a subset of these objects is misalignment, not an
equivalent provider representation.

Providers with overlapping samples require declared cross-provider
summary-score error covariance. Marginal dense, CSR, or block-eigen operators
do not identify this covariance. A future decomposition into shared and
provider-specific effects is a separate multi-cohort model, not an adapter
option.

Computational eigenblocks, statistical covariance regions, and annotation or
source groups are separate objects. Eigenblock traversal must not change the
prior region map, and global or regional prior states are updated once per
complete iteration rather than once per provider or eigenblock.

## 12. One authoritative MCMC schedule

For one complete logical chain, the phase-1 schedule should be:

1. update all marker states and latent/realized effects using the current
   $V_b$ and $V_e$;
2. update global state probabilities or annotation-prior parameters;
3. complete any inactive latent dimensions required by the selected
   augmentation policy;
4. update $V_b$ exactly once;
5. update $V_e$ exactly once when supported and requested;
6. calculate declared derived quantities;
7. retain the completed iteration when required by the common retention rule.

Marker sets, BED files, providers, LD blocks, and eigenblocks may determine
traversal and memory ownership. They must not determine how many times a
global or regional statistical parameter is updated. A global covariance is
updated once per complete iteration; each persistent regional covariance is
also updated once per complete iteration from the markers assigned to that
region.

### 12.1 Retention rule

Define a single completed-iteration retention rule. If $s$ is zero-based and
$s_{\mathrm{burn}}$ is the number of burn-in iterations, retain when

$$
s\geq s_{\mathrm{burn}}
$$

and

$$
(s-s_{\mathrm{burn}})\bmod s_{\mathrm{thin}}=0.
$$

Posterior marker means, state probabilities, covariance means, probability
means, and retained diagnostic draws must all refer to this same retained set
unless a field explicitly declares a different observational schedule.

### 12.2 Fixed-parameter paths

When `updateB`, `updateE`, `updatePi`, or an annotation update is disabled:

- the state remains fixed;
- the path consumes zero RNG draws for that disabled update;
- the trace, if requested, records the fixed value;
- the output states that the parameter was fixed rather than posterior-sampled.

## 13. Genomic covariance and related derived quantities

Three quantities must be kept distinct.

### 13.1 Common-sample genomic covariance

For one common sample,

$$
V_g^{\mathrm{sample}}
=
\frac{1}{N}B^\top X^\top XB.
$$

This is an identified covariance of fitted genetic values in the analyzed
sample, subject to the declared centering and denominator convention.

### 13.2 Target-population genomic covariance

Given a declared compatible target-population genotype covariance
$\Sigma_X^{\mathrm{target}}$,

$$
V_g^{\mathrm{target}}
=
B^\top\Sigma_X^{\mathrm{target}}B.
$$

This requires compatible marker order, alleles, coding, centering,
standardization, population definition, and a positive-semidefinite covariance
operator.

### 13.3 Operator-relative bilinear forms

With trait-specific operators $C_t$ and sample sizes $N_t$, an expression such
as

$$
\frac{1}{2\sqrt{N_tN_s}}
\left(
b_t^\top C_s b_s+b_s^\top C_t b_t
\right)
$$

may be a useful symmetric descriptive quantity. It is not automatically a
paired-sample or target-population covariance.

Such a value must have a separate name and identification status. It must not
be placed in an unqualified covariance field.

### 13.4 Output capability status

Every covariance-like derived output should carry a status such as:

- `identified_common_sample`;
- `identified_target_population`;
- `trait_diagonal_only`;
- `descriptive_operator_relative`;
- `unavailable`.

## 14. Target native architecture

The native architecture should separate the following layers.

```text
validated R model/data specification
                |
                v
      immutable prepared data
                |
                v
   provider/operator collection
                |
                +------------------+
                |                  |
                v                  v
        marker-state policy   covariance policies
                \                  /
                 \                /
                  v              v
                  shared MT chain engine
                           |
                           v
                    named raw result v2
```

### 14.1 Prepared data

Prepared data are immutable and may be shared across logical chains. They own
or reference:

- provider-specific genotype/LD/eigen resources;
- aligned provider-local score vectors and phenotype cross-products;
- global marker universe and local provider maps;
- marker, trait, provider, allele, scale, population, and overlap metadata;
- state descriptors;
- prior descriptors;
- traversal sets;
- operator provenance.

### 14.2 Chain state

Each chain owns mutable plain-C++ state:

- residual or residual-score state;
- latent base effects $U$;
- realized effects $B$;
- marker states $z$;
- component indices;
- $V_b$ and its inverse/factorization;
- $V_e$ and its inverse/factorization;
- probability or annotation parameters;
- counters and accumulators;
- chain-local RNG.

No R API or Rcpp object should be touched inside an OpenMP worker region.

### 14.3 Marker-state policy

The policy owns state enumeration or coordinate definitions, prior
probabilities, marker conditional calculations, state draws, and state
sufficient statistics.

### 14.4 Covariance policies

Separate policies should own:

- $V_b$ prior validation and transition;
- $V_e$ prior validation and transition;
- fixed versus sampled behavior;
- covariance-specific diagnostics;
- prior parameterization metadata.

### 14.5 Result builder

The chain engine should finish in plain C++ structures. Rcpp conversion and raw
schema construction should happen only after parallel execution has ended.

## 15. Raw result schema version 2

An intentional `mtblr_raw` version-2 migration is recommended. It should follow
the package's named namespace architecture while avoiding a simultaneous
single-trait schema overhaul.

Suggested top-level namespaces are:

```text
schema
meta
posterior
draws
final_state
derived
model
data
alignment
diagnostics
chains
```

### 15.1 Posterior summaries

`posterior` should contain posterior means or probabilities only, for example:

```text
posterior$marker$effect_mean
posterior$marker$activity_probability
posterior$marker$component_probability
posterior$covariance$vb_mean
posterior$covariance$ve_mean
posterior$covariance$vg_mean
posterior$probability$state_mean
```

### 15.2 Retained draws and diagnostic traces

`draws` should make the observational unit explicit:

```text
draws$iteration
draws$vb
draws$ve
draws$vg
draws$state_probability
draws$annotation_parameters
draws$selected_markers
```

An iterationwise mean across chains must not masquerade as a posterior draw.
Convergence input should preserve a chain dimension.

### 15.3 Final state

`final_state` contains the terminal state of one declared chain:

```text
final_state$chain_index
final_state$marker$latent_effect
final_state$marker$realized_effect
final_state$marker$state
final_state$covariance$vb
final_state$covariance$ve
final_state$probability
```

If pooled summaries use all chains but final state comes from chain 1, that
policy must be structural metadata rather than an implicit convention.

### 15.4 Derived quantities

`derived` should contain induced quantities and their identification metadata:

```text
derived$genomic_covariance
derived$genomic_covariance_status
derived$operator_relative_bilinear
derived$le_ld_decomposition
```

### 15.5 Formatted-fit stability

The raw schema may improve substantially while the canonical formatted fit
continues to expose stable names where their meanings remain scientifically
valid:

```text
bm, dm, wy, r, b, d,
vbs, vgs, ves, vle, vld,
pi_trace, pi_final, pi_mean,
cov_b_mean, cov_g_mean, cov_e_mean,
cov_b_final, cov_g_final, cov_e_final,
input, data, diagnostics,
convergence, convergence_traces, chains, memory_estimate
```

However, a stable name must not preserve an invalid meaning. If an off-diagonal
`cov_g_mean` is not identified, the migration should return an explicit
unavailable value and status or introduce a clearly named descriptive field.
The exact formatted-fit migration requires a separate field-by-field contract.

## 16. Public API direction

The public entry points should remain:

- `mtblr_bed()`;
- `mtblr_csr()`;
- `mtblr_block_eigen()`.

Their shared model controls should be resolved into explicit internal
specifications before native execution:

```text
data_spec
likelihood_spec
state_policy_spec
marker_covariance_spec
residual_covariance_spec
execution_spec
retention_spec
output_spec
```

### 16.1 Covariance prior arguments

The new API should avoid ambiguous prior arguments. One possible direction is:

```text
vb_init
vb_prior_mean
vb_prior_df
ve_init
ve_prior_mean
ve_prior_df
```

or an explicit structured object:

```text
vb_prior = list(
  family = "inverse_wishart",
  df = ...,
  mean = ...
)
```

The structured approach is more extensible, while explicit scalar arguments
are simpler. This remains an API decision.

### 16.2 State-policy arguments

The current model names should identify the prior family, while a separate
control identifies state geometry when more than one geometry is supported:

```text
state_policy = "shared_component"
```

and later:

```text
state_policy = "trait_specific_component"
```

The package should not assign both different statistical models the same name
without recording the policy.

## 17. Lessons adopted from JWAS

The redesign should adopt or evaluate the following JWAS ideas:

1. explicit separation of latent base effects, realized effects, and marker
   states;
2. one central MCMC schedule for marker, covariance, residual, and probability
   updates;
3. full latent augmentation as a valid reference transition;
4. joint and coordinate MT state samplers;
5. trait-specific BayesR magnitude classes as an optional generalized state
   policy;
6. marker-specific joint priors for annotated MT models;
7. exact one-marker probability tests and empirical transition-frequency
   checks;
8. dense-versus-block transition comparisons.

The redesign should not copy without independent validation:

- implicit covariance-scale conventions;
- two-trait-only restrictions;
- approximate independent-block semantics;
- output conventions that do not distinguish means, draws, and final states;
- any external behavior that conflicts with approved `sblr` theory.

## 18. Testing and validation strategy

### 18.1 Analytical unit tests

Permanent tests should cover:

- state descriptor ordering and normalization;
- one-marker BayesC joint-state probabilities;
- one-marker BayesR shared-component probabilities;
- joint versus coordinate transition frequencies;
- latent-to-realized effect transformation;
- conditional distribution of inactive latent dimensions;
- inverse-Wishart scale and degrees of freedom;
- prior mean/scale conversion;
- full and diagonal $V_e$ conditionals;
- fixed-update zero-RNG behavior.

### 18.2 Augmentation validation

For full-latent and partially collapsed null augmentation:

- derive both kernels independently;
- compare posterior moments on tiny simulated problems;
- compare prior-only stationary moments;
- verify the null-state limit;
- compare autocorrelation and effective sample size in sparse settings;
- ensure active-subspace posterior quantities agree within Monte Carlo error.

### 18.3 Partition invariance

With identical seed and an update scheme constructed for controlled replay,
verify that changing computational sets does not change:

- the number of global covariance updates;
- the posterior target;
- analytical sufficient statistics.

Exact trajectory equality may not be expected if marker traversal order changes.
The test should distinguish target invariance from controlled-replay identity.

### 18.4 Operator reductions

Use exact small fixtures to compare:

- packed BED and its exact cross-product representation under matching
  likelihood contracts;
- dense/full-rank LD and CSR;
- CSR and unfiltered full-rank block eigen;
- filtered block eigen against its explicitly projected likelihood;
- one-trait MT reduction against the corresponding ST model.

### 18.5 Multichain and RNG tests

Verify:

- repeatability with explicit chain seeds;
- serial versus parallel chain equality;
- chain retention independence;
- one prepared immutable operator shared across chains;
- no R RNG or R API use inside workers;
- no RNG consumption by disabled updates or diagnostics.

### 18.6 Known-truth simulations

The first simulation ladder should include:

1. $T=1$ reduction;
2. $T=2$ independent traits with diagonal $V_b$ and $V_e$;
3. $T=2$ correlated marker effects with diagonal $V_e$;
4. $T=2$ correlated residuals with diagonal $V_b$;
5. shared and trait-specific causal patterns;
6. sparse and moderately polygenic architectures;
7. null data;
8. weakly and strongly identified covariance settings.

Evaluation should distinguish:

- posterior calibration;
- covariance recovery;
- marker-state recovery;
- prediction;
- convergence and mixing;
- computational scaling.

Large benchmarking belongs in `sblrbench`, not in the package unit-test suite.

## 19. Phased implementation roadmap

### Phase 0: approve the statistical contract

Deliverables:

- approved $V_b$ augmentation and prior parameterization;
- approved $V_e$ policies;
- approved complete Cheng MT-BayesC$\Pi$ state model;
- approved update and retention schedule;
- approved genomic-covariance naming and identification policy;
- approved raw-schema migration direction;
- updated canonical Methods derivations where required.

The covariance-design prototype resolves the $V_b$ augmentation, complete
$T=2$ state geometry, pattern-probability and joint-transition portions of this
checkpoint. Residual-covariance and public migration decisions remain separate.

### Phase 1: individual-level MT BayesC reference

Phase 1 also fixes the representation-neutral scientific-kernel and likelihood
operator contracts needed by later BED, CSR, and block-eigen adapters. It does
not require the production summary-statistic adapters to be implemented in the
same phase.

Implement:

- packed-BED common-sample likelihood;
- complete centered $Y$;
- all four $T=2$ activity patterns and one declared joint
  $\boldsymbol{\Pi}$ simplex;
- joint BayesC state transitions;
- null collapse with conditional completion for partially active patterns;
- one authoritative inverse-Wishart $V_b$ update;
- the approved initial $V_e$ policy, with full and diagonal policies included
  only after decision D3 is resolved;
- one logical chain;
- raw schema version 2;
- exact analytical tests and tiny known-truth simulations.

This is the scientific production reference. A transparent Cheng full-latent
sampler remains the independent reference formulation used to validate it.

### Phase 2: multichain execution and formatted-fit migration

Implement:

- chain-local RNG;
- static OpenMP chain dispatch;
- pooled posterior summaries;
- chain-preserving convergence draws;
- compact retained chain records;
- final formatted-fit contract and `check_stblr_consistency()` updates;
- memory diagnostics.

### Phase 3: shared-component MT BayesR

Implement:

- deterministic pattern-by-component state descriptor;
- component scaling through a separate $G_j$ while $D_j$ remains the activity
  mask;
- exact small-state joint transition;
- scalable coordinate transition if required;
- component posterior summaries;
- fixed MAF-$S$ multiplier as a separate $q_j$ policy;
- BayesC reduction with one unit positive component.

### Phase 4: MT BayesRC

Port the approved current probability architecture only after BayesR is stable:

- annotation-dependent component probability $P$;
- global conditional non-null sharing $H$;
- current probit-stick coefficient hierarchy;
- exact marker alignment;
- independent tests of annotation and marker transitions.

Do not add annotation-dependent $V_b$, $H$, or $Q$ in this phase unless a
separate theory extension is approved.

### Phase 5: summary statistics without modeled overlap

Implement:

- heterogeneous independent-provider adapters using the operator contract
  fixed in Phase 1;
- provider-local marker maps and dense or CSR likelihood operators;
- diagonal/fixed residual policy;
- trait-specific or shared LD resources;
- explicit no-overlap contract;
- identified diagonal genomic variances;
- clearly named operator-relative cross-trait quantities when useful.

### Phase 6: retained block-eigen likelihood

Implement:

- the same scientific sampler through provider-specific block-eigen
  operators;
- exact unfiltered reductions;
- filtered-likelihood provenance;
- block ownership and memory diagnostics;
- no block-dependent global parameter update frequency.

### Phase 7: generalized state and covariance providers

Candidates include:

- trait-specific BayesR component classes;
- alternative covariance priors;
- structured trait-pattern priors;
- larger-$T$ coordinate or sparse state transitions;
- a persistent region-specific $V_{b,r}$ and $\boldsymbol{\Pi}_r$ reference,
  beginning with independent proper priors and fixed/global reductions, then a
  separately approved partial-pooling production policy;
- fixed covariance-template libraries and region-specific template weights as
  a later mash-like extension.

### Phase 8: overlap-aware summary likelihood

Proceed only after a separate statistical design specifies the required
cross-trait sufficient statistics and identification assumptions.

### Consolidated priority layers

The phase labels above describe implementation dependencies. The long-term
scientific roadmap is grouped as follows so that heterogeneous summary
providers are not postponed behind optional prior extensions.

**Near-term core**

1. MT-BRR as a dense/common-sample covariance reference.
2. Correct Cheng MT-BayesC$\Pi$ with all declared activity patterns.
3. A joint categorical pattern transition.
4. Null collapse and conditional completion of partially active vectors.
5. One authoritative sampled $V_b$ per completed iteration.
6. A common-sample individual-level adapter.
7. Heterogeneous independent-provider summary adapters.
8. Provider-specific dense, CSR, and block-eigen operators, including distinct
   marker maps, sample sizes, block structures, eigenvectors, eigenvalues, and
   retained ranks.

**Controlled extensions**

9. Pattern-by-scale MT-BayesR using the same conditional-completion logic and
   an explicitly selected base- or scaled-effect storage convention.
10. Explicit set-, source-, and region-specific prior models.
11. Fixed or pre-estimated covariance-template libraries.
12. Region-specific weights over shared covariance templates.
13. An overlap-aware summary likelihood with declared cross-provider error
    information.

**Larger-trait extensions**

14. Restricted activity-pattern sets.
15. Structured trait groups.
16. Low-rank plus sparse marker effects,

    $$
    \boldsymbol{\beta}_j
    =
    \Lambda\mathbf f_j+\boldsymbol{\epsilon}_j.
    $$

17. Factor or structured-pattern transitions when enumeration of $2^T$ states
    is infeasible.

Research extensions include global-local shrinkage, sparse precision or
graphical covariance models, multi-population effects, multi-environment or
longitudinal effects, multivariate single-effect fine-mapping, and
nonparametric covariance-template mixtures. These are not immediate
implementation commitments.

## 20. Documentation migration

An approved implementation must update together:

- `docs/methods/` statistical definitions;
- `docs/dev/` architecture, schema, contracts, and migration records;
- `docs/notes/statistical_quantities_and_sblr_outputs.qmd`;
- affected practical MT workflows;
- capability and roadmap pages;
- roxygen documentation;
- raw-schema and formatted-fit tests;
- benchmark scenario definitions in `sblrbench`.

Methods pages must use statistical notation rather than current fit-field names.
Mappings to public arguments and formatted fields belong in Notes and developer
contracts.

## 21. Decision register

### D1. Phase-1 marker-effect covariance prior

**Recommendation:** inverse-Wishart with explicit mean-to-scale conversion.

**Status:** approved by the covariance-design derivation and prototype for the
reference transition.

### D2. Latent augmentation for $V_b$

**Options:**

1. full-latent JWAS-style augmentation over all markers;
2. partially collapsed null augmentation with completed latent vectors for
   non-null markers.

**Recommendation:** retain Cheng full augmentation as the reference; use
partially collapsed null augmentation with conditional completion for every
non-null marker in production.

**Status:** resolved by exact enumeration, numerical reference cases and mixing
comparisons in `mtblr_covariance_design.md`.

### D3. Phase-1 residual covariance

**Recommendation:** support both full inverse-Wishart and explicitly diagonal
traitwise residual models for common-sample BED data.

**Status:** awaiting approval.

### D4. Phase-1 marker-state transition

**Recommendation:** complete joint BayesC pattern draw for the reference
implementation.

**Status:** approved for the complete Cheng pattern space. Coordinate updates
target the same posterior only on a one-coordinate-connected declared state
graph and require explicit mixing qualification.

### D5. BayesR state geometry

**Recommendation:** implement current shared-component `sblr` geometry first;
retain an extension point for JWAS-inspired trait-specific components. The
standalone prototype stores the scaled vector $\boldsymbol{\theta}_j$ and
removes $\gamma_kq_j$ exactly once in the covariance statistic. A production
base-vector implementation is equivalent only if it does not remove that scale
again.

**Status:** pattern-by-scale geometry and scale convention approved by the
covariance design; production implementation remains future work.

### D6. Genomic covariance outputs for trait-specific summary operators

**Recommendation:** return identified diagonal genomic variances; do not label
cross-operator bilinear forms as covariance; provide a separately named
descriptive quantity only if scientifically useful.

**Status:** awaiting approval.

### D7. Raw schema migration

**Recommendation:** intentional `mtblr_raw` version 2 with separate posterior,
draw, final-state, and derived namespaces.

**Status:** awaiting approval.

### D8. Public covariance-prior API

**Options:** explicit `*_prior_mean` and `*_prior_df` arguments, or structured
prior objects.

**Recommendation:** structured internal specification is required; public API
choice remains open.

**Status:** unresolved API decision.

### D9. Region-specific covariance

**Recommendation:** no current support claim. Use persistent independent
$V_{b,r}$ and $\boldsymbol{\Pi}_r$ as the reference, with fixed and
shared-global reductions. Update every persistent regional covariance once per
iteration. Select partial pooling or shared covariance templates before making
unrestricted regional matrices a production default.

The shared-global reduction uses one explicitly declared
$\mathbf a_{\Pi,\mathrm{global}}$; it does not sum replicated regional prior
vectors.

**Status:** statistical transition validated in the standalone prototype;
production and schema design remain future work.

### D10. Covariance-template mixtures

**Recommendation:** place mash-like fixed covariance templates after the Cheng
and regional reference implementations. Singular trait-specific templates must
use active-subspace calculations or explicit masks.

**Status:** planned scientific extension, not current functionality.

### D11. Heterogeneous summary providers

**Recommendation:** make provider-local marker maps and dense, CSR, or
block-eigen operators part of the near-term core. Independent provider
likelihoods may differ in $C_{dt}$ and $N_{dt}$; overlapping providers require
a separate declared summary-error covariance contract.

**Status:** architecture and standalone operator reductions validated in
`mtblr_covariance_design.md`; production adapter design remains future work.

## 22. Acceptance criteria for the design phase

The design phase is complete only when:

1. every phase-1 probability model and covariance conditional has an approved
   derivation;
2. prior parameterizations are unambiguous;
3. the augmentation policy is shown to target the intended posterior;
4. update order and retained-draw semantics are fixed;
5. each likelihood regime has explicit required inputs and identifiable
   outputs;
6. the operator interface is sufficient for BED, CSR, and block eigen without
   representation-specific scientific branches;
7. schema version 2 has a field-level contract;
8. the migration plan identifies reusable, replaced, and retired code;
9. analytical, operator-reduction, RNG, and known-truth tests are specified;
10. heterogeneous provider maps, scale compatibility, independence, and
    provider-specific block-eigen provenance are explicit;
11. computational eigenblocks and statistical covariance regions are distinct;
12. current and proposed capabilities are not conflated.

## 23. Immediate next design task

The covariance research-design checkpoint is complete for shared and complete
$T=2$ Cheng patterns and for the heterogeneous independent-provider operator
abstraction. This does not yet mean that every acceptance criterion for the
full MTBLR design has been met: decisions D3, D6, D7, and D8 still contain
contracts needed by production code.

The next task should specify and implement a narrowly gated production vertical
slice without changing the established posterior target. Its contract stage
must precede edits to the production sampler:

1. approve the public inverse-Wishart mean/scale migration;
2. resolve the initial common-sample $V_e$ policy required by D3;
3. define the version-2 raw-schema representation for actual $V_b$ and
   $\boldsymbol{\Pi}$ draws and settle the minimum public prior interface
   required by D7 and D8;
4. record the D6 naming and identification boundary for genomic covariance
   outputs, even if summary-statistic outputs are implemented later;
5. implement the individual-level joint pattern transition with null collapse
   and conditional completion;
6. retain the full-latent standalone sampler as an independent oracle;
7. define the provider-neutral adapter contract and provider-local maps during
   the vertical slice, then implement dense/CSR and block-eigen adapters in
   their declared phases without assuming common cross-products;
8. add exact pattern, empty-active, BayesR scale, provider-order,
   full-rank-eigen reduction and traversal-partition tests;
9. defer region-specific covariance and mash-like templates to separately
   gated milestones.

## 24. References

- Cheng H, Kizilkaya K, Zeng J, Garrick D, Fernando R. Genomic prediction from
  multiple-trait Bayesian regression methods using mixture priors. *Genetics*.
  2018;209:89--103.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC5937171/>.
- `sblr` canonical Methods pages at the pinned baseline, especially
  `model_theory.qmd`, `multitrait_overlap.qmd`, `mt_bayesr_sbayesr.qmd`, and
  `mt_bayesrc_sbayesrc.qmd`.
- `sblr` developer contracts at the pinned baseline, especially
  `blr_architecture.md`, `blr_output_schema.md`, `blr_model_contracts.md`,
  `blr_convergence_contract.md`, and
  `methods_source_traceability_audit.md`.
- JWAS multi-trait BayesC sampler:
  [`MTBayesABC.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl).
- JWAS multi-trait BayesR sampler:
  [`MTBayesR.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/markers/BayesianAlphabet/MTBayesR.jl).
- JWAS covariance updates:
  [`variance_components.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/variance_components.jl).
- JWAS central MCMC schedule:
  [`MCMC_BayesianAlphabet.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl).
- Barnard J, McCulloch R, Meng X-L. Modeling covariance matrices in terms of
  standard deviations and correlations, with application to shrinkage.
  *Statistica Sinica*. 2000;10:1281--1311.
- Huang A, Wand MP. Simple marginally noninformative prior distributions for
  covariance matrices. *Bayesian Analysis*. 2013;8:439--452.
