# SBayesRC-S C++ implementation and integration record

## Scope and baseline

Phase 4 starts from commit `6528ae15118f1ad4a3fba6661583119e1f9f848b`,
which contains the validated Phase-1, Phase-2, and Phase-3 standalone R
hierarchy. Standard SBayesRC is a separate frozen model. This record is
written before native implementation and is updated only with evidence that
has actually passed.

## Architecture audit before implementation

### Existing ownership

The neutral annotation/probit mechanics live primarily in
`src/st_bayesrc_annotation_prior.h`:

- stable standard-normal probability helpers;
- exact truncated-normal sampling;
- allocation-to-eligible-stick indicators;
- probit-stick q-to-component probability transformation;
- the `std::mt19937` RNG contract;
- scalar Gaussian conditional helpers.

The standard continuous-alpha hierarchy, including `sigmaSqAlpha`, is also in
that header. It is called directly by the BED core
(`blr_bed_bayesrc_core_impl.h`), the CSR/block operator core
(`blr_csr_sbayesrc_core_impl.h`), and multi-trait cores. The marker loops,
effect/component conditionals, residual maintenance, LD operators,
GCTB-compatible block residual policy, OpenMP task ownership, accumulation,
and raw-result assembly remain route-specific existing infrastructure.

### Reuse decision

Phase 4A can reuse unchanged:

- annotation storage and row indexing;
- eligible-stick meaning;
- truncated-normal primitive;
- standard RNG engine;
- q-to-component probability transformation;
- Armadillo small-matrix linear algebra;
- internal Rcpp validation-hook convention.

SBayesRC-S must newly own:

- one shared `delta[j]` across sticks;
- exact-zero excluded slopes;
- global `pi_A` and its beta update;
- stick-specific `tau2[k]` and inverse-gamma update;
- collapsed annotation Bayes factors and immediate slope regeneration;
- selected-model blocked Gaussian redraw;
- annotation-selection traces and summaries.

### Refactoring decision

No standard genomic loop or standard hierarchy is refactored for Phase 4A.
A separate selection-prior header is the smallest safe design. It composes
the neutral helpers above and is exercised through internal, non-exported
from-NAMESPACE Rcpp validation hooks. This keeps standard SBayesRC source and
RNG paths unchanged except for compilation of code that its execution path
does not call.

The existing genomic cores do not yet consume a generic annotation-prior
interface. Phase 4B therefore requires a separate post-4A decision: either a
minimal compile-time hierarchy dispatch inside one chosen core, or a small
separate backend if dispatch would risk the frozen standard path. No Phase 4B
code may be written until every Phase-4A gate passes.

## Phase 4A implementation status

`src/st_bayesrc_annotation_selection.h` implements the validated hierarchy as
an isolated model-specific component. Internal Rcpp hooks expose deterministic
mathematics and an observed-d hierarchy runner for tests; neither is exported
from the package namespace.

Deterministic parity passed for eta, eligible-row s and t, log Bayes factors,
inclusion probabilities, alpha conditional moments, beta and inverse-gamma
conditional parameters, q, and component probabilities. The focused suite
passed 30 assertions. Longer independent-seed posterior comparisons gave:

| fixture | max PIP error | max q error | max component-probability error | max normalized alpha error | max tau mean error |
|---|---:|---:|---:|---:|---:|
| 3 annotations | 0.00573 | 0.00521 | 0.00178 | 0.0262 | 0.00215 |
| 12 annotations | 0.0213 | 0.00620 | 0.00620 | 0.0599 | 0.00339 |

The all-included hierarchy, all-excluded state, zero column, annotation
permutation, and duplicate-column guards passed. The full package test suite
then passed with zero failures, the established opt-in reproducibility skip,
and the established covariance warning.

**SBS4A-R1: the isolated C++ SBayesRC-S annotation hierarchy matches the
validated R reference.**

## Phase 4B integration status

Phase 4B started only after SBS4A-R1. The chosen first route was the CSR
summary-statistic engine because it already owns the reusable BayesRC marker
loop, component/effect conditional, streamed LD operator, residual and
variance updates, OpenMP task ownership, SNP accumulation, diagnostics, and
raw-result assembly. A disabled-by-default policy object was added to that
shared core. Standard SBayesRC passes an inert policy and follows its original
hierarchy function and RNG stream. An internal-only CSR binding composes the
same engine with the selection hierarchy. No public model dispatch or
NAMESPACE entry was added.

The compact smoke, fixed all-included genomic bridge, all-excluded bridge,
probability/output guards, and disabled-policy RNG regression passed. With the
hierarchy frozen, the all-included selection path and standard SBayesRC path
were bit-identical for marker summaries, variance traces, and component
summaries.

The preregistered 160-marker selection-enabled screen then reached a valid
genomic allocation state with an empty later-stick eligible set. The Phase-3
model requires a flat always-included intercept. With zero eligible
observations its conditional is improper: the precision contribution is zero
and there is no proper prior contribution. The implementation therefore
stopped with:

> SBayesRC-S flat intercept is undefined for an empty eligible stick

This is not a likelihood, LD, residual, or numerical failure. It is a support
boundary that the isolated observed-d fixtures did not exercise. Substituting
the standard model's proper intercept prior, retaining an old intercept,
clamping, or sampling from an arbitrary fallback would change the validated
SBayesRC-S posterior, so none was introduced. The selection-enabled
multichain, correlated-annotation, and runtime qualification screens were not
continued after this hard failure.

**SBS4B-R3: genomic integration discrepancy.** Phase 4B is not qualified and
there is no production SBayesRC-S backend. The overall Phase-4 result is
**SBS4-R4**.

The next task must be a mathematical Phase-3 model revision deciding on a
proper intercept prior (and revalidating the R oracle) or another explicitly
proper empty-stick treatment. C++/genomic qualification must not resume until
that revised posterior is separately validated.
