# Annotation-prior architecture: current and proposed capability matrix

## Purpose and status vocabulary

This matrix records the implementation observed at the audit baseline
(`2123699a9cc2e91059e7d81a745420b14eca7f6e`). It is not inferred solely from
function names or documentation: public dispatch, native entry points, raw
schemas, formatter paths, and tests were cross-checked.

Status labels have the following meanings:

- **IMPLEMENTED + TESTED**: a production-facing route exists and focused tests
  exercise its scientific or reduction contract.
- **IMPLEMENTED BUT ROUTE-LIMITED**: implemented and tested, but only for the
  operator/trait combinations shown.
- **IMPLEMENTED BUT NOT FULLY QUALIFIED**: executable code exists, but the audit
  did not find a complete scientific qualification contract for the claimed
  combination.
- **DOCUMENTED / CONCEPTUAL ONLY**: described in theory or development material,
  but not a production capability.
- **PROPOSED**: recommended future architecture or implementation.
- **UNSUPPORTED BY DESIGN**: deliberately excluded from the current model or
  from an explicitly scoped future phase.

`BED`, `CSR`, and `eigen` mean packed individual-level genotypes, sparse
summary-statistics LD, and retained block-eigen summary LD, respectively.
The `S` prefix is a data-route convention: SBayesC and SBayesR use the same
prior families as BayesC and BayesR with a summary-statistics likelihood.

Theory uses \(v_b,v_g,v_e\) for scalar ST variances and \(V_b,V_g,V_e\) for
MT covariance matrices. Package fields remain explicitly named
`vb`/`covb`, `vg`/`covg`, and `ve`/`cove`; case-only R names are not proposed.

## Current model-family and operator support

| Prior family / architecture | ST BED | ST CSR | ST eigen | MT BED | MT CSR | MT eigen |
|---|---|---|---|---|---|---|
| BayesC/SBayesC, global probability and variance | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED |
| BayesR/SBayesR, global mixture probabilities | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED |
| BayesRC/SBayesRC, probit-stick component probabilities | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED |
| BayesC-LV/SBayesC-LV, learned annotation variance | DOCUMENTED / CONCEPTUAL ONLY | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | PROPOSED | PROPOSED | PROPOSED |
| BayesR-LV/SBayesR-LV, learned annotation variance | DOCUMENTED / CONCEPTUAL ONLY | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | PROPOSED | PROPOSED | PROPOSED |
| Fixed marker probability and/or variance | DOCUMENTED / CONCEPTUAL ONLY | IMPLEMENTED BUT ROUTE-LIMITED | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |
| Disjoint-group probability and/or variance | DOCUMENTED / CONCEPTUAL ONLY | IMPLEMENTED BUT ROUTE-LIMITED | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |
| Learned-logistic probability and/or variance | DOCUMENTED / CONCEPTUAL ONLY | IMPLEMENTED BUT ROUTE-LIMITED | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |
| External marker variance multiplier as a first-class provider | DOCUMENTED / CONCEPTUAL ONLY | IMPLEMENTED BUT ROUTE-LIMITED through `fixed_marker` | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |
| Annotation-dependent trait-sharing probabilities | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |
| Annotation-dependent cross-trait covariance | UNSUPPORTED BY DESIGN | UNSUPPORTED BY DESIGN | UNSUPPORTED BY DESIGN | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |

The executable capability resolver is not currently an authoritative registry.
It omits LV and reports MT BayesRC/SBayesRC as unsupported even though public
R routes, native implementations, and tests exist. Conversely, the older
developer capability page omits LV and contains claims that no longer match
dispatch. A generated, test-checked registry is therefore part of the proposed
work.

## Current annotation/prior architecture by quantity

| Mechanism | Probability architecture \(P_j\) | Variance architecture \(Q_j\) | Learning | Identification / bounds | Current scope |
|---|---|---|---|---|---|
| Global BayesC | one global \(\pi\) | \(q_j=1\) | Beta update for \(\pi\) | global scale in \(v_b\) | ST/MT; BED/CSR/eigen |
| Global BayesR | ST: one global component vector; MT: one `joint_pi` simplex over complete joint states, with component/pattern marginals | \(q_j=1\), component ladder \(\gamma_k\) | global simplex update | global scale in \(v_b\) or \(V_b\) | ST/MT; BED/CSR/eigen |
| `fixed_marker` | supplied \(\pi_j\), or centered logistic map from \(A\) and fixed coefficients | supplied \(q_j\), or centered log-linear map from \(A\) and fixed coefficients | none | \(\pi_j\) and \(q_j\) are clamped; log predictor is mean-centered before clamping | ST CSR BayesC only |
| `group` | one learned \(\pi_g\) per disjoint group | fixed or sampled \(q_g\) | Beta and scaled-inverse-chi-square updates | optional marker-count-weighted arithmetic-mean normalization | ST CSR BayesC only |
| `learned_logistic` | \(\operatorname{logit}(\pi_j)=\operatorname{logit}(\pi)+\operatorname{center}(A_j\eta_\pi)\) | \(q_j=\exp\{\operatorname{clamp}(\operatorname{center}(A_j\eta_{vb}))\}\) | random-walk MH for coefficients; logit-scale slice update for global \(\pi\), with an exact zero-offset Beta reduction | unclipped logistic probability; hard variance bounds | ST CSR BayesC only |
| Probit-stick BayesRC/SBayesRC | annotation-dependent component prior probabilities from continuation probits | \(q_j=1\), component ladder \(\gamma_k\) | latent-probit Gibbs updates for \(\alpha\) and variance hyperparameters | optional explicit intercept; probability floor and row normalization | ST/MT; BED/CSR/eigen |
| Log variance (LV) | BayesC global \(\pi\), or BayesR global mixture probabilities | \(q_j=\exp(X_j\theta)\) | elliptical slice sampling | no intercept; centered columns give geometric mean one; no clipping | ST CSR/eigen BayesC and BayesR |
| MAF-S | unchanged | \(q_j=h_j^{S+1}\) | fixed \(S\) or bounded random-walk MH where supported | not log-centered in the historical implementation | summary ST mixtures; fixed MT mixture routes; route-dependent |

## Fixed, learned, and informative information: current support

| Prior quantity | Global | Fixed marker vector | Fixed \(A\)+coefficients | Jointly learned coefficient | Informative learned prior | Fixed baseline + learned correction |
|---|---|---|---|---|---|---|
| BayesC inclusion \(\pi_j\) | IMPLEMENTED + TESTED | IMPLEMENTED BUT ROUTE-LIMITED | IMPLEMENTED BUT ROUTE-LIMITED | IMPLEMENTED BUT ROUTE-LIMITED (`learned_logistic`) | DOCUMENTED / CONCEPTUAL ONLY beyond zero-centered isotropic prior | PROPOSED |
| BayesR/RC component probability | IMPLEMENTED + TESTED | DOCUMENTED / CONCEPTUAL ONLY | IMPLEMENTED + TESTED through fixed `alpha` | IMPLEMENTED + TESTED through probit-stick `alpha` | IMPLEMENTED + TESTED with proper intercept/nonintercept priors, but covariance is limited | PROPOSED |
| Relative variance \(q_j\) | IMPLEMENTED + TESTED | IMPLEMENTED BUT ROUTE-LIMITED | IMPLEMENTED BUT ROUTE-LIMITED | IMPLEMENTED + TESTED in LV; older MH path is route-limited | PROPOSED (`theta_prior_mean`, named `theta_prior_sd`) | PROPOSED |
| MT sharing probabilities | IMPLEMENTED + TESTED globally | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |
| MT covariance \(V_b\) | IMPLEMENTED + TESTED | not a marker-specific quantity | not applicable | DOCUMENTED / CONCEPTUAL ONLY for annotation dependence | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY |

## MAF-S support qualification

MAF-S is an effect-variance architecture, not an annotation probability model.
Current code and tests support fixed and sampled forms on selected ST summary
routes, including scale-aware BayesC/BayesR and SBayesRC paths. MT mixture
models support fixed marker scale derived from MAF across their established
operators; sampled MT \(S\) is explicitly unsupported. BED and BayesC/MT
coverage are not uniform. Until the capability registry is repaired, callers
must rely on route validation rather than treating the presence of a
`maf_effect_s` argument as proof of support.

## Proposed target matrix

The target is not every mathematical combination. It is a small set of
orthogonal providers with explicit scientific qualification.

| Target capability | ST BED | ST CSR | ST eigen | MT BED | MT CSR | MT eigen | Qualification requirement |
|---|---|---|---|---|---|---|---|
| Global BayesC/R | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | trajectory regression controls |
| Common annotation design object | PROPOSED | PROPOSED | PROPOSED | PROPOSED | PROPOSED | PROPOSED | exact cross-route preprocessing fixtures |
| Fixed probability provider | PROPOSED | PROPOSED | PROPOSED | PROPOSED | PROPOSED | PROPOSED | reduction to global and existing fixed-marker routes |
| Logistic inclusion provider | PROPOSED | PROPOSED | PROPOSED | later | later | later | posterior and induced-\(P\) diagnostics |
| Probit-stick component provider | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | preserve established BayesRC trajectories |
| Fixed external variance provider | PROPOSED | PROPOSED | PROPOSED | PROPOSED | PROPOSED | PROPOSED | strict positivity, normalization, calibration, operator parity |
| Canonical LV provider | PROPOSED | IMPLEMENTED + TESTED | IMPLEMENTED + TESTED | later | PROPOSED first | PROPOSED first | informative-prior and MT shared-\(q\) oracles |
| MAF-S variance provider | later | IMPLEMENTED BUT ROUTE-LIMITED | IMPLEMENTED BUT ROUTE-LIMITED | route-dependent | IMPLEMENTED BUT ROUTE-LIMITED fixed form | IMPLEMENTED BUT ROUTE-LIMITED fixed form | explicit normalization/migration decision |
| External \(q\) + LV correction | PROPOSED | PROPOSED | PROPOSED | later | PROPOSED | PROPOSED | exact log-additive composition and calibration |
| Probit-stick probability + LV variance | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | separate identifiability and mixing study; no routine API exposure |
| Annotation-dependent sharing | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | DOCUMENTED / CONCEPTUAL ONLY | later research | later research | later research | new statistical model and calibration theory |
| Annotation-dependent covariance | UNSUPPORTED BY DESIGN | UNSUPPORTED BY DESIGN | UNSUPPORTED BY DESIGN | later research | later research | later research | strong identifiability constraints and dedicated validation |

The first MT-LV target is specifically MT BayesR/SBayesR with its existing
global joint-state simplex unchanged and a shared q. Its component and pattern
probabilities remain marginals, not newly factorized P/H parameters. MT
BayesC-LV is a later possible extension.
BayesRC/SBayesRC plus LV is not part of that target and remains
DOCUMENTED / CONCEPTUAL ONLY pending the separate probability-plus-variance
validation phase.

## Current-versus-target implementation ownership

| Concern | Current ownership | Target ownership |
|---|---|---|
| Annotation preparation | separate fixed/group/logistic, LV, and BayesRC helpers | one immutable `AnnotationDesign`, with model-specific intercept policy declared by the provider |
| Marker probability | embedded in model-specific samplers | a typed probability provider with global, fixed, logistic, and probit-stick implementations |
| Marker variance | fixed-marker, group, learned-logistic, LV, and MAF-S implementations | one composable log-variance provider producing a positive relative \(q\) and sufficient update statistics |
| Trait sharing | global pattern state in MT kernels | a typed sharing provider; global remains the only production implementation initially |
| Cross-trait covariance | MT core state \(V_b\) | explicit covariance architecture owned by the MT core, initially global only |
| Operator | coupled to several public wrapper routes | binding-neutral BED/CSR/eigen likelihood operators consuming the same prior policies |
| Capabilities | duplicated R resolver, docs, wrappers, and tests | one declarative registry tested against dispatch, native registrations, schemas, and docs |
| Raw output | common raw core plus mechanism-specific decoration | typed architecture-result namespaces for probability, variance, sharing, covariance, annotation design, and diagnostics; supplied assumptions remain separate in `model_spec$prior` |

## Required guardrails

1. A disabled provider is a zero-overhead, zero-RNG no-op. Ordinary model
   trajectories remain exact regression controls.
2. `q` is a relative variance architecture. New composed providers use a
   declared normalization over the aligned analysis-marker universe.
3. Initialization is never interpreted as prior information. For informative
   theta, a missing initialization defaults to the prior mean, but both remain
   separately recorded and an explicit initialization overrides that default.
4. Prior component probabilities, posterior allocation probabilities, PIP,
   and posterior beta summaries remain distinct fields.
5. Mathematically composable does not mean scientifically qualified. In
   particular, probit-stick probability plus LV variance remains unavailable
   until a dedicated identifiability and mixing phase passes.
