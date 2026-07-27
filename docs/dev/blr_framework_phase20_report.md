# Phase 20 report: MTBLR BayesRC and SBayesRC

## 1. Executive summary

Phase 20 adds one annotation-informed extension of the Phase 19 MT pattern-by-component model to packed BED, CSR, and block-eigen operators. Annotation-dependent component probabilities are factorized from global conditional trait-pattern probabilities.

## 2. Repository baseline

Work started from `4c8cd09608efdd6f1b148a835c6d8ae5e34017fc` on `master` with a clean tree. R 4.4.1 and Rtools GCC/G++ 13.2 were used on Windows with 12 reported processors. Baseline source tests and built-package checking passed; the check had zero errors, zero warnings, and three established notes (long test path, installed library size, and legacy `std::cout`).

## 3. Existing Phase 18/19 contracts reused

The canonical fit layout, semantic metadata, state descriptor, BayesR conditional kernel, base-B update, chain dispatch, convergence engine, variance traces, memory vocabulary, and raw schema version 1 remain in force.

## 4. BayesRC/SBayesRC model semantics

`bayesrc` is individual-level packed BED; `sbayesrc` is summary-statistics CSR or block eigen. Both use `prior_kernel = "bayesrc"`, `annotation_policy = "annotation_probit_stick"`, and model-semantics version 2.

## 5. Annotation input contract

`annotations` is required and must be a finite numeric matrix or numeric/logical/factor data frame with one or more uniquely named annotation columns.

## 6. Annotation alignment

Explicit unique marker IDs are mandatory. Exact matching reorders annotations to final selected marker order; missing, extra/unused, duplicate, ambiguous, or nonfinite inputs are reported rather than silently intersected.

## 7. Annotation preprocessing

The established ST preprocessing owns intercept detection/addition, continuous standardization, binary centering policy, factor expansion, constant-column failure, and immutable processed metadata. It runs once before chain dispatch.

## 8. Component probit stick-breaking

The established ordered probit stick convention produces component probabilities in Phase 19 component order, with component zero null. `pi_floor` is applied followed by row renormalization.

## 9. Unique null state

The Phase 19 state descriptor retains exactly one null state; active states remain non-null pattern order crossed with positive-component order.

## 10. Factorized state prior

For marker `j`, null prior mass is `theta[j,0]`; active state `(p,k)` has mass `theta[j,k] * omega[p]`. Both probability vectors normalize, so the full state prior normalizes without repeated nulls.

## 11. Global trait-pattern probabilities

`omega` is global and conditional on non-null status. `updatePi` updates it from active-marker pattern counts with a positive Dirichlet prior; annotations never create pattern-specific coefficient models.

## 12. Annotation coefficient prior

One coefficient matrix is shared across traits and controls the common component state. Intercept-flat and inverse-chi-square coefficient-variance conventions are reused from ST SBayesRC.

## 13. Annotation coefficient update

Reached-stick probit augmentation, truncated-normal latent draws, Gaussian coordinate updates, and coefficient-variance updates are chain private and scheduled by `alpha_update_every`. The established coordinate update requires no matrix jitter; jitter diagnostics are therefore zero.

## 14. Effect prior

Conditional on state, the Phase 19 component-scaled base covariance and optional independent MAF-S multiplier are unchanged. Annotations affect state probabilities only.

## 15. Base-B update

The existing update removes known component and optional MAF-S scaling. No component-, pattern-, or annotation-specific B matrices were introduced.

## 16. Independent selection-S policy

`selection_s = NULL` disables MAF scaling, a finite scalar enables fixed scaling, and `-1` is the explicit unit-scale reduction. Sampled MT S remains rejected.

## 17. MAF annotation overlap

Explicit MAF/heterozygosity annotation names combined with `selection_s` set overlap metadata and emit one main-thread advisory. Neither input is removed or disabled.

## 18. Marker update

Marker-specific log prior weights feed the existing Phase 19 state likelihood/prior kernel and log-sum-exp categorical draw. Null effects remain exactly zero.

## 19. Initialization

Phase 19 beta/effective-effect/pattern/component initialization remains. `alpha_init` and positive `sigmaSqAlpha_init` are validated against processed annotation and stick dimensions. Fixed alpha remains unchanged.

## 20. Sampler update order

Summary routes retain their validated Phase 19 B-before-sweep structure; BED retains its Phase 19 marker sweep followed by probability/B/G/E updates. Within each representation, annotation pattern and coefficient updates occur after the completed marker sweep and before completed-iteration variance recording.

## 21. Public signatures

All three MT interfaces expose the same ordered annotation-control block after component controls and before selection-S controls. BayesC/BayesR reject nondefault annotation controls.

## 22. Raw schema

`mtblr_raw` remains version 1 and model semantics remains version 2. A model-specific `annotations` namespace is required only for BayesRC/SBayesRC and is validated additively.

## 23. Public output

All Phase 18/19 marker, covariance, variance-trace, convergence, chain, and memory fields remain. BayesRC-specific results live under `fit$model_parameters$annotations`.

## 24. Prior component probabilities

`prior_component_probabilities` is marker by component, includes null, is row normalized, and is the posterior mean annotation-driven prior probability—not a marker-state posterior.

## 25. Posterior component probabilities

`component_probabilities` retains its Phase 19 posterior marker-state meaning and is separately row normalized.

## 26. Pattern probabilities

Global conditional non-null probabilities use `pattern_pi_final`, `pattern_pi_mean`, and `pattern_pi_trace`. Their names avoid redefining Phase 19 joint-state `pi_*` semantics.

## 27. Core variance traces

`vbs`, `vgs`, `ves`, `vle`, and `vld` remain genuine per-chain post-burn-capable iteration traces, with `vld = vgs - vle`.

## 28. Convergence

The single Phase 18 scalar engine remains unchanged and diagnoses only the five core trait-level quantities. Annotation, probability, component, selection-S, off-diagonal covariance, and marker diagnostics remain deferred.

## 29. MT CSR implementation

`mtblr_csr(method = "sbayesrc")` shares one prepared CSR representation and one processed annotation matrix across deterministic joint chains.

## 30. MT block-eigen implementation

`mtblr_block_eigen(method = "sbayesrc")` shares one reconstructed owner representation, transformed scores, state descriptor, and processed annotations. Exact unfiltered reconstruction reduces to CSR.

## 31. MT packed-BED implementation

`mtblr_bed(method = "bayesrc")` shares one packed genotype preparation and processed annotation matrix and supports full or diagonal residual covariance.

## 32. Multichain execution

Every chain owns RNG, marker/component/pattern state, annotation coefficients/latent workspace, covariance state, and accumulators. Seeds and result order are worker independent; no R/Rcpp work occurs in OpenMP workers.

## 33. Memory

Analytical estimates separately report shared annotations, private per-worker annotation state, per-chain annotation results, and formatted prior-probability output. They are pre-execution estimates, not RSS.

## 34. BayesR reduction

Fixed intercept-only coefficients and fixed probabilities reproduce Phase 19 SBayesR exactly when state priors match and annotation updates are disabled.

## 35. One-trait ST/MT reduction

The shared probit-stick oracle and fixed-prior contract match the established ST convention. Posterior equality is asserted only where scalar/joint covariance priors are matched.

## 36. CSR/block-eigen reduction

The permanent small exact-reconstruction owner compares marker summaries, component/prior probabilities, pattern and annotation summaries, and all five variance traces.

## 37. Individual/summary prior reduction

Packed BED and summary routes use different data levels but the same prior kernel, preprocessing contract, state/component order, and annotation probability oracle; likelihood-dependent posterior equality is not claimed.

## 38. Selection-S reduction

The independent fixed `selection_s = -1` policy is the numerical unit-scale reduction for all three operators while preserving explicit user-intent metadata.

## 39. Annotation alignment/invariance

Marker-row permutation followed by ID realignment is exact. Column, intercept, binary/continuous, and standardization contracts have deterministic validation owners.

## 40. Update-control reductions

Owners cover fixed/updated alpha, fixed/updated pattern probability, fixed covariance controls, scheduled alpha updates, both BED E policies, and one/multiple component cases.

## 41. Reproducibility

Repeated serial fits, deterministic explicit seeds, worker-count changes, compact-chain retention, and convergence-trace retention preserve numerical results after timing normalization.

## 42. Architecture audit

The Phase 20 audit checks 25 guards covering factorization, shared ownership, public activation, semantics, raw schema, convergence, memory, and reductions.

## 43. Mutation sensitivity

Thirty-nine guards cover the required null/order/normalization/factorization, alignment, preparation, RNG, update-control, output, convergence, schema, and Phase 18/19 protection mutations.

## 44. Benchmark

The maintained benchmark records preprocessing, total runtime, shared/private/result/output bytes, fit size, and convergence size for all three operators. Values are regression signals, not production-runtime or speedup claims.

## 45. Public API and documentation

The public method matrix is packed BED `bayesc/bayesr/bayesrc` and summary CSR/block eigen `sbayesc/sbayesr/sbayesrc`. Documentation distinguishes prior from posterior component probabilities.

## 46. Focused tests

The three Phase 20 scientific owners passed 101 focused expectations. Phase 18/19 convergence, model-semantics, public-contract, operator-reduction, raw/formatter, and historical Phase 17 protection owners also passed after intentional native-formal migrations.

## 47. Fast tier

The exact fast CI filter passed with zero failures, errors, or test warnings and one established opt-in peak-RSS skip. It owns the probability oracle, fixed-alpha reduction, exact CSR/block reduction, method/alignment/raw/convergence protection, and package check.

## 48. Full source suite

The final full source suite passed with zero failures, zero errors, and zero test warnings. The two established opt-in skips were extended fresh-process coverage and peak-RSS measurement; ordinary fresh-process owners ran successfully outside the sandbox.

## 49. Installed tests

Built-tarball installed tests passed.

## 50. Package check

R CMD check completed with 0 ERRORs, 0 WARNINGs, and 3 established NOTEs: the long Phase 17C fixture path, installed size (5.5 MB; libs 4.4 MB), and legacy ST native objects using `std::cout`.

## 51. Diff and artifact hygiene

`git diff --check` passes. Generated Rcpp wrappers and registration changed only for the required additive native arguments; NAMESPACE and DESCRIPTION are unchanged. Roxygen changed exactly the three MT public Rd pages. Task-generated objects, DLL, tarball, and check directories were removed; no fixture was regenerated.

## 52. Deviations and blockers

The annotation update reuses the validated coordinate-wise ST probit augmentation and therefore has no coefficient matrix jitter path; zero jitter diagnostics are explicit. Sampled MT `selection_s` and annotation-dependent trait-pattern probabilities remain unsupported by design.

## 53. Recommended next phase

> implement extended STBLR and MTBLR diagnostics for low-dimensional covariance, probability, annotation, and selection-S parameters, followed by explicitly selected marker traces under strict memory guards.

## 54. Readiness marker

PHASE 20 COMPLETE — MTBLR BAYESRC/SBAYESRC ANNOTATION PRIORS ACTIVE
