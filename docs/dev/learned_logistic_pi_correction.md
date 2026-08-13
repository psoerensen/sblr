# Learned-logistic global-pi correction

## Status and scope

This record documents the phase-3A correction of the single-trait CSR BayesC
learned-logistic probability provider. The correction started from branch
`master` at commit `fea4c4389be14c02317d234ba652d42e3f68a42b` (`fea4c43
Phase-2 audit report`). The index and working tree were clean at the start.

The change is intentionally limited to the learned-logistic provider used by
`stblr_csr_annot(annotation_model = "learned_logistic")`. It does not alter
ordinary, fixed-probability, or group-probability CSR BayesC; BayesR or
SBayesRC; log-variance providers; BED or block-eigen fitters; multitrait
models; PIP definitions; or formatted output schemas.

## Authority and evidence read

The work followed the local root `AGENTS.md` (no applicable nested
`AGENTS.md` was present) and the authority classification in
`docs/dev/README.md`. Statistical authority came from the approved BayesC and
annotation-prior Methods specifications, principally
`docs/methods/model_theory.qmd`, `docs/methods/annotation_priors.qmd`, and
`docs/methods/annotation_informed_models.qmd`, together with the explicit
phase-3A decision that the unclipped logistic model is authoritative.

Implementation evidence included
`docs/dev/methods_source_traceability_audit.md`, the maintained BLR
development, architecture, backend, model, output, convergence, and test
contracts, and the maintained annotation-prior architecture audit, matrix,
and implementation plan. The full executable path was traced through
`stblr_csr_annot()`, `stblr_csr_learn_annot()`, the generated Rcpp boundary,
`stblr_cpg_omp_csr_annot()`, and the learned-annotation core and policy types,
including their permanent tests and raw-to-fit path.

## Defect and clipping contradiction

For centered annotation offsets $o_j$, the approved model is

$$
\operatorname{logit}(\pi_j)
=
\operatorname{logit}(\pi)+o_j.
$$

The previous implementation applied hard bounds after the inverse-logit,
including a default upper bound of $0.5$. Marker-state updates and the
annotation-coefficient Metropolis--Hastings target therefore used a clipped
Bernoulli model. In addition, the global active probability was updated with
the ordinary conjugate count draw

$$
\pi\mid d
\mathrel{\overset{\mathrm{old}}{\sim}}
\operatorname{Beta}
\left(a_\pi+\sum_jd_j,\ b_\pi+M-\sum_jd_j\right)
$$

even when offsets were nonzero. That draw is not the full conditional of the
logistic-offset model. Centering the offsets does not restore conjugacy.

The hard cap also made the nominal zero-offset route discontinuous from
ordinary BayesC: once a sampled global probability exceeded the cap, the
marker likelihood no longer used that global probability. The phase-3A
decision therefore removes the hard bounds from the learned provider's
statistical target rather than preserving a special clipped target.

## Corrected statistical target

Marker-state updates, annotation-coefficient updates, and global-$\pi$
updates now use the same probabilities,

$$
\pi_j
=
\operatorname{logit}^{-1}
\left\{\operatorname{logit}(\pi)+o_j\right\},
$$

without a scientifically material lower or upper bound. With
$\pi\sim\operatorname{Beta}(a_\pi,b_\pi)$, the corrected conditional on the
probability scale is

$$
p(\pi\mid d,o)
\propto
\pi^{a_\pi-1}(1-\pi)^{b_\pi-1}
\prod_{j=1}^M
\pi_j(\pi)^{d_j}
\{1-\pi_j(\pi)\}^{1-d_j}.
$$

For $z=\operatorname{logit}(\pi)$, the Jacobian is

$$
\left|\frac{d\pi}{dz}\right|
=
\pi(1-\pi).
$$

Consequently, up to an additive constant, the implemented logit-scale
density is

$$
\log p(z\mid d,o)
=
a_\pi\log\pi
+
b_\pi\log(1-\pi)
+
\sum_{j=1}^M
\left[
d_j\log\pi_j
+
(1-d_j)\log(1-\pi_j)
\right].
$$

The coefficients are $a_\pi$ and $b_\pi$, not $a_\pi-1$ and $b_\pi-1$,
because the transformation contributes one additional power of both $\pi$
and $1-\pi$.

If every effective offset is exactly zero, the same continuous model reduces
analytically to

$$
\pi\mid d
\sim
\operatorname{Beta}
\left(a_\pi+\sum_jd_j,\ b_\pi+M-\sum_jd_j\right).
$$

That exact case retains the existing conjugate gamma-ratio draw and its draw
order. It is a reduction of the corrected target, not a separate target.

## Implementation and numerical safeguards

`src/st_cpg_omp_csr_annot.cpp` now owns stable softplus, log-sigmoid,
inverse-logit, and strict-logit helpers for this provider. Likelihood terms are
evaluated as stable log probabilities rather than as products or logarithms
of rounded probabilities. Extreme finite predictors therefore retain finite
log-likelihood contributions. Open-interval probability values used for
stored marker probabilities use only representability protection at zero and
one; those values are not used to truncate the log posterior.

For nonzero offsets, `samplePi_ST_annot()` uses a univariate slice sampler on
$z$, with finite stepping-out and shrinkage iteration limits. The global
probability retains its natural domain, $\pi\in(0,1)$, and marker probabilities
have no scientifically material cap. The logical chain's existing
`std::mt19937` engine is used with a width of 1, at most 100 stepping-out
expansions, and at most 1,000 shrinkage iterations. These finite limits are
internal termination safeguards; exhaustion raises an informative error.
Prior shapes, states, offsets, the current probability, and all returned
values are validated. All-null and all-active states and prior shapes below or
above one are supported.

The learned core maintains centered offsets, marker logits, and representable
marker probabilities together. It recomputes all three after either $\eta$ or
$\pi$ changes. Marker-state and LD-swap prior odds use logits directly;
`logpost_eta_pi()` uses the same stable, unclipped log probabilities; and the
global $\pi$ update uses the same offsets. `updatePi = FALSE` performs no
global probability draw and retains the supplied global probability.

The exact-zero branch preserves the established conjugate RNG path. The
nonzero branch intentionally consumes slice-sampler uniforms instead of the
obsolete conjugate gamma draws. No R RNG is called from native worker code,
and no shared mutable RNG state was introduced.

A tiny private-copy comparison ran the same zero-offset learned route at the
clean `fea4c438` checkpoint and after the correction. With fixed logical-chain
seed and fixed states, `pi_trace` and `pi_final` were bitwise identical. Thus
the retained branch preserves both the conjugate draw order and its established
trajectory for ordinary finite draws.

## Treatment of `pi_min` and `pi_max`

Hard probability bounds are not controls of the approved learned-logistic
model. The learned provider therefore no longer has native `pi_min` or
`pi_max` parameters. At the public `stblr_csr_annot()` boundary, supplying
either obsolete name, any prefix of those names, or an abbreviated prefix of
`pi_marker` for `annotation_model = "learned_logistic"` now fails clearly
before internal dispatch. The exact name `pi_marker` remains supported as the
initial global probability control, and unrelated unknown arguments retain
R's ordinary unused-argument error. This prevents R partial argument matching
from silently treating an obsolete or ambiguous bound name as `pi_marker`.
Every argument forwarded through the learned-provider `...` path must also
have a unique, nonempty name: unnamed, empty-name, `NA`-name, mixed
named/unnamed, and duplicate-name inputs are rejected before `do.call()` can
bind them positionally or ambiguously.
This is the smallest explicit development-stage API correction: the public
function's formal arguments and fit schema are unchanged because these
controls arrive only through `...`.

Fixed marker-probability calibration continues to use its existing bounds.
The shared R calibration helper now selects clipping explicitly so separating
the learned route does not change fixed-provider semantics.

## Source ownership and files changed

- `R/stblr-csr-annot.R`: rejects obsolete learned-provider bounds, ambiguous
  prefixes, and unnamed, empty-name, `NA`-name, or duplicate forwarded
  arguments before dispatch, while retaining exact `pi_marker`.
- `docs/methods/annotation_priors.qmd`: clarifies that the logistic range is
  natural and that the learned provider has no additional material hard cap.
- `R/stblr-csr-learn-annot.R`: removes learned bound forwarding, validates
  finite initial coefficients, and requests unclipped learned calibration.
- `R/annotation-helpers.R`: separates fixed-provider clipping from learned
  calibration.
- `src/st_cpg_omp_csr_annot.cpp`: stable logistic calculations, coherent
  marker and $\eta$ targets, corrected slice/Beta global update, validation,
  narrow native test entry points, and internal sampler-region parallel
  diagnostics.
- `src/blr_csr_annotation_bayesc_types.h`,
  `src/blr_csr_learned_annotation_bayesc_types.h`, and
  `src/blr_csr_learned_annotation_bayesc_core_impl.h`: remove cap fields and
  maintain offsets, logits, probabilities, and update ordering coherently.
- `R/RcppExports.R` and `src/RcppExports.cpp`: generated signatures and the
  internal analytical-test entry points. They were generated and compared in
  a private copy.
- `tests/testthat/test-learned-logistic-pi-conditional.R`: independent density,
  quadrature, sampler, boundary, and numerical-stability tests.
- `tests/testthat/test-stblr-annotation-interface.R` and
  `tests/testthat/test-stblr-annotation-bayesc-chains.R`: public contract,
  reduction, isolation, and logical-chain reproducibility tests.
- `docs/dev/annotation_prior_architecture_audit.md`,
  `docs/dev/annotation_prior_architecture_matrix.md`, and
  `docs/dev/blr_model_contracts.md`: corrected maintained implementation
  contracts.

No Methods formula, posterior target, or output schema was changed. A focused
post-verification pass clarified one adjacent Methods sentence about the
natural logistic probability range.

## Independent numerical evidence

The reference tests calculate normalized one-dimensional densities by dense
trapezoidal quadrature, independently of the sampler.

For zero offsets with the test's fixed states and prior, the quadrature and
analytical Beta conditional agree in density and give mean `0.4625` and
variance `0.02762153`. The retained conjugate branch reproduces the same
moments and is bitwise deterministic for a fixed seed.

For the centered, unequal offsets

$$
(-4,-2,-1,1,2,4),
$$

which sum to zero, quadrature gives mean `0.257155411`, variance
`0.028646626`, and quartiles `0.123710756`, `0.221965311`, and `0.358760037`.
A deterministic corrected slice run gives mean `0.257926564`, variance
`0.028523170`, and quartiles `0.124606392`, `0.223354173`, and `0.360176716`,
within prespecified Monte Carlo tolerances. The obsolete Beta-count update
would instead have mean `0.3375` and variance `0.02484375`, a material
disagreement.

For the separate unequal, uncentered offsets `(-3,-1,0.5,2,4)`, the same
independent quadrature gives mean `0.4027448`, variance `0.03961301`, and
quartiles `0.2451775`, `0.3860373`, and `0.5458707`. The corresponding old
Beta-count mean is `0.4722222`; the permanent deterministic slice test checks
the corrected mean, variance, and quartiles against this numerical reference.

A separate test-side density calculation shows that the former $0.5$ cap
changes the conditional even at zero offsets. For its example, the clipped
density has mean `0.550129` and variance `0.0467716`, versus the correct Beta
mean `0.4625` and variance `0.02762153`. The old clipped construction is used
only as counterfactual evidence.

The native $\eta$ log posterior is also compared with an independent R
Bernoulli likelihood. Predictors from $-1000$ to $1000$ verify stable finite
log probabilities and valid open-interval representable probabilities.

## Permanent tests and isolation

Permanent tests cover zero-offset density and Beta moments, seeded conjugate
draws, centered non-Beta quadrature, slice mean/variance/quantiles, the old cap
counterexample, the independent $\eta$ target, all-null and all-active states,
extreme finite predictors, prior shapes below and above one, invalid shapes,
invalid states, nonfinite offsets and coefficients, probability integrity,
`updatePi = FALSE`, same-seed reproduction, supported serial/parallel logical
chain equality, and rejection of obsolete or ambiguous public bound names.
The serial/parallel regression uses two traits, verifies that the OpenMP
runtime can provide two workers before running, and then directly checks the
learned sampler region's actual team size and zero-based trait worker IDs. It
requires two distinct sampler worker IDs in each chain while exercising
nonzero offsets and the nonconjugate global $\pi$ update. A separate
public-interface regression confirms that exact `pi_marker` remains valid,
unsafe unnamed or duplicate forwarding is rejected, and fixed-provider bounds
retain their calibration behavior.

Focused existing tests for ordinary, fixed, and group CSR BayesC, learned
annotation, raw fields, prior calibration, LD swap, variance output, dispatch,
and unified reproducibility passed without fixture changes. No frozen expected
hash was updated.

## Validation commands and outcomes

All compilation and mutation-capable checks ran in private copies made from
the starting checkpoint with only the recorded phase-3A changes overlaid.

- `Rcpp::compileAttributes()` in a private copy: passed; generated R and C++
  interfaces matched the working changes exactly and `NAMESPACE` did not
  change.
- R parsing of changed R sources and tests: passed.
- New analytical/quadrature test file: passed.
- Focused learned-annotation, CSR BayesC, schema, LD-swap, output, dispatch,
  RNG, serial/parallel, prior-calibration, and isolation tests: passed.
- Maintained fast scientific/contract suite: passed, with one existing
  synthetic `diag(V)` warning.
- Full ordinary `testthat` source suite: passed; seven expected short-chain or
  synthetic warnings were emitted and one explicitly opt-in extended
  fresh-process test was skipped by its environment gate.
- After the final zero-offset arithmetic and test-ownership refinements, the
  installed corrected package passed the complete analytical/quadrature file
  (48 expectations) and the new unmasked learned-chain reproducibility,
  `updatePi = FALSE`, obsolete-bound, and nonfinite-input contract checks.
- `Rscript tools/audit/blr_generated_interfaces_audit.R .`: passed, with 74
  wrappers and 74 registrations.
- `Rscript tools/audit/blr_architecture_audit.R .`: all 19 checks passed.
- `Rscript tools/audit/blr_documentation_audit.R .`: the Methods and changed
  contracts introduced no failure; the command retained the pre-existing,
  out-of-scope link false positive
  `docs/dev/sbayesrc_s_em_phase5c.md -> docs/dev/delta, alpha`.
- `Rscript tools/check/check_package.R .` in a lean private checkpoint copy:
  source package build, installation, loading, examples, compiled-code checks,
  documentation checks, and 5,405 installed tests passed. The check ended with
  two failures in unrelated opaque retained block-eigen log-variance hashes,
  plus one warning, seven skips, and the installed-size NOTE. Running the two
  hash tests against a separately compiled, privately installed clean
  `fea4c438` checkpoint reproduced the exact same actual hashes
  (`2e62c11f16e0c94b6b0e0f4fd1041ed3` and
  `699c87615c816c7f92ab76be1fd33090`). They are therefore pre-existing
  installed-check/hash consequences, not phase-3A failures; the frozen values
  were not changed.

A focused post-verification pass then exercised the public-name and worker
contracts added after independent review:

- the complete annotation-interface test file passed 157 expectations,
  including every nonempty prefix of `pi_min` and `pi_max`, every non-exact
  prefix of `pi_marker`, exact `pi_marker`, unnamed, empty-name, `NA`-name,
  mixed named/unnamed and duplicate-name inputs, an unrelated unknown
  argument, and unchanged fixed-provider bound calibration;
- the 48-expectation analytical/quadrature file passed unchanged;
- the annotation-chain test file passed 159 expectations with a two-trait
  fixture after the OpenMP probe reported that two workers were available.
  Internal diagnostics from the learned sampler region then recorded actual
  team size two and two distinct zero-based trait worker IDs in every chain.
  Serial repetition and one-worker versus two-worker comparisons were
  identical for the contracted marker, state, probability, annotation,
  variance, chain, and convergence outputs. A stronger synthetic score fixture
  gave both traits positive derived genetic-variance diagonals, and all three
  fits completed without the previous `cov2cor` warning;
- the generated-interface audit again passed with 74 wrappers and 74
  registrations, and the architecture audit again passed all 19 checks;
- the documentation audit introduced no Methods or changed-contract failure
  and retained only the pre-existing out-of-scope
  `docs/dev/sbayesrc_s_em_phase5c.md -> docs/dev/delta, alpha` link finding.

An initial invocation of the build wrapper used a disposable copy that had
inadvertently included ignored `results/local` research artifacts and failed
while copying them because the temporary volume filled. That copy was removed
and replaced with the lean checkpoint-faithful copy above. No repository file
was involved in the storage failure.

## Consequences and remaining limitations

Existing learned-logistic results should be rerun when annotation offsets were
nonzero, when `updatePi = TRUE`, or when the old bounds affected a marker
probability. The correction can change the posterior target, annotation
coefficients, global probability, marker states, effects, and downstream
summaries. Ordinary, fixed, and group BayesC results are not affected.

The slice width and iteration limits are internal numerical controls rather
than public tuning arguments. Exhaustion fails explicitly. The provider still
uses the existing centered linear annotation predictor and random-walk
Metropolis--Hastings coefficient update; this task did not redesign those
approved features. The frozen historical phase-0 record continues to describe
the old implementation as provenance and was not rewritten. The multitrait
$V_b$ issue and every other phase-2 finding remain outside this correction.

## Git state

All phase-3A source, test, generated-interface, contract, and report changes
are left as unstaged working-tree changes. No commit, push, branch operation,
tag, pull request, publication, reset, restore, clean, stash, or unrelated
file modification was performed.
