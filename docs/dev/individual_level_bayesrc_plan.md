# Individual-Level BayesRC Plan

This developer note tracks the staged internal work toward a possible
individual-level BayesRC backend. It does not define or expose a public model.

## Sequence 1 Status

- likelihood-independent probit stick-breaking annotation-prior utilities were
  extracted into `src/st_bayesrc_annotation_prior.h`
- the existing CSR SBayesRC sampler remains the only consumer
- no individual-level BED BayesRC backend was added
- there are no public R API changes and no new statistical model
- the existing SBayesRC component ordering, annotation updates, RNG draw order,
  and formatted output are preserved

## Sequence 2 Status

- added the internal native full-sweep BED BayesRC backend
  `stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc`
- the backend uses the packed-BED BayesR likelihood with the shared probit
  stick-breaking annotation prior
- the public `stblr_bed()` API is not extended in this sequence
- every MCMC iteration visits every selected marker; adaptive null scheduling
  is intentionally absent
- annotation rows are supplied in the final selected BED marker order defined
  by the concatenated `cls` marker map; the low-level caller is responsible for
  applying that same index to `A`
- `selection_s`, LD-swap, annotation-dependent effect variances, and
  multivariate covariance sampling remain deferred

## Sequence 2A Status

- extracted shared packed-BED BayesR-family utilities into
  `src/st_bed_bayesr_common.h`
- ordinary BED BayesR and internal BED BayesRC now include the shared header
- removed the macro-renamed inclusion of the complete BayesR `.cpp` file from
  the BayesRC translation unit
- no public API changes were made
- ordinary BayesR adaptive scheduling remains unchanged
- internal BayesRC remains structurally full-sweep with adaptive skipping off
- raw schemas and formatted outputs remain unchanged

## Sequence 2B Status

- added fixed-prior reduction coverage against ordinary fixed-pi BED BayesR
- added learned annotation-prior and directional enrichment smoke coverage
- added multiple-trait, multiple-chain, retained-chain, repeated-seed, and
  OpenMP thread-count reproducibility checks
- added component probability, `dm`, `dm_component_mean`, `ncomp`, and `pis`
  identity checks
- added annotation output dimension and final marker-prior checks
- added variance trace, CPO, thinning, optional `wy`/residual, native input,
  full-sweep diagnostic, and raw schema validation
- retained regression coverage for ordinary BED BayesR and CSR SBayesRC through
  their existing focused test groups
- annotation rows remain an internal caller responsibility: `A` must already
  follow the final selected BED marker order represented by `cls`; a future
  public wrapper must perform marker-ID alignment before calling native code
- the backend remains internal; public `stblr_bed(method = "bayesrc")` routing
  is not implemented

### Compiled Sequence 2B correction results

- compiled execution confirmed that an `Rcpp::List` initialized from
  `R_NilValue` materialized as an empty list; the optional holder is now an
  `Rcpp::RObject`, so `chains` is actual R `NULL` when retention is disabled
  and remains a trait-by-chain nested list when enabled
- compiled repeated-call testing confirmed the persistent `thread_local`
  `std::normal_distribution` cache as the cause of both repeated-seed and
  thread-count differences; ordinary BED BayesR and BED BayesRC now construct
  uniform and normal distribution state beside each chain-local `mt19937` and
  pass it by reference to marker updates
- chain seeds, marker order, RNG engine, random-update sites, and within-chain
  draw order are unchanged; identical-seed calls, one-core/two-core calls in
  both orders, and reference calls separated by an unrelated fit are identical,
  while a different seed changes stochastic output
- the fixed-prior intercept-only BayesRC run still reduces to fixed-pi,
  full-sweep BayesR at the existing `1e-12` tolerance for marker means,
  inclusion/component probabilities, terminal effects/states, variance traces,
  and CPO diagnostics
- retained-chain reconstruction passes for `bm`, `dm`, component probabilities,
  annotation coefficient means, and `sigmaSqAlpha` means; retained Armadillo row
  vectors remain `1 x m`, while matrix-valued component and annotation fields
  retain their dimensions
- component identities, marker-prior column means, marker-average non-null
  prior probability, `vld = vgs - vle`, CPO normalization, and retained-sample
  counts all pass compiled tests; learned annotation coefficients and variances
  are finite and annotation variances are positive
- the directional fixture now uses 80 individuals and 16 independently
  generated polymorphic markers, with six binary-annotated markers, four
  annotated causal markers, ten unannotated null markers, fixed genotype and
  phenotype seeds, and residual noise SD 0.45; both mean annotated non-null
  prior probability and the posterior mean first-stick enrichment coefficient
  are directionally positive
- the compiled focused file passes 267 expectations with no failures, warnings,
  or skips; the unlimited full package suite passes 3,091 expectations with no
  failures, warnings, or skips, including ordinary BED BayesR and CSR SBayesRC
  regression groups
- `R CMD build .` succeeds; `R CMD check --no-manual` installs, loads, compiles,
  and runs examples successfully but is not clean at this research checkpoint:
  tarball tests that read repository-only `src/` and `docs/dev/` paths fail
  because those paths are absent from the installed test layout, and existing
  `make_credible_sets.Rd` syntax plus namespace/global-function notes remain

## Sequence 3 Status

- added the public `stblr_bed(method = "bayesrc", annotation = ...)` route
- annotation rows are matched by marker ID to the exact final BED order formed
  by concatenating `Glist$rsids[[chr]][cls]`; extra annotation rows are recorded
  and dropped, while missing or duplicated IDs fail clearly
- numeric, logical, integer, and factor data-frame annotations are converted to
  a numeric design, with factors expanded to indicators; intercept addition,
  continuous-column standardization, and optional binary centering follow the
  established CSR SBayesRC preparation controls
- baseline `pi` initializes the probit stick-breaking intercept coefficients;
  non-intercept coefficients default to zero and `sigmaSqAlpha` defaults to one
  per stick
- the validated `bed_bayesrc` `stblr_raw_v1` output is formatted through the
  shared raw-schema formatter and existing SBayesRC annotation aliases
- public component, annotation, chain, variance, and CPO outputs are exposed;
  BayesRC remains unscheduled and full-sweep
- the public fixed-alpha, intercept-only reduction matches fixed-pi BED BayesR
  to `1e-12` for `bm`, `dm`, and component probabilities; shuffled annotation
  rows, chromosome/column subsets, public-to-native equivalence, multiple traits,
  retained chains, variance identities, and CPO outputs pass their compiled tests
- the unlimited Sequence 3 package suite passes 3,181 expectations with no
  failures, warnings, or skips; the focused public BED interface file passes
  183 expectations and the BED block-eigen regression file passes 193
- `R CMD build .` succeeds; `R CMD check --no-manual` installs, compiles, loads,
  checks examples, and checks compiled code successfully, but remains non-clean
  because installed-tarball tests read repository-only `src/` and `docs/dev/`
  paths, alongside the pre-existing `make_credible_sets.Rd` and namespace notes
- selection-S, LD-swap, adaptive scheduling, BED block-eigen fitting,
  annotation-dependent effect-size variances, and covariance adjustment remain
  disabled
