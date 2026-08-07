# sblr 0.2.0

* Added development-only, block-SBayesRC mixing references without changing
  defaults: an exact latent-probit PX/sandwich update and an exact conditional
  particle-Gibbs block allocation/effect reference. The PX transition passed a
  tiny exact posterior oracle and improved short-run occupancy agreement, but
  did not resolve joint alpha/occupancy convergence. Conditional SMC passed
  exact tiny-block enumeration and remained diverse at 100 and 500 markers;
  its block-local active-count jumps were too small to establish a solution to
  the global alpha regime problem, so no particle kernel was added to the
  production sampler.

* Completed an exact-method design audit for a coordinated BayesRC/SBayesRC
  alpha--allocation transition. A beta-inclusive subset move passed independent
  tiny-posterior and detailed-balance checks, but fixed global alpha moves lost
  acceptance extensively with marker count while acceptance-preserving scaling
  made alpha movement vanish. No production transition, interface, default, or
  sampler RNG path was changed; the result is retained as a development
  reference for future collapsed or transport-based method work.

* Added an official-compatible block-specific residual-variance policy for
  retained block-eigen SBayesR/SBayesRC. `residual_policy = "gctb_block"` is
  now the retained-route default, with the pinned `allMixVe`, `mixVe`,
  `samVe`, and `fixVe` modes. `fixed_block` provides a transparent fixed-Ve
  sensitivity route. The prior global projected update remains available as
  `global_projected_legacy` for exact historical reproduction.
* Under block policies, `ves` is the mean block residual variance and the
  explicit summary heritability trace is `sum(Vg_b) / Vy`; neither is claimed
  to be a BED-equivalent individual-level residual quantity. Block-level
  posterior means, final values, resampling counts, and official minimum-ratio
  resets are retained as compact diagnostics.

* Added RNG-neutral compact allocation histories for BayesR/BayesRC extended
  diagnostics: component counts, realized active count, and sequential-stick
  eligible/continue/stop counts. Packed BED, CSR, and retained block-eigen
  scalar routes share the public draw-by-chain contract, avoiding full
  marker-by-iteration component storage.

* Added development-only scalar BayesRC/SBayesRC controls for auditing fixed
  repetitions of the existing allocation/model and annotation-hierarchy Gibbs
  blocks. Defaults remain one update each and preserve the historical RNG
  trajectory. Study 06 development evidence found that extra hierarchy updates
  improve conditional alpha exploration but no tested schedule resolves joint
  occupancy mixing; no package default changed.

* Repaired selected-marker component traces for packed-BED fits that use a
  fitted `cls` marker subset. Trace identifiers are now resolved against the
  fitted marker order, and native BED backends reject out-of-range trace
  indices before entering parallel chain execution. Trace capture remains RNG
  neutral and does not change sampler transitions.

* Corrected the shared BayesRC/SBayesRC annotation update for a non-intercept
  design column that is identically zero on a non-empty stick-eligible subset.
  The coefficient is now drawn from its exact hierarchical-prior conditional
  instead of retaining a stale value. A focused mixing review found that the
  remaining Study 06 development bottleneck is component-allocation movement,
  not evidence supporting a blocked or non-centred alpha update.

* Derived and validated an exact two-marker BayesRC allocation conditional.
  Study 06 development runs found the fixed-pair move too sparse and ineffective
  to retain in production; a larger collapsed allocation move requires a
  separate design review. Public samplers and interfaces are unchanged.

* BayesRC/SBayesRC now use a proper normal prior for every probit-stick
  intercept by default. Prior means reproduce the resolved global mixture and
  the default SD is 1. `alpha_init` remains an initialization only. The former
  flat prior is available only through explicit `intercept_flat = TRUE` and
  fails for empty or completely separated sticks.
* Empty annotation sticks now receive posterior-valid prior-only Gibbs updates
  instead of retaining stale coefficients. BED, CSR, retained block-eigen, and
  multitrait annotation routes share the same resolved intercept-prior kernel.

- Aligned prior variance calibration across scalar and multivariate BayesC,
  BayesR, and BayesRC routes. The requested `h2` now targets the initial/prior
  expected genetic variance under resolved component probabilities, component
  multipliers, marker MAF-S scales, annotation/group variance multipliers, and
  MT trait-pattern probabilities.
- Separated initial mixture weights from Beta/Dirichlet prior-mean weights and
  added transparent prior-calibration metadata to fitted-model inputs.
- Kept explicit variance-prior overrides authoritative and added a validation
  error for unsafe automatic MT full-covariance conversion.
- Corrected scalar BayesR prior variance initialization to account for
  component variance multipliers, so `h2` has the same prior expected-genetic-
  variance interpretation as BayesC.
- Added the retained low-rank scalar block-eigen LD operator for SBayesC,
  SBayesR, and SBayesRC and made it the default `stblr_block_eigen()`
  representation.
- Retained historical reconstructed-dense behavior explicitly through
  `representation = "dense_reconstructed"`.
- Made cumulative positive-eigenvalue-mass retention the low-rank policy,
  with `eigen_prop = 0.995` by default.
- Kept retained-factor sampling single-trait only; MT block-eigen sampling
  remains reconstructed dense.
- The retained low-rank operator follows the GCTB/SBayesRC eigenspace
  likelihood strategy in sblr cross-product units. Its default residual policy
  is now the official-compatible block-specific contract; the former global
  projected contract is explicit legacy behavior.
- Classified scalar LD operators in fit metadata: packed BED is the supplied-
  genotype reference, complete CSR construction is the summary-statistics
  reference, thresholded/windowed CSR is explicitly approximate, and retained
  block eigen is the scalable projected SBayesRC route.
- Added opt-in, failure-only SBayesRC residual-scale diagnostics. Invalid
  scales still fail without clamping, absolute-value substitution, stale-
  variance fallback, or an operator switch.
- Corrected extreme-tail truncated-normal draws used by BayesRC annotation
  updates and connected the public low-rank residual-rebuild interval to the
  scalar SBayesRC retained-eigen kernel.

# sblr 0.1.0

- Stabilized the seven canonical STBLR and MTBLR fitting interfaces and the
  individual-versus-summary model semantics for BayesC, BayesR, and BayesRC.
- Aligned CSR sparse-LD, block-eigen, and packed-BED operators under common
  output, seed, multichain, memory, and warning contracts.
- Added deterministic logical-chain execution, core and extended convergence
  diagnostics, and explicit selected-marker diagnostics.
- Added annotation-informed ST and MT models and independent fixed or sampled
  `maf_effect_s` support where scientifically implemented.
- Consolidated the framework around permanent scientific owners and removed
  obsolete compatibility, migration, and phase-oriented scaffolding.
- Repaired and expanded the package help, Notes, Methods, workflows, and
  GitHub Pages documentation.
