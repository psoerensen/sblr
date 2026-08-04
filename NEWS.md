# sblr 0.2.0

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
  likelihood strategy, represented in sblr cross-product units with a global
  projected residual-variance contract. GCTB's block-specific variance
  procedure is not reproduced.
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
