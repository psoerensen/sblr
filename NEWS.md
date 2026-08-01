# sblr 0.2.0

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
