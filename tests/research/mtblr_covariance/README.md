# MTBLR covariance research prototype

This directory contains a standalone, inspectable R prototype for the first
MTBLR covariance-design phase. It is research code, not a package backend or a
permanent `testthat` contract. The prototype separates latent slab effects
`beta`, realised effects `alpha`, marker states, and covariance draws. Its
shared-component BayesR transition stores the scaled component effect `theta`;
for an active all-traits component, `alpha` equals `theta`.

Files:

- `mtblr_reference_model.R`: validation, inverse-Wishart convention, marker
  conditionals, and diagnostics;
- `mtblr_exact_reference.R`: exact discrete-state enumeration for known
  covariance and an independent low-dimensional covariance grid;
- `mtblr_full_latent.R`: full latent augmentation;
- `mtblr_active_only.R`: active-only BayesC and a scaled-effect BayesR
  transition with marker-specific positive multipliers `q`;
- `mtblr_current_hybrid.R`: diagnostic reconstruction of the current MTBLR
  covariance transition;
- `compare_samplers.R`: multi-chain comparison and mixing summaries;
- `test_reference_cases.R`: deterministic identities, numerical-reference
  refinement, reproducibility, unequal-`q` BayesR enumeration, and Monte Carlo
  checks;
- `mtblr_pattern_reference.R`: four-pattern state ordering, conditional
  completion, exact configuration enumeration, and fixed regional references;
- `mtblr_pattern_samplers.R`: Cheng full-latent and null-collapsed,
  conditionally completed samplers with joint or coordinate state updates;
- `mtblr_sampled_residual.R`: independent degrees-of-freedom/scale residual
  inverse-Wishart conditional and completed-active Cheng sampler with sampled
  $\boldsymbol{\Pi}$, $V_b$, and $V_e$;
- `test_sampled_residual.R`: Phase 4b conditional, posterior-mean, residual-
  reconstruction, and sampled-chain reference checks;
- `mtblr_general_t_reference.R`: independent complete $2^T$ pattern
  enumeration, arbitrary active-subset marker conditionals, Schur completion,
  and $V_b$/$V_e$ sufficient-statistic calculations for Phase 5A tests;
- `mtblr_regional.R`: persistent regional covariance prototypes and a
  source-faithful reconstruction of the current set-loop hybrid;
- `compare_pattern_samplers.R`: learned-pattern and joint/coordinate mixing
  comparisons;
- `test_pattern_reference_cases.R`: deterministic and Monte Carlo four-pattern
  and two-region reference checks, including coordinate-graph validation and
  the exact shared-global reduction;
- `mtblr_provider_operators.R`: provider-local marker maps and transparent
  dense, CSR, and block-eigen summary-likelihood operators;
- `test_provider_operators.R`: heterogeneous-provider compatibility,
  exact-representation, retained-rank, provider-local/order invariance, and
  overlap-boundary checks.

From the repository root, run:

```powershell
Rscript tests/research/mtblr_covariance/test_reference_cases.R
Rscript tests/research/mtblr_covariance/compare_samplers.R
Rscript tests/research/mtblr_covariance/test_pattern_reference_cases.R
Rscript tests/research/mtblr_covariance/compare_pattern_samplers.R
Rscript tests/research/mtblr_covariance/test_provider_operators.R
Rscript tests/research/mtblr_covariance/test_sampled_residual.R
```

Only base R is required. Seeds are explicit. The inverse-Wishart convention is
documented in `docs/dev/mtblr_covariance_design.md`.

For shared-component BayesR, the prototype uses
$\boldsymbol{\theta}_j\mid k_j,V_b\sim
N_T(\mathbf 0,\gamma_{k_j}q_jV_b)$ and stores
$\boldsymbol{\alpha}_j=\boldsymbol{\theta}_j$ for an active all-traits
component. Consequently, its covariance statistic divides
$\boldsymbol{\theta}_j\boldsymbol{\theta}_j^\top$ by
$\gamma_{k_j}q_j$ exactly once. The default $q_j=1$ is a reduction, not the
complete scale contract. Trait-specific component assignments require a
separate derivation and are not implemented here.

Regional pattern weights and a shared-global pattern state have different
prior contracts: the former receives one Dirichlet prior vector per region;
the latter receives one explicitly declared global vector, which is never
constructed by summing replicated regional priors. Coordinate state updates
are accepted only when the declared pattern graph is connected by
one-coordinate moves; joint categorical updates remain available for valid
disconnected restricted spaces.
