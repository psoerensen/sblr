# Scalar block-eigen backend contract

> **Status: CURRENT CONTRACT.** This file was intentionally reintroduced and
> revised with the retained low-rank implementation after the historical
> cleanup manifest. Current source and tests use the contract identifiers
> below. See [`README.md`](README.md) for the authority hierarchy.

`stblr_block_eigen()` defaults to `representation = "low_rank"`, contract
`block_low_rank_v1`, policy `cumulative_positive_mass`, and
`eigen_prop = 0.995`. It uses the general cross-product scale documented in
`stblr_low_rank_operator_design.md`, keeps SNP effects and priors in marker
space, shares immutable factors across chains, and gives each chain its own
double reduced residual.

`representation = "dense_reconstructed"` is the unchanged historical packed
backend, contract `block_dense_reconstructed_v1`. Its absolute threshold and
ridge policies remain reference/reproducibility routes. Frozen dense numerical
regressions must request it explicitly.

The retained operator is scalar-only. The obsolete MT reconstructed
block-eigen covariance-hybrid route is no longer exported in Phase 5B;
corrected MT block-eigen likelihoods remain unavailable. The scalar route does
not discover blocks, represent cross-block LD, or select a hybrid operator
automatically.

See `stblr_low_rank_operator_design.md`, `stblr_low_rank_gctb_crosswalk.md`, and
`stblr_low_rank_performance.md` for mathematics, external scale mapping, and
validation evidence.
