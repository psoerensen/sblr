# Study 12 — SBayesR-LV annotation evaluation, Aspect 1

Status: preregistered before any new SBayesR-LV comparative result was
inspected (2026-08-09).

## Scientific question

On the immutable Study 10 informative-annotation simulation, does qualified
SBayesR-LV use annotation information relative to ordinary SBayesR, and does
its downstream inference resemble SBayesRC while its low-dimensional
annotation posterior is more stable?

The three and only primary arms are ordinary SBayesR (R), SBayesR-LV (LV),
and SBayesRC (RC). This is one simulation and not a model-ranking claim.

## Frozen implementation and data

Every arm must load `sblr` 0.2.0 from the clean sibling source at Git SHA
`2123699a9cc2e91059e7d81a745420b14eca7f6e`. An isolated study-local library
is allowed; an unrelated installed copy is not. A dirty or different sibling
tree is a hard stop.

The data are referenced by `canonical_annotation_simulation_v1`, which points
to the existing immutable Study 10 semantic checkpoint. No phenotype,
genotype, truth, annotation, summary statistic, sample split, marker order, or
block is regenerated. The registered coordinates are 4,500 training people,
500 validation people, 35,000 markers, 35 contiguous 1,000-marker blocks,
target h2 0.5, target counts 34,000/500/300/200 and realized counts
33,989/502/305/204. All checkpoint and semantic hashes in the manifest must
match before fitting.

Raw annotations are `Intercept`, `enriched_binary`, `continuous_signal`, and
`null_annotation`, aligned to the frozen marker IDs. LV receives the three
nonintercept raw biological columns: production preprocessing centers binary
columns and centers/scales continuous columns to SD one exactly once. RC
receives the established four-column Study 10 matrix without further
standardization or binary centering. Numerical theta and alpha are not common
estimands; induced architectures are compared.

## Common operator and model settings

All arms use the same frozen summary statistics, working Glist, block starts,
marker order, and retained block-eigen route. The fixed operator contract is
low-rank cumulative-positive-eigenvalue-mass retention at 0.995,
`residual_policy="gctb_block"`, `block_ve_mode="allMixVe"`,
`resam_thresh=1.1`, minimum Ve ratio 0.7, block Ve histories retained, and
residual rebuild every 100 iterations. Retained ranks are not optimized by
arm and must match after fitting.

All arms use gamma `(0, 0.01, 0.1, 1)`, the target count proportions as
starting Pi, initial h2 0.5, and common B/E semantics. LV uses the qualified
`annotation_model="log_variance"` interface with `theta_prior_sd=0.7`, zero
theta initialization, and production annotation preprocessing. RC uses the
established Study 10 configuration: `sigmaSqAlpha_a=2`,
`sigmaSqAlpha_b=2`, proper normal intercept priors centered on the initial
reverse-stick intercepts with SD one, and otherwise the standard joint update.
This RC choice was made before new LV results: Study 11 showed that a=4,b=4
did not resolve mixing and the exact flat-intercept match separated.

## Frozen MCMC schedule

Each arm uses four genuinely independently seeded chains, 9,000 total
iterations, 3,000 burn-in iterations, 6,000 retained iterations, thinning one,
and four worker cores. Arm seeds and chain seeds are frozen in `spec.R`.
No arm is extended selectively after results are observed.

## Endpoints and convergence interpretation

The primary scientific path is annotation architecture -> PIP/component
architecture -> beta -> prediction. Required outputs include theta and q for
LV; raw alpha and induced component/active priors for RC; common PIP,
component occupancy, beta, variance, validation, convergence, and runtime
metrics for all arms. LV theta has a preferred R-hat target <=1.05, but R-hat
alone is insufficient; bulk/tail ESS, MCSE, and q chain stability are required.
RC interpretation follows raw alpha -> induced prior -> occupancy -> PIP ->
beta/prediction. Poor mixing is a result, not a gate failure, if outputs remain
scientifically assessable.

The seven interpretation rules in `spec.R` are frozen. In particular, stable
theta with unstable q is not success; poor alpha with stable induced priors is
reported as a distinction; and similarity of beta does not prove identical
latent architecture. A hard stop applies only to provenance, canonical input,
matched-operator, missing-capability, or post-result design-change failures.
