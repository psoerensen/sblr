# Study 12 — SBayesR-LV annotation evaluation, Aspect 2

Status: preregistered on 2026-08-10 before generating the new phenotype,
summary statistics, or inspecting any Aspect-2 model fit.

## Scientific question

Aspect 2 is the reciprocal of completed Aspect 1. It asks whether qualified
SBayesR-LV recovers and exploits a known log-variance annotation architecture,
and how ordinary SBayesR and the preregistered SBayesRC model approximate that
architecture. The only fitted arms are R, LV, and RC. This single scenario is
not a universal model ranking.

## Frozen implementation and common backbone

Every arm uses `sblr` 0.2.0 at clean sibling SHA
`2123699a9cc2e91059e7d81a745420b14eca7f6e`, loaded from the Aspect-1 isolated
library. The common backbone is `canonical_annotation_simulation_v1`: the
same 5,000 people and 4,500/500 split, 35,000 ordered markers, 35 contiguous
1,000-marker blocks, raw annotation realization, training allele frequencies,
working Glist, and LD reference used in Aspect 1. Genotypes, samples,
annotations, block definitions, and LD are not regenerated.

Aspect-1 sources and local evidence are protected by the SHA-256 inventory
`results/local/12_logvar_annotation_evaluation/provenance/aspect1_preservation_snapshot.csv`.
Aspect-2 outputs use a separate `aspect2/` subtree and cannot overwrite them.

## Frozen LV truth

Production LV preprocessing is applied exactly once to the same three raw
nonintercept annotations: binary is centered only, continuous columns are
centered and divided by their sample SD, and there is no intercept. The true
coefficient vector is `(log(3), log(2), 0)` for enriched binary, continuous
signal, and null annotation. Thus `eta_true=X theta_true`, `q_true=exp(eta_true)`,
the intended conditional variance ratios are 3, 2, and 1, and mean eta / the
geometric mean of q must be numerically zero / one.

The global BayesR gamma vector is `(0, 0.01, 0.1, 1)`. Exact component counts
are 34,000/500/300/200. A new seed, 1215412, randomly permutes this fixed label
multiset independently of annotations; annotation values never enter the
allocation probabilities. Before fitting, construction independence and
descriptive component/annotation correlations are recorded. The registered
sanity limits are maximum absolute component/annotation correlation 0.10 and
maximum absolute active-versus-null standardized annotation mean difference
0.20; these are guards against alignment/coding mistakes, not hypothesis tests.

With effect seed 1215424, each active raw effect is drawn independently as
`N(0, gamma[c_j] q_true[j])`. A single global multiplier then sets training
genetic variance to one. This preserves gamma and q ratios and, together with
the exact frozen Study-10 standardized residual vector (original seed 1005438,
training variance one), gives realized h2 0.5. The Aspect-2 phenotype and
training summary statistics are generated once and shared by every arm.

## Operator, priors, chains, and seeds

All arms use the same retained block-eigen route as Aspect 1: cumulative
positive eigenvalue mass 0.995, expected 984 modes in every block,
`residual_policy="gctb_block"`, `block_ve_mode="allMixVe"`, resampling
threshold 1.1, minimum Ve ratio 0.7, block histories, and residual rebuild
every 100 iterations. Rank is never optimized by arm.

Each arm uses four independent chains, 9,000 total iterations, 3,000 burn-in,
6,000 retained draws, thinning one, and four cores. R2 chain seeds are
1226141/1226242/1226343/1226444; LV2 seeds are
1227141/1227242/1227343/1227444; RC2 seeds are
1228141/1228242/1228343/1228444. LV fixes `theta_prior_sd=0.7`. RC exactly
reuses Aspect 1's proper-intercept, a=2/b=2 production/reference configuration.
No arm is shortened, extended, restarted, or tuned after results.

## Endpoints and interpretation

Primary LV endpoints are theta truth coverage/bias and marker-level q truth
recovery on raw and log scales, separately from q cross-chain stability.
Common endpoints preserve Aspect-1 definitions for PIP/components, beta,
prediction, variance quantities, MCMC stability, and runtime. RC is assessed
through alpha -> induced prior -> occupancy -> PIP -> beta/prediction.

If directly available, RC's marker expected prior second moment is
`sum_{k>0} Pr(c_j=k|A) gamma_k`, normalized to arithmetic mean one before
comparison with q truth. No substitute proxy is allowed.

The seven Outcome A--G rules in `spec_aspect2.R` are frozen. In particular,
plausible theta without q recovery is not success; q recovery with uncertain
theta is reported as induced-architecture recovery; reciprocal own-model
advantages imply architecture dependence, not universal superiority. Poor
mixing is a result. Hard stops are limited to provenance, frozen-backbone,
truth-construction, matched-input, operator, or genuine execution failures.

For the registered descriptive screen, q recovery requires raw and log-scale
Pearson/Spearman at least 0.95 and log-q RMSE at most 0.25. A beta-truth or
validation-genetic-value correlation gain of at least 0.02 is called material;
absolute differences at most 0.01 are called essentially equal. Exact values,
uncertainty, chain diagnostics, and all downstream endpoints remain the
scientific evidence; these screens are not universal calibration thresholds.
