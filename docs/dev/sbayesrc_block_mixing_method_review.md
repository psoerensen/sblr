# Block SBayesRC mixing method review

## Scope and provenance

This review concerns learned-alpha scalar SBayesRC with the retained
block-eigen likelihood only. The likelihood remains
`residual_policy = "gctb_block"` with `block_ve_mode = "allMixVe"`.
BED, CSR, block residual-variance policies, mixture priors, and annotation
priors are outside scope.

The development baseline is `sblr` commit
`740b837ed5f794dc39d57f0706f5a47f718f76e4`. Primary source was inspected at:

- standalone SBayesRC v0.2.6, commit
  `b95d3fcbad8ff358290922a58fff893439296138`, especially `src/AnnoProb.cpp`,
  `src/AnnoProb.h`, `src/SBayesRC.cpp`, and `R/sbr_eigen.r`;
- GCTB commit `cc7fa7d765c83a89c6375946cf77fe50ba1a317e`, especially
  `scr/model.cpp`, `scr/model.hpp`, and `scr/gctb.cpp`.

The source snapshots are local read-only evidence under the ignored
`results/local/reference_sources/` directory. No sibling repository was
modified.

## Standalone SBayesRC v0.2.6

`AnnoProb::annoEffect_sample_Gibbs` updates each sequential probit stick
separately. For stick `j`, all markers are used at the first stick and only
markers continuing through stick `j - 1` are used later. It samples one
truncated-normal latent response per eligible marker, then performs a scalar
Gibbs sweep over the intercept and non-intercept annotation coefficients.
`AnnoProb::annoSigmaSS_sample` samples a separate variance per stick from

\[
 \sigma^2_{\alpha,j}=
 \frac{\sum_{k>0}\alpha_{kj}^2+4}{\chi^2_{q+4}},
\]

where `q = numAnno - 1`. `samplePi` then reconstructs marker probabilities.

This is the same basic Albert--Chib/stick-wise transition as current `sblr`,
not a hidden coordinated allocation move. Important model/robustness
differences already documented elsewhere remain: the standalone code uses a
flat intercept, sets an empty stick to intercept -10 with all other
coefficients zero, and uses the production variance prior above. Those are not
sampler improvements for the current proper-prior `sblr` target.

## GCTB ApproxBayesRC

`ApproxBayesRC::sampleUnknowns` performs a SNP allocation/effect update, then
calls `AnnoEffects::sampleFromFC_Gibbs` when annotation probabilities are
learned. That function is again stick-wise Albert--Chib with scalar
coefficient updates. It shuffles non-intercept coefficient order, but does not
jointly update an annotation across sticks and does not make a global or
blockwise allocation/hierarchy move. The optional TGS code is the already
separate high-LD SNP transition and is not an annotation-coupling kernel.

Therefore the ApproxBayesRC implementation provides no same-posterior hidden
block/global allocation transition to transfer.

## GCTB `sampleFromFC_joint` and `sampleFromFC_indep`

The named joint and independent methods belong to `ApproxBayesRD`, not
`ApproxBayesRC`. Both add BayesC-style point-mass inclusion indicators for
non-intercept annotation effects. `sampleFromFC_indep` samples a separate
inclusion state for each annotation-by-stick coefficient. The active
`sampleFromFC_joint` samples one inclusion state for an annotation jointly
across populated sticks and sets that annotation to zero across all sticks
when excluded.

This changes the prior from current `sblr`'s always-included hierarchical
normal coefficients to a spike-and-slab annotation-selection model. It is a
different statistical model, not a posterior-preserving sampler for BayesRC.
It may be a methodological clue that shared annotation selection regularizes
the hierarchy, but it must not be copied into this task.

## Classification

| Mechanism | Classification | Reason |
|---|---|---|
| Stick-wise Albert--Chib update | Same model family; already present | No coordinated allocation movement |
| Random coefficient order | Same posterior; minor scan change | Previous repeated hierarchy updates did not solve joint mixing |
| High-LD TGS | Different transition purpose | Local SNP localization, not global annotation/occupancy feedback |
| ApproxBayesRD independent inclusion | Different model | Adds coefficient-specific BayesC indicators |
| ApproxBayesRD joint inclusion | Different model | Adds one shared annotation indicator across sticks |
| Empty-stick intercept -10 | Heuristic/different contract | Not a draw from the current proper posterior |

## Method sequence

The cheapest defensible same-posterior candidate is a stick-wise sandwich
step in the latent probit augmentation, conditional on the current
`sigmaSqAlpha`. It can be followed by an exact blocked Gaussian alpha draw and
then the unchanged variance update. This directly targets augmentation-induced
alpha persistence without altering `p(c | alpha)` or either prior.

Only if this exact PX/sandwich transition fails the registered 1,500-marker
block mixing gate should a conditional SMC/particle-Gibbs allocation/effect
transition be designed. The different-model GCTB annotation-selection clue is
not a reason to skip that hierarchy.
