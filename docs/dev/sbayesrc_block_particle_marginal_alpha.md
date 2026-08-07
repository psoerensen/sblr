# Block particle-marginal alpha target

The canonical derivation and validation plan for the block-factorized
particle-marginal alpha target is in
[`sbayesrc_particle_marginal_alpha_design.md`](sbayesrc_particle_marginal_alpha_design.md).

That design conditions on the current block residual variances, effect
variance, annotation variance, and remaining global state; marginalizes each
block's allocation and effect state; uses an unbiased likelihood-scale SMC
normalizing-constant estimate; and retains a selected terminal path in the
standard PMMH extended target. It is a development reference only and does not
change the production block SBayesRC sampler.

