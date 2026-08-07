# Particle-marginal global alpha design for retained block SBayesRC

## Scope

This is a feasibility design for the retained block-eigen learned-alpha route.
The block likelihood, `gctb_block` residual variances, effect prior, probit
sticks, proper intercept prior, and `sigmaSqAlpha` prior are unchanged. During
one alpha transition the current block residual variances, common effect
variance, `sigmaSqAlpha`, and all remaining global parameters are conditioned
on. Block allocations and marker effects are marginalized.

## Block target and relative normalizing constant

For block `b`, retained factor `Q_b`, transformed score `w_b`, block residual
variance `Ve_b`, common variance `vb`, marker probabilities `pi_i(alpha)`, and
mixture multipliers `gamma`, define

```text
gamma_b(c,beta; alpha) =
  exp{-||w_b - Q_b beta||^2 / (2 Ve_b)}
  prod_i pi_i,c_i(alpha) p(beta_i | c_i,vb,gamma).
```

The block likelihood is its sum/integral. The factor
`exp{-||w_b||^2/(2 Ve_b)}` is independent of alpha and cancels from every alpha
ratio. The implementation therefore estimates the relative normalizer

```text
Z_b(alpha) = L_b(alpha) / exp{-||w_b||^2/(2 Ve_b)}.
```

Blocks are conditionally independent under the fitted block operator, hence
`Z(alpha) = prod_b Z_b(alpha)`.

## Prefix proposal and unbiased estimator

Visit markers in the fixed block order. With residual `r` after the prefix,
adding `(c_j,beta_j)` contributes

```text
pi_j,c_j(alpha) p(beta_j | c_j)
exp{(beta_j q_j' r - beta_j^2 q_j'q_j / 2) / Ve_b}.
```

Normalizing this expression gives the already validated fully adapted
one-marker component/effect proposal. Its normalizer is `G_j`. A complete
sequential importance path has weight `W = prod_j G_j`. If `P` independent
paths are drawn,

```text
Zhat_b(alpha,u_b) = P^-1 sum_p W_p
```

is unbiased because each path weight is the exact target/proposal ratio. This
reference deliberately does not resample: it avoids discontinuous ancestry
bookkeeping, makes unbiasedness transparent, and permits a simple correlated
randomness representation. It is sequential importance sampling, a valid SMC
normalizing-constant estimator. Independent `u_b` across blocks makes the
product estimator unbiased.

## Auxiliary randomness and selected paths

Every component draw, Gaussian effect draw, and terminal path-selection draw
is generated from a fixed standard-normal auxiliary vector `u_b`; normal CDFs
provide uniforms. Its density is `m_b(u_b)`. The selected terminal path is
drawn with probability proportional to its importance weight. The augmented
estimator identity is

```text
Zhat_b = P^-1 sum_p W_p,
Pr(J_b=p | paths) = W_p / sum_l W_l.
```

Consequently the standard particle-marginal extended target has marginal
`p(alpha | y,theta)`, and conditional on an accepted extended state the
selected paths have the correct allocation/effect distribution. On acceptance
all selected block paths replace `c`, `beta`, and the block residuals. On
rejection alpha, `u`, `Zhat`, and selected paths all remain unchanged.

## PMMH acceptance

For proposal `q(alpha'|alpha)` and independent auxiliary draw `u' ~ m`,

```text
log r = log p(alpha' | sigmaSqAlpha)
      - log p(alpha  | sigmaSqAlpha)
      + log Zhat(alpha',u') - log Zhat(alpha,u)
      + log q(alpha|alpha') - log q(alpha'|alpha).
```

The alpha prior includes each resolved proper intercept normal density and all
non-intercept normal densities with the current stick-specific
`sigmaSqAlpha`. A symmetric random walk cancels the proposal terms.

## Correlated pseudo-marginal proposal

Stack the independent block auxiliaries into `u ~ N(0,I)`. The proposal

```text
u' = rho u + sqrt(1-rho^2) epsilon,  epsilon ~ N(0,I)
```

is reversible with respect to `m(u)` and preserves every block's independent
standard-normal marginal. Using it jointly with an alpha proposal therefore
leaves the extended PMMH target invariant. Correlation is assessed through
the variance of the log-estimator difference, not merely marginal estimator
variance.

## Hard gates

Before any scientific transition: exact enumeration must validate `Z_b`, the
likelihood-scale estimator must be unbiased, selected-path marginals must
match the exact target, and estimator/proposal noise must be practical for all
15 small blocks. No production integration is permitted unless meaningful
stick-wise alpha moves have non-negligible predicted acceptance and selected
paths reorganize global occupancy at defensible cost.
