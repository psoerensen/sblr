# Coordinated BayesRC alpha--allocation transition

## Scope

This note derives the current scalar BayesRC/SBayesRC posterior terms needed
for an exact coordinated annotation/allocation transition. The likelihood,
spike-and-slab prior, probit sticks, alpha priors, `sigmaSqAlpha` hierarchy,
ordinary Gibbs transitions, and route-specific residual operators are fixed.
The development baseline is `sblr` commit
`0c89234273389e14112ba0e08ef9d11d3e1819dc`.

No coordinated transition may be integrated until its support, proposal
density, and stationary target pass the gates below. In particular, the
derivation distinguishes an allocation-label proposal from a valid move on
the mixed discrete/continuous state `(c, beta)`.

## Current ordered probit-stick prior

There are `K` mixture components indexed `0, ..., K - 1` and `K - 1` sticks.
For marker `i` and stick `j`, with the annotation row `A_i` including its
all-ones intercept column,

\[
p_{ij}=\Phi(A_i\alpha_j).
\]

The native implementation treats `p_ij` as continuation probability. Before
the package's documented positive probability floor and row normalization,

\[
\Pr(c_i=k\mid\alpha)=
\begin{cases}
(1-p_{ik})\prod_{h<k}p_{ih}, & k<K-1,\\
\prod_{h=0}^{K-2}p_{ih}, & k=K-1.
\end{cases}
\]

The stick response is

\[
y_{ij}=I(c_i>j),
\]

and marker `i` is eligible for stick `j` exactly when `j = 0` or
`c_i > j - 1`, equivalently `c_i >= j`. Thus a marker stops at its component
index and is absent from every later stick likelihood. This is the exact
zero-based contract in `st_bayesrc_annotation_prior.h`.

## Relevant posterior factorization

Let `theta` collect residual and marker-effect variances and all remaining
state. With the spike represented by a point mass at zero,

\[
\begin{aligned}
 p(\beta,c,\alpha,\sigma_\alpha^2,\theta\mid y)
 &\propto p(y\mid\beta,\theta)
 \prod_i p(\beta_i\mid c_i,\theta)\Pr(c_i\mid\alpha)\\
 &\quad\times\prod_j p(\alpha_j\mid\sigma_{\alpha,j}^2)
 p(\sigma_{\alpha,j}^2)p(\theta),
\end{aligned}
\]

where

\[
p(\beta_i\mid c_i=0)=\delta_0(\beta_i),\qquad
p(\beta_i\mid c_i=k>0)=N(0,v_b\gamma_k s_i).
\]

`s_i` is one unless the already existing MAF-effect scaling contract is
active. If only stick `j` changes, the alpha target contribution is its proper
intercept density times the existing zero-mean Gaussian densities for the
non-intercept coefficients. `sigmaSqAlpha_j` is held fixed by the proposed
transition. Every allocation probability changes through stick `j`, including
markers outside any proposed allocation subset.

For a proposed state `(alpha_j', c', beta')` with all other state fixed, the
non-cancelling target terms are therefore:

- the likelihood change caused by `beta'`;
- the spike or Gaussian effect-prior terms for changed `(c_i, beta_i)`;
- allocation probabilities under `alpha_j'` for every marker;
- the proper alpha prior density for stick `j`.

The `sigmaSqAlpha` prior, other alpha sticks, other variance priors, and fixed
route/operator terms cancel.

## Candidate A: alpha plus allocation labels, beta fixed

Suppose `alpha_j'` has evaluable density
`q_alpha(alpha_j' | alpha_j, c)` and component labels in subset `S` have
evaluable density `q_c(c_S' | alpha_j', beta, c_-S)`. The formal log ratio is

\[
\begin{aligned}
\log r_A={}&\log p(\alpha_j'\mid\sigma_{\alpha,j}^2)
-\log p(\alpha_j\mid\sigma_{\alpha,j}^2)\\
&+\sum_i\{\log\Pr(c_i'\mid\alpha')-
               \log\Pr(c_i\mid\alpha)\}\\
&+\sum_{i\in S}\{\log p(\beta_i\mid c_i')-
                       \log p(\beta_i\mid c_i)\}\\
&+\log q_\alpha(\alpha_j\mid\alpha_j',c')-
  \log q_\alpha(\alpha_j'\mid\alpha_j,c)\\
&+\log q_c(c_S\mid\alpha_j,\beta,c_{-S})-
  \log q_c(c_S'\mid\alpha_j',\beta,c_{-S}).
\end{aligned}
\]

This expression exposes a support failure, not merely poor tuning. If
`beta_i != 0`, the reverse target density for `c_i = 0` is zero. If
`beta_i = 0` under the point-mass component, proposing an active component
while retaining exactly zero lies in a set of zero measure under its Gaussian
effect prior. Consequently a beta-fixed proposal cannot reversibly cross the
null/non-null boundary. Restricting it to labels `c_i > 0` is valid but
preserves active count and cannot address the demonstrated bottleneck.

Candidate A is therefore mathematically pathological for the requested
occupancy transition and is rejected before implementation.

## Candidate B: collapse alpha

For one stick the allocation likelihood is a probit-regression likelihood,

\[
L_j(\alpha_j)=\prod_{i:c_i\ge j}
 \Phi(A_i\alpha_j)^{I(c_i>j)}
 \{1-\Phi(A_i\alpha_j)\}^{I(c_i=j)}.
\]

Multiplying by the proper Gaussian alpha prior does not yield a closed-form
integral. Albert--Chib augmentation makes `alpha_j | z_j` Gaussian, but
integrating the truncated latent variables requires a multivariate-normal
orthant probability whose dimension is the eligible marker count. It does not
factor marker by marker after alpha is integrated. Exact quadrature is useful
for a tiny oracle but is not scalable to the Study 06 annotation dimension and
eligible counts.

Candidate B has no tractable exact production form under the unchanged model.

## Candidate C: alpha plus allocations and effects

Let `S` be selected independently of the current state. Remove `beta_S` from
the fitted value and write `G_S` and `s_S` for the fitted-operator submatrix and
corrected score. For each component configuration `d in {0,...,K-1}^|S|`, let
`D_d` be the active subset and `V_d` its diagonal effect-prior covariance.
Conditional on `alpha`, the collapsed weight is

\[
\log w_d(\alpha)=\sum_{i\in S}\log\Pr(c_i=d_i\mid\alpha)
-\tfrac12\{\log|V_d|+\log|P_d|\}
+\tfrac12h_d^TP_d^{-1}h_d,
\]

with

\[
P_d=G_d/v_e+V_d^{-1},\qquad h_d=s_d/v_e.
\]

The null configuration has zero determinant and quadratic contributions.
After sampling `d` from normalized weights, active effects are sampled from
`N(P_d^{-1}h_d, P_d^{-1})`; null effects are exactly zero. This is the current
one- and two-marker conditional generalized to `|S|` markers.

Consider the proposal:

1. choose stick `j` and fixed subset `S` with state-independent probability;
2. propose `alpha_j' ~ q_alpha(alpha_j' | alpha_j)`;
3. draw `(c_S', beta_S')` from the exact collapsed conditional under
   `alpha_j'`, holding the outside state fixed.

Let

\[
Z_S(\alpha)=\sum_d\exp\{\log w_d(\alpha)\}.
\]

All likelihood, effect-prior, and allocation-prior terms inside `S` cancel
against their exact conditional proposal. The exact MH ratio is

\[
\begin{aligned}
\log r_C={}&
\log p(\alpha_j'\mid\sigma_{\alpha,j}^2)
-\log p(\alpha_j\mid\sigma_{\alpha,j}^2)\\
&+\sum_{i\notin S}
 [\log\Pr(c_i\mid\alpha')-\log\Pr(c_i\mid\alpha)]\\
&+\log Z_S(\alpha')-\log Z_S(\alpha)\\
&+\log q_\alpha(\alpha_j\mid\alpha_j')-
  \log q_\alpha(\alpha_j'\mid\alpha_j).
\end{aligned}
\]

For a symmetric fixed-covariance Gaussian random walk, the final line is zero.
This construction is exact, has reversible null/non-null support, and includes
beta consistently.

## Scalability and selection gate

Exact evaluation requires `K^|S|` component configurations and a Cholesky
factorization for each active subset. For the Study 06 four-component model:

| subset size | configurations |
|---:|---:|
| 1 | 4 |
| 2 | 16 |
| 4 | 256 |
| 8 | 65,536 |
| 10 | 1,048,576 |
| 20 | about 1.1 trillion |

Small `S` is computationally feasible but leaves the allocation-prior change
for nearly every marker in the outside-state sum. It therefore cannot absorb a
global alpha/sparsity shift. A large enough `S` to coordinate the demonstrated
occupancy regimes is exponentially intractable. State-dependent selection of
only discordant markers does not remove this problem and would additionally
require its exact reverse selection probability.

The exact Candidate C construction will be implemented only as an independent
tiny-model feasibility oracle. Production integration is conditional on
showing non-negligible acceptance and mode crossing at a subset size whose
enumeration cost is scalable. No public sampler control is added merely to
expose an exact but unusable transition.

## Validation plan

The independent R oracle will verify:

- component probability and eligible-stick mappings;
- Candidate A's null/non-null support failure;
- Candidate C collapsed weights, forward/reverse proposal densities, and MH
  ratios;
- detailed balance for fixed state pairs;
- stationary allocation, alpha, and beta moments against enumeration plus
  quadrature;
- acceptance and cost as marker count grows while exact subset size remains
  fixed.

If the feasibility gate fails, the result is AA-R5 rather than a heuristic
production transition. Ordinary BED, CSR, and block-eigen trajectories then
remain bit-for-bit unchanged by construction.
