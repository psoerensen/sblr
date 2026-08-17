# Phase 7 pattern-by-scale MT-BayesR checkpoint

## Status and literature position

Phase 7 adds a factorized pattern-by-scale extension of the corrected Cheng
MT-BayesC$\Pi$ transition. Phase 7A fixed and independently checked the target,
Phase 7B qualified one shared implementation, and Phase 7C exposes it as
`method = "bayesr"` through `mtblr_bed()`, `mtblr_csr()`, and
`mtblr_block_eigen()`.

The model is a synthesis of the complete binary activity-pattern construction
used by Cheng-style multi-trait Bayesian alphabet models and a BayesR positive
scale ladder. It is not claimed as a novel model, the unrestricted joint-weight
form sometimes called MT-BayesR/BayesR3, or a covariance-template mixture such
as mr.mash. The factorization is deliberate: activity describes trait sharing,
whereas scale describes effect magnitude conditional on a non-null pattern.

Primary positioning sources are Cheng et al.'s
[general multi-trait mixture-prior paper](https://doi.org/10.1534/genetics.118.300650),
the [BayesR3 implementation paper](https://doi.org/10.1038/s42003-022-03624-1),
and the [mr.mash model paper](https://doi.org/10.1371/journal.pgen.1010539).

## Target and state order

For marker $j$, activity pattern $\boldsymbol\delta_j$, and positive scale
$k_j$, the realised effect is

$$
\boldsymbol\alpha_j=D_j\boldsymbol\theta_j,
\qquad
\boldsymbol\theta_j\mid k_j,V_b
\sim N_T(\mathbf 0,\gamma_{k_j}q_jV_b).
$$

There is exactly one null state. Its prior mass is $\Pi_{\mathbf 0}$. Every
non-null state $(\boldsymbol\delta,k)$ has prior mass
$\Pi_{\boldsymbol\delta}\omega_k$, with

$$
\boldsymbol\Pi\sim\operatorname{Dirichlet}(\mathbf a_\Pi),
\qquad
\boldsymbol\omega\sim\operatorname{Dirichlet}(\mathbf a_\omega).
$$

Patterns retain first-trait-fastest Phase 5 ordering. The joint categorical
workspace is ordered as the null state followed by every non-null pattern in
pattern order, with positive scale changing fastest. Pattern and component
assignments remain separate in retained output. No duplicate null component is
created.

## Stored effects and marker conditional

Native chain state stores the unscaled complete base vector

$$
\boldsymbol\beta_j=
\frac{\boldsymbol\theta_j}{\sqrt{\gamma_{k_j}q_j}},
\qquad
\boldsymbol\beta_j\mid V_b\sim N_T(\mathbf0,V_b).
$$

The raw latent-effect field contains the scientifically scaled completed vector
$\boldsymbol\theta_j$; null latent values remain unavailable. For active subset
$A$, likelihood information $h_A,L_{AA}$, and
$g_{jk}=\gamma_kq_j$, the base-coordinate conditional uses

$$
Q_{A,k}=V_{b,AA}^{-1}+g_{jk}L_{AA},
\qquad
u_{A,k}=\sqrt{g_{jk}}h_A,
$$

$$
\boldsymbol\beta_{j,A}\mid-
\sim N(Q_{A,k}^{-1}u_{A,k},Q_{A,k}^{-1}).
$$

The non-null integrated weight is

$$
\log\Pi_{\boldsymbol\delta}+\log\omega_k
-\frac12\log|V_{b,AA}|-\frac12\log|Q_{A,k}|
+\frac12u_{A,k}^{\mathsf T}Q_{A,k}^{-1}u_{A,k}.
$$

Inactive base coordinates are completed under $V_b$ and the complete vector is
then multiplied once by $\sqrt{g_{jk}}$. The null state draws no vector.

## Authoritative updates

The update order is marker sweep, $\boldsymbol\Pi$, $\boldsymbol\omega$, $V_b$,
optional BED $V_e$, convergence capture, then retained capture. Scale counts use
non-null markers only. A deterministic single-component $\boldsymbol\omega$
does not draw from a Dirichlet distribution and therefore consumes no extra RNG.

Because native state stores base vectors, the authoritative statistic is

$$
S_b=\sum_{j:\boldsymbol\delta_j\ne\mathbf0}
\boldsymbol\beta_j\boldsymbol\beta_j^{\mathsf T}
=\sum_{j:\boldsymbol\delta_j\ne\mathbf0}
\frac{\boldsymbol\theta_j\boldsymbol\theta_j^{\mathsf T}}
{\gamma_{k_j}q_j}.
$$

Thus scale is removed exactly once before

$$
V_b\mid-\sim\operatorname{IW}_T(
\nu_{b0}+M_+,\Psi_{b0}+S_b).
$$

BED keeps its fixed or sampled full $V_e$ policy. Independent summary providers
keep fixed positive residual scales and do not infer a full $V_e$.

## Providers, memory, execution, and output

The BED route reuses the packed owner and residual invariant. CSR and
block-eigen routes reuse the single Phase 6A provider kernel. Complete CSR,
full-rank block eigen relative to its declared block-diagonal operator, and
retained-rank reconstruction semantics are unchanged.

Memory adds $M$ component assignments, $M\times K$ marginal component
probabilities, retained assignments, scale-simplex traces, length-$K$ compact
diagnostics, and the simultaneously live pattern-scale conditional workspace.
For each non-null pattern with active-set size $a_\delta$, preflight includes
state tables, active-vector/matrix containers, and numeric storage proportional
to $K\sum_{\delta\ne0}(a_\delta+a_\delta^2)$, multiplied by concurrent chain
workspaces. R preflight applies the fit limit before provider construction and
sampling; a matching native guard protects the immediate allocation. It does
not allocate a mandatory $M\times2^T\times K$ posterior array or a dense joint
transition matrix.

Raw v2 adds component assignments and marker-by-component sub-probabilities
$\Pr(\boldsymbol\delta_j\ne\mathbf 0,k_j=k\mid D)$. These entries sum to the
marker non-null probability, and raw-v2 validation enforces that identity from
the unique declared null activity pattern; a null draw is not assigned a
positive scale.
Raw v2 also adds retained/final/convergence $\boldsymbol\omega$, component
scales, and compact
scale occupancy/change diagnostics. Under the Phase 7 policy, raw-v2 validates
the exact compact pattern and scale fields, axes, integer counts, event totals,
and required absent dense transition matrix. The activity-pattern fields and
$\boldsymbol\Pi$ retain their Phase 5 meaning. Formatted BayesR fits expose
`component_assignments`, `component_probabilities`, `component_scales`,
`omega_trace`, `omega_final`, and `omega_mean`; no ambiguous `pi`, `pis`, or
`pim` alias is added.

## Gate evidence and reductions

The Phase 7A independent R reference enumerates the single null plus all
non-null pattern-by-scale states without calling production code. It verifies
weights, active moments, conditional completion, one-time scale removal,
the one-scale Cheng reduction, the null-plus-all-active multivariate BayesR
reduction, and the mathematical $T=1$ BayesR reduction.

Phase 7B focused qualification checks the native conditional, the exact
$K=1$, $\gamma_1=1$, $q_j=1$ Cheng trajectory, $V_b$ statistic, general-$T$
execution, CSR/full-rank and retained-rank targets, BED/summary diagonal
residual reduction, and serial/parallel reproducibility. Phase 7C checks public
BED, CSR, and block-eigen calls against the same internal routes.

## Deferred work

Unrestricted joint pattern-by-scale weights, trait-specific scale assignments,
overlap-aware summary likelihoods, sampled provider scales, annotation-informed
MT-BayesR, MT-BayesRC/SBayesRC, regional or template covariance models, factor
models, restricted large-$T$ pattern transitions, and marker scheduling remain
deferred.
