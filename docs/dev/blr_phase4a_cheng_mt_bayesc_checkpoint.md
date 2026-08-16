# Phase 4a Cheng MT-BayesC$\Pi$ checkpoint

## Status and boundary

**Status:** `IMPLEMENTED - READY FOR INDEPENDENT VERIFICATION`

This checkpoint records one qualification-only vertical slice. It is not a
public MTBLR replacement. The internal route supports exactly two traits, one
common-sample packed-BED provider, complete phenotypes, and a supplied fixed
full residual covariance. The legacy public MT routes and their covariance
hybrid are unchanged and are not called by this implementation.

The slice reuses:

- the Phase 2 global-marker, BED-resource, and non-factorized multi-trait
  provider contracts;
- `PackedBedMatrix`, `BedPackedGenotypeView`, and the established MT BED
  decoder and standardization maps;
- the coherent inverse-Wishart draw helper, not the legacy MT transition;
- the Phase 3 joint-multitrait chain task, uint32 seed, retention, static
  scheduling, convergence, and worker-diagnostic contracts;
- the Phase 1 resolved specification and `blr_raw` version-2 envelope.

No complete genotype matrix or cross-product matrix is materialized. Phase 4a
does not implement sampled $V_e$, summary-statistic MT likelihoods, MT-BayesR,
regional covariance, or a public route.

## Scientific target

For common individuals and two ordered traits,

$$
Y=XA+E,
\qquad
E_{i\cdot}\sim N_2(\mathbf 0,V_e),
$$

where $V_e$ is fixed, finite, symmetric, positive definite, and ordered by the
declared trait IDs. For marker $j$,

$$
\boldsymbol{\alpha}_j=D_j\boldsymbol{\beta}_j,
\qquad
D_j=\operatorname{diag}(\boldsymbol{\delta}_j).
$$

The canonical activity-pattern order is fixed as

$$
(0,0),\quad(1,0),\quad(0,1),\quad(1,1).
$$

For every non-null pattern,

$$
\boldsymbol{\beta}_j\mid V_b\sim N_2(\mathbf0,V_b),
$$

and

$$
\Pr(\boldsymbol{\delta}_j=\boldsymbol{\delta}^{(k)})=\Pi_k,
\qquad
\boldsymbol{\Pi}\sim\operatorname{Dirichlet}(\mathbf a_\Pi).
$$

The marker scale is $q_j=1$ for every marker. There is no component scale,
annotation multiplier, region, or covariance template.

## Marker transition

For marker $j$, the sampler first restores its current realised contribution,
giving

$$
R_j=Y-\sum_{\ell\ne j}
\mathbf x_\ell\boldsymbol{\alpha}_\ell^\top.
$$

It then forms

$$
c_j=\mathbf x_j^\top\mathbf x_j,
\qquad
\mathbf s_j=R_j^\top\mathbf x_j,
\qquad
\Omega_e=V_e^{-1}.
$$

For active coordinates $A$,

$$
P_A=V_{b,AA}^{-1}+c_j\Omega_{e,AA},
\qquad
\mathbf h_A=(\Omega_e\mathbf s_j)_A.
$$

The active effect conditional is

$$
\boldsymbol{\beta}_{j,A}\mid-
\sim N(P_A^{-1}\mathbf h_A,P_A^{-1}),
$$

and the joint categorical log weight, up to a common constant, is

$$
\log w_{\boldsymbol{\delta}}
=
\log\Pi_{\boldsymbol{\delta}}
-\frac12\log|V_{b,AA}|
-\frac12\log|P_A|
+\frac12\mathbf h_A^\top P_A^{-1}\mathbf h_A.
$$

The null log weight is $\log\Pi_{(0,0)}$. Weights are normalized with
log-sum-exp and one joint categorical state is drawn. The completed realised
effect is then removed from the residual. Markers are traversed in declared
BED/resource order.

The deterministic native marker contract agrees with an independently coded R
calculation of all four probabilities, conditional means, and conditional
covariances to at most $10^{-13}$ on the qualification fixture. Replacing the
off-diagonal $V_e$ by its diagonal changes the probabilities, establishing
that the full residual covariance participates in the transition.

## Null collapse and conditional completion

For a partial pattern with inactive set $I$, Phase 4a completes the latent
vector using

$$
\boldsymbol{\beta}_{j,I}\mid\boldsymbol{\beta}_{j,A},V_b
\sim N\left(
V_{b,IA}V_{b,AA}^{-1}\boldsymbol{\beta}_{j,A},
V_{b,II}-V_{b,IA}V_{b,AA}^{-1}V_{b,AI}
\right).
$$

For $(1,1)$ it samples the complete likelihood-informed vector. For $(0,0)$
it does not draw a latent vector: the realised effect is exactly zero, the
latent draw is represented as unavailable (`NA`) on its fixed raw array axes,
and the marker is excluded from the covariance update. A zero vector is never
reported as a sampled null latent vector.

For every non-null marker the raw validator checks

$$
\boldsymbol{\alpha}_j=D_j\boldsymbol{\beta}_j.
$$

## Authoritative marker-covariance transition

The sole marker-covariance update is

$$
V_b\mid-
\sim\operatorname{IW}_2\left(
\nu_0+M_+,
\Psi_0+
\sum_{j:\boldsymbol{\delta}_j\ne\mathbf0}
\boldsymbol{\beta}_j\boldsymbol{\beta}_j^\top
\right).
$$

The parameterization is degrees of freedom followed by scale, with density
proportional to

$$
|V_b|^{-(\nu+3)/2}
\exp\left\{-\frac12\operatorname{tr}(\Psi V_b^{-1})\right\}.
$$

The prior is proper when $\nu_0>1$; the stronger finite-mean condition is
$\nu_0>3$ and is not silently imposed. The prior scale, initial state, and
every sampled state are checked for finite symmetry and strict positive
definiteness.

The draw produced here is the authoritative chain state. It governs the next
marker sweep and is captured without replacement in retained draws,
convergence traces, posterior summaries, and the final state. Diagnostic
sufficient-statistic checks reconstruct the last update from completed
non-null latent effects. No heuristic covariance is calculated later.

## Pattern-probability transition and update order

After each complete marker sweep,

$$
\boldsymbol{\Pi}\mid-
\sim\operatorname{Dirichlet}(\mathbf a_\Pi+\mathbf n_\delta),
$$

where counts follow the fixed four-pattern order. The frozen iteration order
is:

1. joint marker transitions in declared marker order;
2. the Dirichlet pattern-probability update;
3. the inverse-Wishart marker-covariance update;
4. completed-iteration convergence capture;
5. retained-draw capture when the post-burn transition index is retained.

Internal controls may fix $V_b$ or $\boldsymbol{\Pi}$ for qualification. Their
resolved prior metadata distinguishes fixed from sampled states.

## Execution and raw output

The analysis mode is `joint_multitrait`; a logical task is one whole chain.
Parallelism is across chains only. Each task owns its uint32 seed, one
`std::mt19937`, effects, activity states, residual, $V_b$, $\boldsymbol{\Pi}$,
captures, and diagnostics. The OpenMP loop uses static scheduling and task
results are assembled in canonical chain order. Worker code does not call the
R API or construct Rcpp objects.

The qualification `blr_raw` v2 object preserves these principal axes:

| Quantity | Axes |
|---|---|
| realised and latent effects | draw $\times$ chain $\times$ marker $\times$ trait |
| joint activity state | draw $\times$ chain $\times$ marker |
| traitwise activity | draw $\times$ chain $\times$ marker $\times$ trait |
| $V_b$ | draw $\times$ chain $\times$ trait-row $\times$ trait-column |
| $\boldsymbol{\Pi}$ | draw $\times$ chain $\times$ activity-pattern |
| predictions | draw $\times$ chain $\times$ observation $\times$ trait |
| final realised and latent effects | chain $\times$ marker $\times$ trait |
| final $V_b$ | chain $\times$ trait-row $\times$ trait-column |
| final $\boldsymbol{\Pi}$ | chain $\times$ activity-pattern |

The object also contains markerwise pattern probabilities, traitwise PIPs,
pleiotropic probabilities, exact retained and convergence indices, task IDs,
task seeds, actual sampler-worker diagnostics, fixed $V_e$, prior metadata,
and explicit `qualification_only` provenance. There is no ambiguous `pi`,
`pis`, or `pim` field and no sampled $V_e$ draw.

Here `posterior$pleiotropic_probabilities` is required and has the marker axis
in declared global-marker order. Every value is finite and lies in $[0,1]$.
It equals the `posterior$activity_pattern_probabilities` column whose declared
activity-pattern ID is `1_1`, corresponding to $(1,1)$; validation aligns by
that ID rather than assuming a column number and applies the raw-object
consistency tolerance $10^{-12}$. This conditional raw-v2 field is not required
for models where a pleiotropic joint activity pattern is scientifically
undefined.

## Validation evidence

The focused permanent suite provides deterministic algebra, raw-schema,
covariance-coherence, validation, fixed-control, and concurrency checks. In a
one-marker fixed-$V_b$, fixed-$\boldsymbol{\Pi}$ case, independent exact pattern
probabilities were

$$
(0.74356433,\ 0.09665351,\ 0.07082904,\ 0.08895312),
$$

while 12,000 native retained draws gave

$$
(0.74525000,\ 0.09966667,\ 0.07016667,\ 0.08491667),
$$

with maximum absolute difference $0.00403645$. The independently maintained
research suites continue to reproduce their exact/numerical Cheng targets,
including the documented learned-pattern posterior mean
$(0.547150,0.093526,0.092124,0.267200)$.

In a separate 12,000-draw comparison with both $V_b$ and
$\boldsymbol{\Pi}$ sampled, the completed-active R reference and native BED
sampler gave $V_b$ posterior means differing by at most $0.002177213$ and
$\boldsymbol{\Pi}$ posterior means differing by at most $0.003880714$. This is
a Monte Carlo comparison between independent implementations, not an exact
oracle and not a shared-RNG trajectory test.

Repeated runs and one-worker versus two-worker runs with identical final task
seeds produce identical scientific arrays and final states. Keeping or
discarding convergence traces changes no draw or final state. Retained
covariance draws are finite, symmetric, and positive definite; pattern
simplexes normalize; and null latent values follow the explicit unavailable
contract.

## Deferred work

Phase 4b has added sampled full inverse-Wishart $V_e$ as a separate
qualification policy; see
[the Phase 4b checkpoint](blr_phase4b_sampled_residual_covariance_checkpoint.md).
Independent verification and later public-promotion gates still precede any
maintained replacement. Summary
operators, heterogeneous or overlapping providers, MT-BayesR/MT-BayesRC,
regional covariance, covariance-template mixtures, structured large-trait
states, and public migration remain separately gated.
