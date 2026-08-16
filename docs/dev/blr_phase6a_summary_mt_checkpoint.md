# Phase 6A independent-summary Cheng MT checkpoint

## Status

**Status:** `CURRENT_QUALIFICATION`.

Phase 6A implements an internal, qualification-only general-$T$ Cheng
MT-BayesC$\Pi$ sampler for declared independent summary providers. It does not
export `mtblr_csr()` or `mtblr_block_eigen()`, alter public `mtblr_bed()`, or
make an overlap-aware likelihood available.

The implementation reuses the Phase 2 global-marker, provider, resource, and
operator contracts; the Phase 3 task, seed, retention, convergence, and worker
contracts; and the Phase 5A complete-pattern Cheng state and covariance
transitions.

## Likelihood contract

Provider $p$ owns one trait $t(p)$, a local score $s_p=X_p^\top y_p$, a
declared cross-product operator $C_p=X_p^\top X_p$, a positive fixed residual
scale $\phi_p$, and a provider-local-to-global marker map. Its contribution is

$$
\log L_p=-\frac{1}{2\phi_p}
\left(\alpha_p^\top C_p\alpha_p-2\alpha_p^\top s_p\right)
+\text{constant}.
$$

Provider contributions add. Every provider declares
`likelihood_regime = "independent_summary"`, owns exactly one trait, and uses
`residual_contract = "fixed_provider_residual_scale"`. A non-`NULL` overlap
group is rejected because Phase 6A does not model cross-provider score-error
covariance.

Provider residual scores are maintained as

$$
r_p=s_p-C_p\alpha_p.
$$

For global marker $j$, the old realised effect is restored in every provider
that contains the marker. The likelihood information for trait $t$ is

$$
h_{jt}=\sum_{p:t(p)=t,\,j\in p}\frac{r_{p,j}^{(-j)}}{\phi_p},
\qquad
d_{jt}=\sum_{p:t(p)=t,\,j\in p}\frac{C_{p,jj}}{\phi_p}.
$$

A marker absent from a provider contributes neither score nor precision. It
is not converted into a zero-effect observation.

## Cheng transition

The complete canonical $2^T$ binary pattern order is unchanged: the first
declared trait changes fastest. For active set $A$,

$$
Q_{j,A}=V_{b,AA}^{-1}+\operatorname{diag}(d_{j,A}),
\qquad
m_{j,A}=Q_{j,A}^{-1}h_{j,A}.
$$

The integrated log weight is

$$
\log w_{j,A}=\log\Pi_A
-\frac12\log|V_{b,AA}|
-\frac12\log|Q_{j,A}|
+\frac12 h_{j,A}^\top Q_{j,A}^{-1}h_{j,A}.
$$

The null log weight is $\log\Pi_\varnothing$. One joint categorical pattern is
drawn. Active coordinates are drawn from
$N(m_{j,A},Q_{j,A}^{-1})$; inactive coordinates of a selected non-null marker
are completed from the established Gaussian Schur conditional under $V_b$.
The null remains collapsed and consumes no latent-effect RNG.

After each sweep, $\boldsymbol\Pi$ is updated from its Dirichlet conditional
and the sole marker covariance state is updated as

$$
V_b\mid-\sim\operatorname{IW}_T
\left(\nu_{b0}+M_+,\Psi_{b0}+
\sum_{j:\delta_j\ne0}\beta_j\beta_j^\top\right).
$$

Only completed non-null latent vectors enter the statistic. Provider residual
scales are fixed; Phase 6A does not sample $V_e$ or provider variances.

## Shared operator implementation

One native chain kernel serves both representations. Immutable operator
resources are constructed once before the OpenMP chain region and shared by
all tasks. Chain-local provider residual vectors are the only mutable operator
state.

- CSR uses the maintained `SparseLdCsrStorage`/`SparseLdCsrView` off-diagonal
  action and its separately validated positive diagonal. CSR thresholding is
  recorded as an approximation.
- Block eigen retains provider-supplied $U\sqrt{\Lambda}$ factors and applies
  columns within their declared blocks. It does not densify the operator.
  Full rank is exact for the declared block-diagonal operator. Retained rank
  targets the reconstructed retained operator.

Providers may have different marker subsets, local order, sample size,
population, CSR sparsification, block partition, and retained rank. Resources
are not summed into a global matrix and are not copied per chain.

## Execution and memory

One complete joint chain is one Phase 3 logical task. Static scheduling,
uint32 task seeds, post-burn divisible retention, unthinned convergence, and
sampler-region worker diagnostics are unchanged. Diagnostics and capture
consume no RNG.

The checked Phase 6A preflight starts from the Phase 5A general-$T$ output
contract and adds immutable operator storage, provider score/map metadata,
marker-provider references, chain-private residual-score vectors, and the
native input-conversion peak. The mandatory marker-by-pattern output continues
to scale as $M2^T$. Checked R and native arithmetic reject overflow or a memory
limit violation before native operator construction.

## Raw-v2 semantics

The internal entry point is
`.blr_cheng_mt_bayesc_summary_qualification()`. It returns validated `blr_raw`
version 2 with dynamic draw, chain, marker, trait, covariance, and activity-
pattern axes. It records provider metadata, local/global maps, fixed residual
scales, operator approximation status, exact task seeds and indices, compact
occupancy/change diagnostics, and the authoritative $V_b$ states.

Individual-level predictions, residual-covariance draws, and residual-
covariance posterior summaries are required-present `NULL`. Approximate
operator quadratics are not reported as SSE or genomic variance.

## Qualification evidence

Focused permanent checks cover:

- independent general-$T$ marker weights, means, and covariances;
- complete CSR and full-rank block-eigen reductions;
- retained-rank agreement with its reconstructed operator;
- heterogeneous marker subsets, local order, and provider order;
- additive splitting of genuinely independent same-trait providers;
- the diagonal-residual common-sample marker-conditional and retained-summary
  reduction;
- $T=2$ and $T=3$ short chains;
- serial/two-chain parallel identity;
- overlap rejection, memory accounting, and raw-v2 semantics.

## Deferred boundary

Phase 6A does not expose public MT CSR or block-eigen functions. Public
promotion is a separate Phase 6B decision. Also deferred are overlap-aware
summary likelihoods, common-sample covariance recovery from marginal
summaries, sampled provider scales, missing or overlapping phenotypes,
MT-BayesR/BayesRC, annotations, $q_j$, MAF-S, regional covariance, covariance
templates, and restricted large-$T$ pattern spaces.
