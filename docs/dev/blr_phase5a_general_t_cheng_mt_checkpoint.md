# Phase 5A general-$T$ Cheng MT checkpoint

Status: implemented, qualification-only, pending independent verification.

## Scope and source reconstruction

Phase 5A generalizes the corrected common-sample packed-BED Cheng
MT-BayesC$\Pi$ route from exactly two traits to a dynamically declared modest
$T\geq2$. The route remains internal through
`.blr_cheng_mt_bayesc_bed_qualification()`. No public MT or ST dispatch,
default, export, or legacy result meaning changes. Public replacement is Phase
5B work.

The legacy `.mtblr_models()` helper already enumerates all binary patterns with
the first declared trait changing fastest. That ordering and its stable IDs are
reused. The legacy covariance transition, heuristic final covariance
replacement, regional/set behavior, and ambiguous latent-versus-realised
output are not reused. Phase 5A instead reuses the qualified packed-BED view,
Phase 2 resource/provider contract, generic degrees-of-freedom/scale
inverse-Wishart RNG, Phase 3 execution contract, resolved specification, and
raw-v2 validator.

The Phase 4 kernel's fixed two-element arrays, four-state containers,
two-column residual loops, and two-dimensional serializers are replaced by one
dynamic kernel. A duplicate two-trait sampler is not retained.

## Activity-pattern contract

For marker $j$,

$$
\delta_j\in\{0,1\}^T,\qquad
D_j=\operatorname{diag}(\delta_j),\qquad
\alpha_j=D_j\beta_j.
$$

State $k\in\{0,\ldots,2^T-1\}$ uses bit $t$ for declared trait $t$, so the
first trait changes fastest. The $T=2$ order remains exactly `(0,0)`, `(1,0)`,
`(0,1)`, `(1,1)`. The binary matrix, row IDs, Dirichlet parameters, native
states, raw axes, and summaries are validated as one aligned contract, with
exactly one null and one all-active row.

## Marker conditional and completion

With marker-excluded residual $R_j$, $c_j=x_j^\top x_j$,
$s_j=R_j^\top x_j$, and $\Omega_e=V_e^{-1}$, active subset $A$ uses

$$
P_A=V_{b,AA}^{-1}+c_j\Omega_{e,AA},\qquad
h_A=(\Omega_es_j)_A.
$$

The active draw has mean $P_A^{-1}h_A$, covariance $P_A^{-1}$, and log weight

$$
\log w_A=\log\Pi_A-
\frac12\log|V_{b,AA}|-
\frac12\log|P_A|+
\frac12h_A^\top P_A^{-1}h_A.
$$

The null weight is $\log\Pi_0$. Log-sum-exp normalization precedes one joint
categorical draw. Singleton, multidimensional partial, and all-active subsets
use Cholesky operations. For inactive subset $I$,

$$
\beta_{j,I}\mid\beta_{j,A},V_b\sim N\left(
V_{b,IA}V_{b,AA}^{-1}\beta_{j,A},
V_{b,II}-V_{b,IA}V_{b,AA}^{-1}V_{b,AI}
\right).
$$

Production completion uses factorized solves and validates the symmetric Schur
complement without arbitrary jitter. Null markers draw no latent vector and
store unavailable latent entries. Completed non-null vectors alone enter the
$V_b$ statistic.

## Authoritative covariance transitions and schedule

The marker covariance update is

$$
V_b\mid-\sim\operatorname{IW}_T\left(
\nu_{b0}+M_+,
\Psi_{b0}+\sum_{j:\delta_j\ne0}\beta_j\beta_j^\top
\right).
$$

When enabled, the completed-sweep residual update is

$$
V_e\mid-\sim\operatorname{IW}_T\left(
\nu_{e0}+N,
\Psi_{e0}+(Y-XA)^\top(Y-XA)
\right).
$$

Both priors use degrees of freedom and scale. Propriety requires
$\nu_0>T-1$; a finite prior mean requires $\nu_0>T+1$, but the latter is not
imposed. Each sampled matrix immediately becomes the next-sweep state and the
same state is captured in convergence, retained draws, summaries, and final
output. Fixed $V_e$ remains a zero-residual-covariance-RNG path.

The iteration order remains marker sweep, Dirichlet update, $V_b$ update,
optional $V_e$ update, unthinned convergence capture, then retained capture.

## Computational boundary

Complete enumeration is exponential. The implementation permits
$2\leq T\leq12$, but this guard does not imply that every marker count, chain
count, or retained-output request up to $T=12$ is feasible. Sampler work and
the mandatory markerwise pattern posterior both scale as $M2^T$.

Before native sampling, one R estimator calculates conservative peak
incremental memory for the fit. Its named breakdown includes pattern metadata
and workspaces, chain-private state, compact diagnostics, retained and
convergence capture, native results coexisting with Rcpp copies, raw-v2 bound
arrays and material binding temporaries, predictions, final state, provider
metadata, both contracted marker-by-pattern probability matrices, and
the packed-BED owner. For $N$ selected samples, that shared owner contributes

$$
M\times\operatorname{round\_up}\left(\left\lceil N/4\right\rceil,64\right)
$$

bytes. It is counted once, not once per chain. If rows are selected or
reordered, BED decoding also temporarily owns
$\left\lceil N_{\mathrm{source}}/4\right\rceil$ source-row bytes; the source
count is distinct from the selected count. The constructor, native guard, and R
preflight use the same checked 64-byte stride contract. The configured limit
applies to this estimated
incremental fit allocation, not total process RSS. The default qualification
limit is 256 MiB; `NULL` and positive `Inf` retain their resolved-spec meanings.
Failures report $M$, $T$, $K=2^T$, chain, retained, and convergence counts plus
the selected/source sample counts and largest named components, and occur
before provider construction or native sampling. For example, $M=10{,}000$
and $N=200{,}000$ require a 50,048-byte aligned stride and therefore
500,480,000 bytes for the packed owner alone.

The earlier mandatory dense pattern-transition diagnostic scaled as
$K^2=4^T$ and was not part of the posterior target. It is no longer allocated.
Diagnostics now retain one length-$K$ pattern-occupancy vector and one total
pattern-change count per chain. The historical `transition_counts` field is
present with value `NULL`; no zero matrix is fabricated. Compact diagnostics
use exact chain and declared activity-pattern IDs. Occupancy totals record one
event per marker update across burn-in and sampling sweeps, while each chain's
pattern-change count cannot exceed that event total. Raw-v2 validation enforces
this `compact_occupancy_v1` policy. Compact diagnostics are zero-RNG and do not
affect scheduling or scientific trajectories.

All byte products and sums are checked before arithmetic exceeds R's exact
integer range. Native code separately checks the dimensions of its immediate
allocations using checked `size_t` products. Thus the $T\leq12$ ceiling is a
computational implementation guard, not a model-theory limit or a promise that
large $M$, $T$, chain, and retained-draw combinations will fit in memory.

## Resolved specification and raw output

`prior$probability$activity_patterns` is the explicit binary source of truth;
its rows equal `model$state_space` and columns equal declared trait IDs. Raw
validation uses it for masks, PIPs, and the unique all-active probability
rather than decoding `1_1`.

All raw arrays retain dynamic named draw, chain, marker, trait,
activity-pattern, observation, trait-row, and trait-column axes. The
pleiotropic probability is the unique $(1,\ldots,1)$ column. Size-one draw and
chain axes are retained.

## Qualification evidence

- The complete Phase 4a/4b qualification suite remains green.
- The independent $T=3$ reference agrees for all eight probabilities, active
  means, active covariances, and representative Schur completions. Maximum
  errors were $4.17\times10^{-17}$ for probabilities and
  $5.56\times10^{-17}$ for conditional moments.
- $T=3$ fixed and sampled runs validate dynamic raw axes, $Y-XA$ statistics,
  and inverse-Wishart degrees of freedom and scales.
- $T=4$ fixed and sampled runs cover all 16 patterns, four-trait raw axes, SPD
  covariance draws, and serial/two-worker chain identity.
- Permuting traits, covariances, priors, and pattern masses produces the
  correspondingly permuted deterministic target.
- R/native packed-stride parity covers sample counts around byte and 64-byte
  alignment boundaries. The $M=10{,}000$, $N=200{,}000$ counterexample is
  rejected before provider construction and native entry, while a separate
  native-limit test rejects before packed-owner construction.
- Compact diagnostic adversarial tests reject missing or fabricated dense
  `transition_counts`, invalid chain/pattern axes, and invalid occupancy/change
  counts. Finite-limit and unlimited runs have identical scientific output and
  compact counts under the same seeds.
- Against clean Phase 4b, the $T=2$ deterministic conditional, states, pattern
  draws, seeds, and indices are exact. Dynamic trait loops change floating-
  point accumulation by at most $2.78\times10^{-17}$ in fixed mode and
  $4.44\times10^{-16}$ in sampled mode, with no RNG or state divergence. A
  duplicate two-trait kernel was not added solely for bitwise identity.

## Deferred work

Phase 5B owns public replacement and legacy MT migration. Summary-statistic MT
likelihoods, heterogeneous operators, overlap, missing phenotypes, fixed
effects, MT-BayesR/RC, regional covariance, templates, factors, and scalable
restricted-pattern methods remain deferred.
