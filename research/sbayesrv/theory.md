# SBayesRV version-1 research specification

Status: research contract and deterministic oracle. SBayesRV is the canonical
research name for the variance-annotation family descended from the maintained
BayesR-LV implementation. It is not a literature novelty claim, a new public
method, or authorization for `method = "sbayesrv"`.

SBayesRC and SBayesRV act on different prior quantities. SBayesRC changes
marker-specific component-membership probabilities through coefficients
$\alpha$. SBayesRV retains global BayesR component probabilities and changes
active-effect variance through coefficients $\theta$ and relative multipliers
$q_j$. The BayesR multiplier $\gamma_k$, global active-effect scale $v_b$, and
marker multiplier $q_j$ remain distinct.

## Annotation design and identification

Let $A$ be the raw marker-by-annotation matrix and let $X$ be its processed
version. Marker and annotation IDs are explicit and ordered. There is no
annotation intercept. A binary $0/1$ column is centered without scaling. A
continuous column is centered and divided by its sample standard deviation.
Every processed column therefore has mean zero. Non-finite, constant,
duplicate, or rank-deficient columns are rejected.

For marker $j$,

$$
\eta_j=X_j\theta,\qquad q_j=\exp(\eta_j).
$$

Because every column of $X$ is centered,

$$
\frac{1}{M}\sum_{j=1}^M\eta_j=0,
\qquad
\left(\prod_{j=1}^M q_j\right)^{1/M}=1.
$$

Thus $q$ describes relative marker variance and $v_b$ remains the global
scale. This normalization removes the direct intercept-like confounding
between $q$ and $v_b$; it does not make individual coefficients well
identified when annotations are strongly correlated. Numerical evaluation does
not clip $q$. Finite inputs whose exponential is outside the positive finite
double range are rejected explicitly.

## Version-1 prior

Let $c_j\in\{0,1,\ldots,K\}$ be the BayesR component, with global,
annotation-independent probabilities

$$
\Pr(c_j=k)=\pi_k,\qquad \sum_{k=0}^K\pi_k=1.
$$

The null component is exactly zero, $\gamma_0=0$ and $\beta_j=0$ when
$c_j=0$. For $k>0$, $\gamma_k>0$ and

$$
\beta_j\mid c_j=k,\theta,v_b
\sim N\!\left(0,v_b\gamma_kq_j\right).
$$

The proper coefficient prior is

$$
\theta\sim N(0,\sigma_\theta^2I),
$$

with maintained version-1 default $\sigma_\theta=0.7$. There is no annotation
intercept and no annotation-dependent $\pi_k$ in version 1.

## Hard-state conditional theta target

Given current effects and allocations, let
$\mathcal A=\{j:c_j>0\}$. Up to a constant independent of $\theta$, the
conditional log posterior is

$$
\ell_{\mathrm{cond}}(\theta)=
-\frac12\sum_{j\in\mathcal A}
\left\{
\eta_j+
\frac{\beta_j^2}
{v_b\gamma_{c_j}\exp(\eta_j)}
\right\}
-\frac{1}{2\sigma_\theta^2}\theta^\mathsf T\theta.
$$

Define

$$
r_j=\frac{\beta_j^2}
{v_b\gamma_{c_j}\exp(\eta_j)}.
$$

With $X_{\mathcal A}$ ordered as the active markers, the analytic gradient and
Hessian are

$$
\nabla_\theta\ell_{\mathrm{cond}}=
\frac12X_{\mathcal A}^\mathsf T(r-\mathbf 1)
-\frac{\theta}{\sigma_\theta^2},
$$

$$
\nabla_\theta^2\ell_{\mathrm{cond}}=
-\frac12X_{\mathcal A}^\mathsf T\operatorname{diag}(r)X_{\mathcal A}
-\frac{I}{\sigma_\theta^2}.
$$

The maintained sampler does not use this gradient. It applies elliptical slice
sampling (ESS) with the Gaussian prior as the ellipse and the active-state term
as the likelihood. If $\mathcal A$ is empty, the conditional posterior is the
Gaussian prior and production draws $\theta$ directly from it. ESS is an exact
MCMC transition for this conditional target, subject to its finite-state
guards; it is not an optimization or Laplace approximation.

## Independent-summary collapsed oracle

For a deterministic independent or negligible-LD summary fixture, write the
marker likelihood in score/precision units as

$$
\log L_j(\beta_j)=h_j\beta_j-\frac12d_j\beta_j^2+C_j,
$$

where $d_j\geq0$ and $C_j$ is common to all components. For active component
$k$, define

$$
\tau_{jk}^2=v_b\gamma_kq_j,
\qquad
D_{jk}=1+d_j\tau_{jk}^2.
$$

The active Gaussian conditional is

$$
\beta_j\mid h_j,d_j,c_j=k,\theta
\sim N(m_{jk},s_{jk}^2),
$$

$$
s_{jk}^2=\left(d_j+\tau_{jk}^{-2}\right)^{-1}
=\frac{\tau_{jk}^2}{D_{jk}},
\qquad
m_{jk}=s_{jk}^2h_j.
$$

After dropping $C_j$, the null and active log weights are

$$
w_{j0}=\log\pi_0,
$$

$$
w_{jk}=\log\pi_k
-\frac12\log D_{jk}
+\frac12\frac{h_j^2\tau_{jk}^2}{D_{jk}},\qquad k>0.
$$

Let

$$
a_j=\operatorname{logsumexp}_{k=0}^K w_{jk},
\qquad
\rho_{jk}=\exp(w_{jk}-a_j).
$$

Then $\rho_{jk}$ are normalized collapsed responsibilities, the collapsed
marker likelihood is $\exp(C_j+a_j)$, and the complete collapsed theta log
posterior, up to $\sum_jC_j$, is

$$
\ell_{\mathrm{coll}}(\theta)
=\sum_{j=1}^M a_j
-\frac{1}{2\sigma_\theta^2}\theta^\mathsf T\theta.
$$

For $k>0$,

$$
g_{jk}:=\frac{\partial w_{jk}}{\partial\eta_j}
=\frac12\left{
-\frac{d_j\tau_{jk}^2}{D_{jk}}
+\frac{h_j^2\tau_{jk}^2}{D_{jk}^2}
\right\},
$$

and $g_{j0}=0$. Therefore

$$
\nabla_\theta\ell_{\mathrm{coll}}=
X^\mathsf T
\left[
\begin{array}{c}
\sum_k\rho_{1k}g_{1k}\\
\vdots\\
\sum_k\rho_{Mk}g_{Mk}
\end{array}
\right]
-\frac{\theta}{\sigma_\theta^2}.
$$

Component weights, marker likelihoods, and responsibilities are evaluated with
log-sum-exp. The scientific target has no probability or variance clipping.
This factorized collapsed likelihood is an oracle only for independent or
negligible-LD summaries under the stated score/precision model. It is not an
exact appreciable-LD likelihood and is not an LD-block or retained-eigen
posterior.

## Reductions and research variants

- $\theta=0$, equivalently $q_j=1$ for every marker, gives ordinary SBayesR.
- Holding a supplied positive $q$ fixed gives the fixed marker-variance-
  multiplier model.
- One positive component gives the corresponding BayesC-style variance model.
- With no active markers, the hard-state conditional theta posterior is its
  Gaussian prior.
- Without a normalization convention, multiplying every $q_j$ by a constant
  is confounded with $v_b$.
- Correlated annotations can make raw $\theta$ weakly identified even when the
  induced $\eta$ and $q$ surface is stable.

The following are separate research designs:

1. **Joint learned SBayesRV:** update $\theta$ inside the genomic model.
2. **Frozen-q SBayesRV:** estimate $q$ externally or by cross-fitting and hold
   it fixed in the genomic fit.
3. **Oracle-q SBayesRV:** supply the data-generating $q$ for qualification.
4. **Probability-plus-variance models:** allow both $\alpha$ and $\theta$;
   deferred beyond version 1 because probability and variance signals require
   separate identification and mixing evidence.

Bayesian-MAGMA-style or other marker summaries may inform frozen or
cross-fitted $q$. Same-marker posterior evidence must not be recycled into its
own prior without an explicit design that prevents invalid data reuse.

## Maintained implementation crosswalk

The source code and tests are authoritative when prose differs.

| Theoretical quantity or transition | Current representation and source | Schedule and retained output | Qualification or discrepancy |
|---|---|---|---|
| Processed $X$ | `.stblr_preprocess_logvar_annotations()` in `R/stblr-logvar-annotations.R`, called by `R/stblr-csr-annot.R` and `R/stblr-block-eigen.R` | Constructed once before native dispatch; transform metadata becomes `annotation_transform` | Focused preprocessing and rejection tests in `tests/testthat/test-logvar-interface.R` |
| $\theta$ initialization and $N(0,\sigma_\theta^2I)$ prior | `.stblr_logvar_theta_init()` and `theta_prior_sd`; native constant/default and draws in `src/st_logvar_annotation_prior.h` | Initialized before chain start; one vector per trait/chain | Default 0.7, finite positive validation, empty-active prior-draw test |
| $\eta=X\theta$ and $q=\exp(\eta)$ | `calculate_eta()` and `calculate_prior_scale_from_eta()` in `src/st_logvar_annotation_prior.h` | Rebuilt at initialization and immediately after each theta update | Deterministic native/R checks; overflow guard at finite log scale, no clamp |
| Global component probabilities $\pi$ | Common BayesR state in `src/blr_csr_bayesr_core_impl.h` | Global Dirichlet update after the residual update in executable code; `pi_trace`, `pi_mean`, `pi_final` | Probabilities do not depend on annotations. Older implementation-contract prose places this update earlier; that is an update-order documentation discrepancy, not a target change |
| Component/effect update using $v_b\gamma_kq_j$ | `sampleBetaR_ST_csr()` in `src/st_cpg_omp_csr_bayesr.cpp`, called by the common engine with the LV policy's `prior_scale()` | Marker/component sweep first each iteration | Theta-zero and fixed-q trajectory reductions in `tests/testthat/test-logvar-bayesr.R` |
| Scale-aware $v_b$ | `sampleB_bayesr_ST_csr()` in `src/st_cpg_omp_csr_bayesr.cpp` | After marker/optional LD-swap updates and before theta; uses $\sum_{j:c_j>0}\beta_j^2/(\gamma_{c_j}q_j)$ | Deterministic fixed-q trajectory and learned-oracle evidence |
| Hard-state theta target and transition | `theta_log_likelihood_bayesr()` and `elliptical_slice_update()` in `src/st_logvar_annotation_prior.h`; attached by `CsrLogvarBayesRPolicy` in `src/blr_csr_logvar_bayesr_core_impl.h` | Post-$v_b$ hook; exact ESS using current hard effects/components, then rebuild $q$ | Production does **not** use collapsed responsibilities and does not use the analytic gradient |
| CSR likelihood | `src/st_cpg_omp_csr_logvar_bayesr.cpp` plus common CSR BayesR engine | Optional CSR LD-swap; scalar residual update follows theta | Full conditional target uses the declared CSR operator, not the factorized collapsed oracle |
| Retained block-eigen likelihood | `R/stblr-logvar-block-eigen.R`, `src/st_cpg_omp_block_eigen_logvar_bayesr.cpp`, and the same common engine/policy | Same prior and theta transition; operator and residual policy are block-eigen specific | Effectively exact fixtures compare CSR/eigen; retained rank remains an approximate declared likelihood |
| Native raw and formatted theta/q | Raw-v1 `annotation` namespace is decorated by `src/st_logvar_annotation_rcpp.h`; formatting in `R/stblr-logvar-annotations.R` | `theta` is a retained-sample posterior mean; `theta_trace` is retained when requested; `marker_prior_scale` is posterior mean $q$; `annotation_variance_ratio=\exp(\mathbb E\theta)$ | The variance ratio is not $\mathbb E\{\exp(\theta)\}$. Marker-by-iteration q histories are not retained |
| Component/effect/$v_b$ outputs | Common raw-v1 and formatter paths | Posterior component assignment probabilities, posterior mean effects, final state, and variance traces follow the common BayesR contract | These are posterior genomic quantities, not substitutes for $\theta$ or prior $q$ |
| Raw-v2 bridge | `R/blr-raw-v2.R` recognizes `annotation_log_variance` and keeps native diagnostics | Eligible one-chain fits can receive raw-v2 common BayesR quantities; multi-chain v1 lacks every final chain state | Current raw-v2 has no dedicated theta/q draw namespace; authoritative research cross-checks use raw-v1 annotation fields and formatted fields without changing encapsulation |

### Explicit audit answers

1. The maintained BayesR-LV prior target is mathematically the version-1
   SBayesRV target above.
2. Production learns theta jointly from current hard component and effect
   states, not collapsed responsibilities.
3. The theta transition is Gaussian-prior elliptical slice sampling, an exact
   conditional MCMC update; it is not an approximate gradient update.
4. Centered $X$ fixes geometric-mean $q$ at one, separating relative $q$ from
   global $v_b$. Rank checks remove exact design nonidentifiability, while
   correlated columns can still be weakly identified.
5. Theta draws are the sampled annotation parameters; q is their deterministic
   marker-level transform. `exp(theta)` is a conditional coefficient-scale
   variance ratio. Global component probabilities, effects, allocations, and
   $v_b$ are separate sampled genomic quantities.
6. Theta traces are retained only under the requested chain-output contract;
   q is retained as posterior marker means, not a marker-by-draw history.
7. CSR and retained block-eigen use the same prior, scale-aware $v_b$, and theta
   transition relative to different declared likelihood operators.
8. General model limitations include relative-scale identification, weak
   coefficient identification under correlated annotations, fixed version-1
   theta prior structure, and operator-relative likelihoods. Study 12's finite
   scenarios, short chains, weak q-amplitude recovery, and lack of downstream
   gain are evidence limitations of that study, not mathematical definitions
   of SBayesRV.
