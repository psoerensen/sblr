# BLR prior-variance calibration audit

> **Status and provenance.** This is a historical before-change audit anchored
> at commit `bd8e2c8148a0d9540dc20716455706beeb16fa86`, followed by the
> calibration contract selected and subsequently implemented. The
> “Scalar route audit” and “Multivariate route audit” tables record defects at
> that audited commit; they are not claims about current package behavior or
> current route support. The final “Calibration contract selected after audit”
> section records the current calibration design. For current annotation-prior
> capabilities and future calibration work, use
> [the annotation-prior architecture audit](annotation_prior_architecture_audit.md),
> [capability matrix](annotation_prior_architecture_matrix.md), and
> [implementation plan](annotation_prior_architecture_implementation_plan.md).

## Scope and native effect-scale contract

This audit was performed from commit `bd8e2c8148a0d9540dc20716455706beeb16fa86`
before changing prior construction. It covers the public scalar and multivariate
BLR routes and the R helpers that construct their native variance arguments.

The native scalar kernels use

```text
Var(b_jt | state k) = B_t * gamma_k * q_jt * v_jt
```

for active states. BayesC has `gamma = 1`; ordinary global models have
`q = v = 1`. The CSR MAF-S kernels use
`q_jt = {2 f_j (1-f_j)}^(maf_effect_s + 1)`. Fixed-marker, group, and learned
BayesC kernels multiply the active-effect variance by their resolved marker or
group variance multiplier `v_jt`. BayesR and BayesRC kernels multiply it by the
selected component multiplier. These multipliers also divide the sufficient
statistic used by the native `B` update, so `B` is the unmultiplied base
marker-effect variance. No output rescaling changes this contract.

The native MT BayesR/BayesRC kernel applies a joint state's scalar
`gamma_s * q_j` to the marker covariance and restricts it to the state's binary
trait pattern `z_s`. BayesC is the same contract with `gamma_s = q_j = 1`.
Consequently,

```text
W_tu = sum_j sum_s p_js gamma_s q_j z_st z_su
E[G]_tu = B_tu * W_tu
```

where the last multiplication is elementwise. This is the authoritative native
effect-scale contract used below.

Phenotype variance is `colSums(y^2)/(n-1)` for centered packed-BED phenotypes
and `stats$yy/(n-1)` for summary-statistic routes. The target is
`vy_t * h2_t`. Residual defaults remain `vy_t * (1-h2_t)`.

## Scalar route audit

| Public route | Internal prior route | Initial/prior probabilities | Native multipliers | Formula at audited commit | Does `h2` equal expected prior genetic variance? |
|---|---|---|---|---|---|
| `stblr_bed(method="bayesc")` packed BED | `.stblr_bed_impl()` -> `.make_stblr_priors()` | `pi_vb_init`; separately `pi_prior_mean` from Beta controls | none | `B=vy*h2/(m*pi_vb_init)`; `ssb=((nub-2)/nub)*vy*h2/(m*pi_prior_mean)` | Yes, global special case |
| `stblr_csr(method="sbayesc")` CSR | `stblr_csr()` -> inline prior block | same separate BayesC probabilities | optional fixed or sampled MAF-S scale | divides only by `m*pi`; ignores `q` | Only without MAF-S |
| `stblr_block_eigen(method="sbayesc")` | block-eigen wrapper -> inline prior block | same | optional fixed or sampled MAF-S scale | divides only by `m*pi`; ignores `q` | Only without MAF-S |
| `stblr_bed(method="bayesr")` | `.stblr_bed_impl()` -> `.make_stblr_bayesr_priors()` | initial `pi`; Dirichlet `alpha` also accepted | `mixture_var` | both `B` and `ssb` divide by `m*sum(pi*gamma)` | `B`: yes; `ssb`: only if normalized `alpha` has mean `pi` |
| `stblr_csr(method="sbayesr")` | CSR BayesR dispatch -> shared BayesR helper | initial `pi`; custom `alpha` | component and optional MAF-S | both use initial mixture weight; ignore `alpha` mean and `q` | Only for `B` without MAF-S and for default matching `alpha` |
| `stblr_block_eigen(method="sbayesr")` | block-eigen BayesR -> same helper | as CSR | component and optional MAF-S | same as CSR | Same defect |
| `stblr_bed(method="bayesrc")` | `.stblr_bed_impl()` -> annotation processing -> `.make_stblr_priors()` | marker probabilities resolved by processed `A` and initial annotation coefficients | component multipliers | divides by total active probability only | No; component multipliers and marker heterogeneity omitted |
| `stblr_csr_annot(method="sbayesrc")` | `stblr_csr_sbayesrc_generic()` -> `.stblr_make_csr_variance_priors()` before annotation coefficients | global architecture probabilities initially; actual marker probabilities resolved later | component and optional MAF-S | divides by scalar active probability | No |
| block-eigen `sbayesrc` | same generic wrapper with block-eigen native function | same | same | same | No |
| fixed-marker annotation BayesC | `stblr_csr_prior_annot()` -> variance helper **before** resolving marker priors | resolved `fixed_pi_marker` or annotation-derived `pi_marker`; Beta mean remains separately available | resolved `fixed_vb_multiplier` or annotation multiplier | divides by scalar `m*pi_vb_init` / `m*pi_prior_mean` | No for heterogeneous marker priors or multipliers |
| group BayesC | `stblr_csr_group_annot()` -> variance helper before group expansion | initial group probabilities, with group Beta priors | initial group variance multipliers | divides by scalar global probabilities | No for heterogeneous groups |
| learned-annotation BayesC | `stblr_csr_learn_annot()` -> variance helper before initial annotation transform | marker probabilities produced from initial annotation coefficients | annotation-derived variance multipliers | divides by scalar global probabilities | No for heterogeneous annotations |

The scalar BayesRC annotation link is probit stick breaking. Processed annotation
rows `A_j` and the actual initialized coefficient matrix determine marker
probabilities exactly via `sbayesrc_marker_pi(A, alpha_init, gamma)`. There is no
closed-form prior expectation over subsequently sampled annotation coefficients.
The safe prior-scale policy is therefore the named
`resolved_initial_annotation_probabilities` policy; it is exact for the initial
state and deliberately does not claim to integrate over annotation coefficients.

For sampled MAF-S, native initialization uses `maf_effect_s_init`; the audited R
construction instead behaves as if all `q_j=1`. Calibration must use
`q_j={2f_j(1-f_j)}^(maf_effect_s_init+1)` once, and must not renormalize `B`
during MCMC.

## Multivariate route audit

| Public route/operator | Internal prior route | Probability source | Multipliers | Formula at audited commit | Assessment |
|---|---|---|---|---|---|
| MT BayesC, packed BED | `mtblr_bed()` -> `.mtblr_models()` -> local defaults | joint trait-pattern `p_s` | none | all traits divide by `m*(1-p_null)` | Incorrect when trait marginal or pairwise pattern weights differ |
| MT BayesC, CSR | `mtblr_csr()` -> same pattern helper | same | none | same total-active divisor | Same defect |
| MT BayesC, block eigen | `mtblr_block_eigen()` -> same | same | none | same | Same defect |
| MT BayesR/SBayesR, all three operators | `.mtblr_bayesr_spec()` expands pattern x component states; wrappers construct defaults | initial `joint_pi`; prior mean `joint_pi_prior/sum` | joint component `gamma_s`; fixed MAF-S `q_j` | wrappers still divide by total non-null probability | Omits trait patterns, components, MAF scale, and prior-mean distinction |
| MT BayesRC/SBayesRC, all three operators | BayesR state expansion plus `.mtblr_bayesrc_controls()` | marker-state probabilities from processed annotations and initialized component coefficients; conditional pattern probabilities initialized separately | component and optional fixed MAF-S | wrappers divide by total non-null global probability | Omits marker annotation probabilities, patterns, components, and MAF scale |

MT BayesRC uses annotation-dependent component probabilities and conditional
active trait-pattern probabilities. The resolved initial marker-state
probability is their product (with the null component mapped to the null state).
As in scalar BayesRC, the prior-mean annotation policy must be explicitly named
`resolved_initial_annotation_probabilities`; no unavailable analytic expectation
over annotation coefficients is invented.

The automatic target supplied by `h2` is diagonal. It is uniquely and safely
calibrated by `B_tt=vy_t*h2_t/W_tt`. Explicit user `vb` remains authoritative.
An automatic full marker covariance from a user-supplied full genetic covariance
would require elementwise division by `W`; this need not preserve positive
semidefiniteness. Therefore no such automatic conversion is safe: a full `vg`
with omitted `vb` must receive a clear validation error and an explicit `vb`
must be supplied. Explicit `vb`, `ssb_prior`, and residual priors are validated
but not recalibrated.

## Calibration contract selected after audit

For scalar trait `t`, resolve marker/component probabilities and every native
variance multiplier first, then compute

```text
weight_initial_t = sum_j sum_k p_initial_jkt gamma_k q_initial_jt v_jt
weight_prior_t   = sum_j sum_k p_prior_jkt   gamma_k q_initial_jt v_jt
B_t              = vy_t h2_t / weight_initial_t
ssb_prior_t      = ((nub-2)/nub) vy_t h2_t / weight_prior_t
```

Initial probabilities and Beta/Dirichlet prior means are distinct inputs.
Global BayesC reduces exactly to `m*pi`; global BayesR reduces exactly to
`m*sum_k(pi_k gamma_k)`. For annotation models without an analytic coefficient
prior expectation, initial resolved probabilities are also the documented
prior-mean policy. Fixed and sampled-initial MAF-S use the resolved initial
scale, with no dynamic recalibration.

For MT models, construct separate initial and prior-mean matrices `W` using the
native state patterns and multipliers. Default diagonal `vb` and `ssb_prior`
use `diag(target)/diag(W_initial)` and
`((nub-2)/nub)*diag(target)/diag(W_prior)`, respectively. Matrix-valued metadata
remain matrices; scalar mixture-weight labels are not used for MT fits.
