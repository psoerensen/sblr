# Methods-to-source traceability audit

**Status:** current audit record

**Audit date:** 2026-08-12

**Theory checkpoint:** `44f84f21b6a84b6d2986d81509e93d5c6857a41d`

**Theory checkpoint subject:** `Finalize canonical statistical methods documentation`

**Audited branch:** `master`

## 1. Executive summary

This audit compared the approved statistical specification in
`docs/methods` with the executable public R dispatch, R and C++ source,
generated interfaces, output formatting, permanent tests, and maintained
developer contracts at the theory checkpoint above. The repository was clean
when the audit began. No implementation, Methods, test, contract, or generated
file was changed.

The principal single-trait BayesC and BayesR marker conditionals agree with the
approved theory. The exact null state, active-state conditional mean and
variance, normalized state probabilities, residual-score update, weighted
global-scale update, MAF-S multiplier, and posterior marker summaries were
traced through the BED, CSR, and retained block-eigen routes. The SBayesRC
continuation-probit construction, annotation-coefficient conditional, and
scaled-inverse-chi-squared hierarchy also agree algebraically with the current
relative-scale parameterization. The retained block-eigen construction obeys
the approved `Q_b` and `w_b` identities.

Two confirmed implementation errors require correction before every reported
quantity can be treated as a draw from the documented model:

1. Multi-trait marker updates use a latent inverse-Wishart draw of the shared
   marker-effect covariance, but each iteration later replaces that draw with
   a heuristic diagonal/correlation construction and stores the replacement
   as `V_b`. The replacement is not the conditional distribution stated by
   the latent inverse-Wishart code. It does not appear to alter the subsequent
   marker-effect transition because the proper latent draw overwrites it before
   the next marker block, but the stored/final `V_b` and its traces are not
   posterior draws from the documented covariance model.
2. The learned-logistic single-trait CSR BayesC route defines marker-specific
   inclusion probabilities through an annotation offset from a global
   `pi`, then updates that global `pi` with the ordinary conjugate Beta count
   update. Once annotation offsets are nonzero, the Beta update is not the
   full conditional. This changes the posterior target when `updatePi = TRUE`.

There are also user-facing output-definition risks. CSR SBayesRC reports
posterior means in fields named `alpha_final` and `sigmaSqAlpha_final`;
potentially indefinite CSR quadratic forms are exposed under unqualified
variance/SSE names; and multi-trait summary routes expose a descriptive
cross-operator bilinear form as `cov_g_*` even when a paired common-sample
genotype covariance is unavailable. The canonical formatted fit deliberately
uses `pi_trace`, not `pis`; `pi_trace` is a trace of prior state mass and is
not a PIP.

The audit found one confirmed posterior-target error: the learned-logistic
global-`pi` update. No confirmed posterior-target error was found in the
ordinary single-trait BayesC, BayesR, or SBayesRC marker kernels, the ordinary
multi-trait joint-state marker kernel, the current multi-trait `P x H` state
factorization, or the retained block-eigen likelihood transformation.

### Finding counts

Each item in the prioritized table has one primary classification.

| Primary classification | Count |
|---|---:|
| `CONFIRMED_IMPLEMENTATION_ERROR` | 2 |
| `LIKELY_IMPLEMENTATION_ERROR_REQUIRING_TEST` | 0 |
| `EQUIVALENT_PARAMETERIZATION` | 1 |
| `INTENTIONAL_NUMERICAL_APPROXIMATION` | 1 |
| `IMPLEMENTATION_SPECIFIC_SAFEGUARD` | 1 |
| `METHODS_DOCUMENTATION_MISMATCH` | 0 |
| `DEVELOPER_CONTRACT_MISMATCH` | 2 |
| `OUTPUT_SCHEMA_MISMATCH` | 4 |
| `MISSING_OR_INADEQUATE_TEST` | 3 |
| `UNSUPPORTED_INTERNAL_CODE` | 1 |
| `UNRESOLVED_SCIENTIFIC_DECISION` | 3 |
| `NO_ACTION` | 1 |
| **Total** | **19** |

| Severity | Count |
|---|---:|
| Critical | 0 |
| High | 7 |
| Moderate | 9 |
| Low | 2 |
| Informational | 1 |

## 2. Audited checkpoint and repository state

The precondition commands were:

```text
git status --short
git log -1 --oneline
git diff --cached --name-only
git branch --show-current
git rev-parse HEAD
```

They established:

- branch `master`;
- current and theory-checkpoint commit
  `44f84f21b6a84b6d2986d81509e93d5c6857a41d`;
- no modified, staged, or untracked files;
- no pre-existing user change to incorporate into the audit state.

The complete root `AGENTS.md` was read. A repository-wide search found no
additional `AGENTS.md` applying to `R`, `src`, `tests`, `tools`, `docs/dev`, or
`docs/methods`. The authority classifications in `docs/dev/README.md` were
applied. Historical, superseded, decision-record, and research-experimental
documents were not used as evidence of current public support.

All compilation and mutation-capable validation was run in a private copy made
from `git archive` of the exact checkpoint. It contained the local root
`AGENTS.md` and no working-tree overlay because the original tree was clean.

## 3. Authority and theory reference

The audit read the package metadata and build/documentation configuration, all
12 Methods QMD files, all maintained current contracts listed in the phase-2
request, the relevant annotation and block-eigen qualification records, the
public R interfaces and formatters, their native implementation paths, and the
permanent tests identified by the maintained test-ownership contract.

The statistical reference was the approved Methods checkpoint. In particular:

- BayesC has one exact null state and an active Normal prior with variance
  `v_b q_j`.
- BayesR has one exact null component and positive components with variance
  `gamma_k q_j v_b`.
- ordinary multi-trait BayesR uses one full joint pattern-by-component simplex;
  it is not a `P x H` model.
- current multi-trait BayesRC uses marker-specific component probabilities `P`
  times global conditional non-null sharing probabilities `H`.
- fixed MAF-S on standardized effects uses `q_j(S) = h_j^(S+1)`.
- exact genotype variances require an actual Gram matrix or compatible
  positive-semidefinite operator; indefinite sparse-operator expressions are
  algebraic operator-relative quadratics.
- cumulative marginal-PIP sets are not joint causal-configuration posterior
  probabilities.

Named external-method cross-checks used the following primary sources:

- Habier et al. (2011), BayesC, BMC Bioinformatics,
  <https://doi.org/10.1186/1471-2105-12-186>.
- Erbe et al. (2012), BayesR, Journal of Dairy Science,
  <https://doi.org/10.3168/jds.2011-5019>.
- Zeng et al. (2018), BayesS, Nature Genetics,
  <https://doi.org/10.1038/s41588-018-0101-4>.
- Lloyd-Jones et al. (2019), SBayesR, Nature Communications,
  <https://doi.org/10.1038/s41467-019-12653-0>.
- MacLeod et al. (2016), BayesRC, BMC Genomics,
  <https://doi.org/10.1186/s12864-016-2443-6>.
- Zheng et al. (2024), SBayesRC, Nature Genetics,
  <https://doi.org/10.1038/s41588-024-01704-y>, together with the authors'
  official reference repository, <https://github.com/zhilizheng/SBayesRC>.

## 4. Public model and interface inventory

Support was determined from public validation and dispatch, executable native
code, coherent output formatting, documentation, and permanent tests. The
standalone `.blr_model_capability_matrix()` was not used as a support oracle.

| Public interface | Public models | Data/operator | Native boundary and owner | Formatter/schema | Explicit restrictions |
|---|---|---|---|---|---|
| `stblr_bed()` | `bayesc`, `bayesr`, `bayesrc` | ST individual-level packed BED | scheduled BED BayesC entry; `stblr_cpg_omp_bed_marker_scheduled_chains_bayesr()`; `stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc()`; cores in `blr_bed_*_core_impl.h` | `.as_stblr_fit()` then `.blr_finalize_fit()` | no MAF-S; no learned log-variance route |
| `stblr_csr()` | `sbayesc`, `sbayesr` | ST summary-statistic CSR | `st_cpg_omp_csr.cpp`; `st_cpg_omp_csr_bayesr.cpp`; `run_csr_bayesc_engine()` and `run_csr_bayesr_engine()` | named `stblr_raw`, `.as_stblr_fit()`, `.blr_finalize_fit()` | sparse operator need not be PSD |
| `stblr_csr_annot()` | BayesC probability providers | ST CSR | learned/fixed/group policy paths in `st_cpg_omp_csr_annot.cpp` and `blr_csr_learned_annotation_bayesc_core_impl.h` | canonical ST formatter | not a general BayesR/BayesRC provider interface |
| `stblr_csr()` SBayesRC route | `sbayesrc` | ST CSR | `st_sbayesrc_omp_csr.cpp` and `st_bayesrc_annotation_prior.h` | canonical ST formatter plus annotation fields | current relative grid and priors differ from published defaults |
| `stblr_block_eigen()` | `sbayesc`, `sbayesr`, `sbayesrc` | ST retained low-rank or reconstructed dense block operator | `st_block_low_rank.cpp`, `blr_block_low_rank.h`, then the ST CSR-family policies | canonical ST formatter with block diagnostics | low-rank residual policy is model/argument dependent |
| `mtblr_bed()` | `bayesc`, `bayesr`, `bayesrc` | MT common-individual BED | `mtblr_bed_chains_internal()` and `blr_mt_bed_core_impl.h` | `mtblr_legacy_to_raw()` then MT formatters | residual covariance may be full or diagonal |
| `mtblr_csr()` | `sbayesc`, `sbayesr`, `sbayesrc` | MT summary CSR | `mtblr_csr_chains_raw_internal()` and `blr_mt_default_core_impl.h` | named MT raw then `.mtblr_bayesr_format_fit()` / `.mtblr_bayesrc_format_fit()` | no general sample-overlap likelihood; trait operators may differ |
| `mtblr_block_eigen()` | `sbayesc`, `sbayesr`, `sbayesrc` | MT reconstructed-dense block operator | `mtblr_block_eigen_chains_raw_internal()` and `run_mt_block_eigen_core()` | named MT raw and MT formatters | diagonal residual treatment; not the ST retained `Q/w` runtime |

Single-trait BayesC-V/SBayesC-V and BayesR-V/SBayesR-V are publicly available
only for CSR and retained block-eigen routes. BED-LV and all MT-LV routes are
rejected. Fixed marker probabilities, group probabilities, and learned
logistic probabilities are public only on the single-trait CSR BayesC provider
route. Fixed MAF-S is supported on the ST CSR/retained-block BayesC, BayesR,
and SBayesRC families and on the applicable MT BayesR/BayesRC families; sampled
`S` is ST only. `maf_effect_s` is independent of the summary-model `s` prefix.

## 5. Native implementation ownership

The principal ownership map is:

| Operation | Source owner and symbol |
|---|---|
| CSR algebra and residual-score column updates | `src/st_ld_operator.h`: `CsrOperator` |
| ST BayesC CSR transition | `src/blr_csr_bayesc_core_impl.h`: `run_csr_bayesc_engine()`, `run_csr_bayesc()` |
| ST BayesR CSR transition | `src/blr_csr_bayesr_core_impl.h`: `run_csr_bayesr_engine()`, `run_csr_bayesr()` |
| ST BED BayesC | `src/blr_bed_scheduled_bayesc_core_impl.h`: `run_bed_scheduled_bayesc_chain()` |
| ST BED BayesR | `src/blr_bed_bayesr_core_impl.h`: `run_bed_bayesr_chain()` |
| ST BED BayesRC | `src/blr_bed_bayesrc_core_impl.h`: `run_bed_bayesrc_chain()` |
| SBayesRC stick prior | `src/st_bayesrc_annotation_prior.h`: probability and hierarchy helpers; `src/st_sbayesrc_omp_csr.cpp`: sampler/binding |
| Learned logistic inclusion | `src/st_cpg_omp_csr_annot.cpp`: `make_pi_from_annotation()`, `logpost_eta_pi()`, `samplePi_ST_annot()`; `src/blr_csr_learned_annotation_bayesc_core_impl.h` |
| Learned log-variance | `src/st_logvar_annotation_prior.h`; `src/st_cpg_omp_csr_logvar_bayesc.cpp`; `src/st_cpg_omp_csr_bayesr_logvar.cpp` |
| Retained block construction/runtime | `src/st_block_low_rank.cpp`; `src/blr_block_low_rank.h` |
| MT state and covariance helpers | `src/mtblr.cpp`: `sampleB_latent()`, `sampleBset()`, `sampleB()`, `computeG()` and joint-state helpers |
| MT CSR/block sampler | `src/blr_mt_default_core_impl.h`; raw adapters in `src/mtblr.cpp` |
| MT BED sampler | `src/blr_mt_bed_core_impl.h`; chain aggregation in `src/blr_mt_bed_chains_*` |
| ST raw schema/formatting | `R/sparse_ld_bed_helper.R`: `.as_stblr_fit()`; `R/blr-unified.R`: `.blr_finalize_fit()` |
| MT formatting | `R/mtblr-bayesr.R`; `R/mtblr-bayesrc.R`; route-specific R wrappers |
| Fine mapping and credible sets | `R/credible_sets.R`; `R/finemap-stblr-csr.R`; `R/extract-stblr-finemap-loci.R` |

Legacy positional adapters remain in source as conversion boundaries, but the
public paths return named raw objects and are formatted through the canonical
raw-to-fit layer. No silent public positional fallback was found.

## 6. Theory-to-source traceability matrix

`C_jj` below is the marker diagonal on the operator's cross-product scale and
`s_j` is the current partial-residual score.

| Theoretical operation | Methods target | Public/native path | Implemented operation and output | Status / consequence |
|---|---|---|---|---|
| BayesC active conditional | `V_j = v_e/(C_jj + v_e/(v_b q_j))`, `mu_j = V_j s_j/v_e` | BED, CSR, retained block; BayesC cores | algebraically identical precision, mean, and Normal draw | `EXACT_MATCH` |
| BayesC state odds | active prior mass `pi_j`; exact null | same | log determinant and score terms contain the required halves; log-scale normalization is stable | `EXACT_MATCH` |
| BayesC residual score | `r <- r - C_.j Delta b_j` | CSR/block policies | sparse column update or `Q`-space update | `EXACT_MATCH` relative to fitted operator |
| BayesC `v_b` | weighted active sum `sum b_j^2/q_j` | CSR and BED policies | inverse-chi-squared scale uses active count and weighted sum | `EXACT_MATCH` |
| BayesR component conditional | exact zero plus `gamma_k q_j v_b` | BED/CSR/block BayesR | correct component precision, log weight, normalization, state draw, and effect draw | `EXACT_MATCH` |
| BayesR global scale | `sum b_j^2/(gamma_cj q_j)` | BayesR policies | same weighted quadratic over non-null states | `EXACT_MATCH` |
| BayesR PIP | `1 - Pr(c_j=0 | D)` | component counts to formatter | non-null retained-state frequency; component probabilities retained separately | `EXACT_MATCH` |
| SBayesRC probabilities | continuation probit with final residual stick | `st_bayesrc_compute_snp_pi()` | probabilities are nonnegative, floored where configured, and renormalized to one | `EXACT_MATCH` plus safeguard |
| SBayesRC latent sticks | eligible markers at stick `k`; success/failure truncation directions | annotation update helpers | correct eligibility and truncated-Normal directions; empty sticks retain the prior contribution | `EXACT_MATCH` |
| SBayesRC alpha variance | `(sum alpha^2+b)/chi^2_(r+a)` | CSR/MT/individual BayesRC hierarchy | same convention; current default `(a,b)=(2,2)` | `EXACT_MATCH`; not published `(4,4)` |
| Published SBayesRC grid | `(0,10^-5,10^-4,10^-3,10^-2) v_g^SNP` | current relative ladder | `(0,.001,.01,.1,1) v_b` | `EQUIVALENT_PARAMETERIZATION` only if `v_b=.01 v_g^SNP`; current updating does not enforce this absolute equality |
| Learned logistic `eta` | Bernoulli likelihood under marker-specific logit probabilities | ST CSR annotation BayesC | MH log posterior includes marker Bernoulli terms and coefficient prior | `EXACT_MATCH` conditional on global `pi` and clipping definition |
| Learned logistic global `pi` | nonconjugate conditional under offset logits | same | ordinary Beta active/null count update | `CONFIRMED_IMPLEMENTATION_ERROR`; posterior target changes |
| Log-variance multiplier | `q_j=exp(z_j^T theta)`, identifiable relative scale | ST CSR/block BayesC-V/R-V | centered predictor, no intercept, geometric-mean-one normalization, exponent guard | `EXACT_MATCH` plus numerical safeguard |
| Log-variance global scale | `sum b_j^2/(gamma_cj q_j)` as applicable | LV policies | correct weighted quadratic | `EXACT_MATCH` |
| MAF-S | `q_j=h_j^(S+1)` | claimed ST and fixed-MT paths | exact exponent; `S=-1` yields unit multiplier | `EXACT_MATCH` |
| Ordinary MT state prior | one full joint state simplex, unique null | MT BayesR paths | joint pattern/component state list and Dirichlet update | `EXACT_MATCH`; not `P x H` |
| MT BayesRC state prior | marker `P` times conditional non-null global `H` | MT BayesRC paths | marker component prior multiplied by pattern prior then normalized | `EXACT_MATCH` |
| MT marker conditional | active subspace covariance `gamma_k q_j V_b,p` | `mt_joint_marker_kernel()` / BED analogue | masked effective effects with latent full effect vector | `EXACT_MATCH` for declared diagonal/common-operator likelihood conditions |
| MT covariance used by markers | inverse-Wishart latent augmentation | `sampleB_latent()` before marker blocks | proper inverse-Wishart draw from latent base effects | `EXACT_MATCH` under the latent augmentation |
| MT covariance stored/reported | posterior draw of shared `V_b` | `sampleB()` after marker sweep | heuristic marginal variance draws plus shrunk/clamped empirical correlations and eigenvalue repair | `CONFIRMED_IMPLEMENTATION_ERROR`; reported `V_b` is not the covariance draw used by the marker kernel |
| Exact BED genomic covariance | `B_effect^T X^T X B_effect/N` | ST/MT BED | residual/fitted-value identities on common individuals | `EXACT_MATCH` |
| CSR genomic quadratic | operator-relative `b^T C_tilde b/N` unless PSD | `CsrOperator`, `computeG()` | algebraic quadratic stored under variance/covariance names | calculation matches operator; label/interpretation does not |
| Retained block operator | `C_tilde=Q^T Q`, `s_tilde=Q^T w` | low-rank block route | stored orientation is dimensionally equivalent; identities hold | `EXACT_MATCH` up to retained-rank/float approximation |
| Residual policies | `gctb_block`, `fixed_block`, or global projected contract | low-rank policy dispatch | source dispatch agrees with maintained contract | `EXACT_MATCH`; not complete GCTB reproduction |
| `v_LE`, `v_LD` | diagonal contribution and algebraic remainder | ST/MT paths | `v_LD=v_g-v_LE`; negative values are permitted | `EXACT_MATCH`; names need operator qualification for indefinite CSR |
| Posterior mean effect | `E(b_j|D)` | native accumulated `bm` | accumulated effective effect divided by retained draws | `EXACT_MATCH` |
| Draw-level `v_g` | `E[b^T C b/N|D]` estimated draw by draw | `vgs` traces | computed each iteration, not from posterior mean effects | `EXACT_MATCH` |
| Cumulative-PIP set | marginal-PIP ranking/set construction | credible-set helpers | marginal-PIP cumulative sets and approximate LD-conditioned extensions | `EXACT_MATCH` interpretation; not joint causal configurations |

## 7. Conditional-update audit

### 7.1 BayesC

For an active marker with prior variance `tau_j = v_b q_j`, the code's log
Bayes factor is algebraically

```text
0.5 log(v_e / (v_e + C_jj tau_j))
+ 0.5 s_j^2 tau_j / {v_e (v_e + C_jj tau_j)}.
```

This equals the integrated Normal likelihood ratio. The active conditional has
precision `C_jj/v_e + 1/tau_j`, mean `s_j/(C_jj+v_e/tau_j)`, and variance
`v_e/(C_jj+v_e/tau_j)`. The code uses active probability, not null
probability, for `pi`. Residual changes have the correct sign. Beta updating of
ordinary BayesC `pi` uses active and null counts correctly. No-active and
all-active edges are protected without changing the ordinary interior
conditional.

### 7.2 BayesR

The null component is state zero. Positive state `k` uses
`gamma_k q_j v_b`; the determinant, quadratic, and prior-mass terms are
present before log-sum-exp normalization. C++ state zero maps to R component
label one only at the representation boundary. PIP is the posterior frequency
of all nonzero components. Mixture-probability updates use full component
counts. The global-scale update removes both `gamma_k` and `q_j`, as required.

### 7.3 BayesRC and SBayesRC

The continuation-probit probabilities sum to one after construction. The last
component receives the remaining stick. Latent variables above and below zero
use the correct truncation directions, and only markers eligible at a stick
enter that stick's regression. Intercept and non-intercept priors are separated.

The current relative component ladder is internally coherent with the current
`v_b` update, but it is not automatically the published absolute SBayesRC
prior. Exact equivalence to the published grid requires
`v_b = 0.01 v_g^SNP`. Current initialization calibrates a scale and subsequent
posterior updating treats it as unknown; it does not preserve that fixed
absolute identity. Current `(2,2)` alpha-variance hyperparameters also differ
from the published `(4,4)` values under the displayed convention.

## 8. Annotation architecture audit

### Probability architecture `P`

- Fixed marker and group probabilities on the ST CSR BayesC route are validated
  and used directly.
- Learned logistic probabilities are
  `logit^-1(logit(pi)+center(A eta))`, with numerical clipping. Centering the
  annotation predictor does not make the global-`pi` full conditional Beta.
  `samplePi_ST_annot()` is therefore invalid when `eta` is nonzero. The
  smallest validating test is a one-dimensional grid or slice sampler for
  `pi` compared with long-run sampler moments at fixed `eta` and states.
- SBayesRC/BayesRC continuation-probit coefficients are updated with the stated
  Gaussian conditionals and scaled-inverse-chi-squared hierarchy.

### Variance architecture `Q`

The maintained learned log-variance path uses centered annotations, excludes
an intercept from the relative predictor, exponentiates with guards, and
normalizes to geometric mean one. Fixed/group/external multipliers are checked
for positivity. The global variance update uses `b_j^2/q_j` for BayesC and
`b_j^2/(gamma_cj q_j)` for BayesR. No unweighted active-effect sum was found on
the maintained LV paths. Caps and exponent guards change the allowable
numerical parameter region and must remain documented as implementation
safeguards rather than mathematical necessities.

### Sharing `H` and covariance `V_b`

Ordinary MT BayesR has one joint simplex. MT BayesRC alone factorizes
component probabilities `P` and conditional non-null sharing `H`. Annotation-
dependent `H`, annotation-dependent `V_b`, and MT learned `Q` are not current
public implementations.

## 9. Likelihood-operator audit

### BED

The public BED routes decode packed PLINK genotypes, apply the route's declared
allele coding, centering, and scaling, impute missing genotypes consistently,
and maintain individual residuals. Marker diagonals are genotype sums of
squares on the individual likelihood scale. MT BED uses common individuals;
its fitted-value covariance is therefore an exact common-sample covariance.

### Sparse CSR

The CSR object represents a symmetric sparse cross-product-scale operator with
explicit diagonal and marker ownership/alignment checks. Omitted edges are
treated as zero by the fitted approximate operator. The runtime does not
establish positive semidefiniteness and does not repair or warn about an
indefinite matrix. `CsrOperator` nonetheless names its quadratic and residual
calculations as genetic variance and SSE, and final output uses `vgs`, `vle`,
`vld`, and `ves`. The algebra is coherent for the fitted operator, including
`vld = vg-vle` and possible negative `vld`, but the quantities are not
guaranteed biological variances when the operator is indefinite.

### Block eigen

For a block with positive marker diagonal `D_b`, the implementation whitens to
the correlation scale, retains eigenpairs under the configured proportion
rule, and constructs a runtime factor equivalent to

```text
Q_b = Lambda_(b,k)^(1/2) U_(b,k)^T D_b^(1/2),
w_b = Lambda_(b,k)^(-1/2) U_(b,k)^T D_b^(-1/2) s_b.
```

The source orientation makes `Q_b` retained-rank by marker, so
`Q_b^T Q_b=C_tilde_b` and `Q_b^T w_b=s_tilde_b`. Zero/nonpositive input
diagonals are rejected before whitening. The canonical statistical denominator
is the common global `N`. The score step is an orthogonal projection in
whitened coordinates, not a Euclidean projection of the raw score.

The `gctb_block`, `fixed_block`, and `global_projected_legacy` residual-policy
dispatch agrees with the maintained contract. Eigen truncation, block-boundary
separation, and float storage are intentional approximations. The GCTB-like
eigenspace likelihood does not establish complete GCTB/SBayesRC equivalence:
the residual policy, scale calibration, priors, defaults, and output contracts
also have to agree.

MT block eigen reconstructs block operators for the MT core; it is not the ST
retained `Q/w` scalar runtime and should not be used as evidence of retained
MT low-rank support.

## 10. Variance, covariance, and output audit

The ST and MT `vgs` traces are calculated from each current draw, not from the
posterior mean effect. Thus the implementation does not confuse
`E[b^T C b/N|D]` with `E[b|D]^T C E[b|D]/N`.

For BED and retained block operators, the genomic quadratic has a genuine
positive-semidefinite interpretation. For CSR it is operator-relative unless
PSD has separately been established. `v_LE` uses the fitted diagonal and
`v_LD` is the algebraic difference; a negative `v_LD` is valid and is not
clipped.

For MT CSR/block routes, `computeG()` forms diagonal trait quadratics and an
off-diagonal symmetric bilinear expression normalized by `sqrt(N_t N_s)`.
When traits use different samples/operators, this is descriptive scaling; it
does not recover paired genotype cross-products or identify a residual/genetic
covariance. The final names `cov_g_mean/final` do not disclose this limitation.

The marker-effect covariance issue is distinct. `sampleB_latent()` draws the
covariance used in the marker kernel, but `sampleB()` later constructs new
diagonal inverse-chi-squared draws, empirical shared-marker correlations,
`n_shared/(n_shared+20)` shrinkage, a `[-0.95,0.95]` clamp, and an eigenvalue
floor. The latter value is accumulated into `vbs`/`cov_b_*` and returned as
final `V_b`. It is not a draw from the inverse-Wishart conditional documented
immediately above `sampleB_latent()`.

CSR SBayesRC explicitly assigns posterior means to `alpha_final` and
`sigmaSqAlpha_final` because final draws are not retained. BED BayesRC retains
actual chain-final values. Shape/finite tests do not detect this semantic
difference.

The canonical public formatter renames raw `pis` to `pi_trace`. A formatted fit
does not consistently contain `fit$pis`, by design of the current schema.
`pi_trace` is family-dependent prior state mass: ordinary BayesC active mass,
BayesR total non-null mixture mass, or BayesRC mean prior non-null mass. It is
not a marker PIP. Marker PIPs are the canonical `dm` quantity; posterior
component frequencies are separate.

## 11. Single-trait and multi-trait comparison

| Property | ST | Ordinary MT BayesR | MT BayesRC |
|---|---|---|---|
| State prior | null/active or component simplex | full joint pattern-by-component simplex | marker component `P` times conditional non-null sharing `H` |
| Exact null | yes | one unique global null | one unique global null |
| Effective effects | scalar effect | latent vector masked by activity pattern | latent vector masked by `P x H` state |
| Effect covariance | scalar `v_b` | shared `V_b` with component multiplier | shared `V_b` with component multiplier |
| MAF-S | fixed and selected sampled routes | fixed only where dispatched | fixed only where dispatched |
| Summary overlap | not applicable | no arbitrary unidentified overlap inference | no arbitrary unidentified overlap inference |

The probability reductions are correct: BayesR PIP marginalizes all non-null
states; trait PIPs marginalize states active for that trait; component and
pattern marginals are not substituted for one another. MT BayesC is not
silently recoded as `P x H`, and ordinary MT BayesR is not harmonized with the
BayesRC factorization.

## 12. RNG, parallel, and convergence audit

Logical chains receive deterministic resolved seeds and own separate
`std::mt19937` engines. OpenMP uses static task ownership and does not share an
R RNG in worker regions. Chain aggregation occurs after workers complete.
Permanent reproducibility tests cover serial/parallel equality for maintained
small fixtures. No shared mutable marker state or thread-unsafe RNG access was
found in the maintained public paths.

Burn-in, thinning, retained-marker summaries, covariance traces, probability
traces, selected-marker diagnostics, and annotation traces were checked against
the convergence contract. One avoidable RNG-consumption issue exists: the
heuristic `sampleBset()` draw is immediately overwritten by
`sampleB_latent()` before marker updates. This does not by itself change the
posterior target, but it changes draw ordering and complicates reproducibility.

The extended fresh-process reproducibility harness did not run to completion
because its callr subprocess invoked `testthat::test_path()` without an active
testthat context. This is a harness failure, not evidence of a sampler race.

## 13. Fine-mapping audit

`make_credible_sets()`, `make_credible_sets_from_ld()`, and
`make_multisignal_credible_sets_from_ld()` operate on marginal PIPs. Their
cumulative sets estimate expected active-marker mass, not a joint probability
that a set contains a causal variant. The multi-signal helper is an approximate
LD-conditioned partitioning procedure, not SuSiE-style effect-specific
credible sets. `finemap_stblr_csr()` performs local refitting with explicit
marker and LD alignment; `extract_stblr_finemap_loci()` preserves the local
fit/post-processing distinction.

One input-contract gap remains: the public credible-set helpers replace
nonfinite PIPs with zero but do not reject finite values outside `[0,1]`.
Permanent tests should independently enforce the probability domain.

## 14. Test ownership and validation assessment

Permanent tests cover BayesC selection math, BayesR component probabilities,
SBayesRC stick probabilities and alpha hierarchy, MAF-S, log-variance math,
retained-block identities and reductions, residual policies, MT state models,
output schemas, logical-chain reproducibility, credible sets, and local fine
mapping. The focused tests passed in source-load mode.

Important weaknesses are:

- no independent target test for the learned-logistic global-`pi` conditional;
- no test that the reported MT `V_b` is the same posterior draw used by the
  marker kernel;
- CSR SBayesRC tests check only dimensions/finiteness of fields named `final`;
- no permanent test exercises indefinite CSR output terminology/diagnostics;
- no probability-domain test rejects finite PIPs outside `[0,1]`;
- architecture audit consumes the acknowledged-stale capability matrix, so its
  pass does not independently prove public support;
- two installed-package frozen MD5 gates fail while the corresponding
  source-load functional tests pass. The hashes cover serialized nested output
  rather than independently derived scientific quantities and therefore do not
  identify the numerical cause.

Passing a fixture generated by the same implementation was not treated as
proof of a conditional formula. Analytical/reference tests were weighted more
strongly than self-generated trajectory fixtures.

## 15. Prioritized findings

| ID | Severity | Primary classification | Affected route | Evidence and consequence | Smallest independent validation / correction owner |
|---|---|---|---|---|---|
| F01 | High | `CONFIRMED_IMPLEMENTATION_ERROR` | all public MT BED/CSR/block BayesC/R/RC with `updateB` | `sampleB_latent()` supplies the covariance used by marker updates, but `sampleB()` overwrites and reports a heuristic covariance. Reported `V_b` is not a posterior conditional draw; existing `cov_b_*`, `vbs`, and final `V_b` may be affected. Marker posterior target does not appear affected. | Freeze marker/latent effects and compare reported draws with inverse-Wishart theory. Correct `mtblr.cpp` and both MT cores; results/schema migration note required. |
| F02 | High | `CONFIRMED_IMPLEMENTATION_ERROR` | ST CSR learned-logistic BayesC, `updatePi=TRUE` | marker probabilities use offset logits, but `samplePi_ST_annot()` uses conjugate Beta counts. This changes the posterior target and may affect existing fits. | Fixed-state numerical conditional for global `pi` versus quadrature/slice sampling. Correct learned-annotation core and tests. |
| F03 | Moderate | `OUTPUT_SCHEMA_MISMATCH` | CSR SBayesRC | `alpha_final` and `sigmaSqAlpha_final` are assigned posterior means in `st_sbayesrc_omp_csr.cpp`. | Retain actual chain finals or rename/remove the fields through an intentional schema migration. |
| F04 | Moderate | `OUTPUT_SCHEMA_MISMATCH` | ST/MT CSR | potentially indefinite operator quadratics are exposed as unqualified variance/SSE fields. | Construct a symmetric indefinite fixture and verify labels/diagnostics. Correct schema/formatter/docs after choosing a PSD policy. |
| F05 | High | `OUTPUT_SCHEMA_MISMATCH` | MT CSR/block with trait-specific operators/samples | `computeG()` returns a descriptive `sqrt(N_t N_s)` bilinear under `cov_g_*`; paired covariance is not identified. | Two incompatible trait operators with known target covariance. Qualify or split output semantics. |
| F06 | Moderate | `OUTPUT_SCHEMA_MISMATCH` | all formatted fits | requested `fit$pis` consistency is absent; canonical output is `pi_trace`, whose meaning varies by family and is prior mass, not PIP. | Contract test across families. Decide whether the stable schema is sufficient or an intentional API migration is wanted. |
| F07 | Moderate | `DEVELOPER_CONTRACT_MISMATCH` | ST retained block | `R/stblr-block-eigen.R` roxygen still broadly contrasts package global projected residual variance with GCTB although current defaults select `gctb_block` for retained SBayesR/SBayesRC. | Update roxygen/public docs after source corrections; no sampler change. |
| F08 | Moderate | `DEVELOPER_CONTRACT_MISMATCH` | support auditing | `.blr_model_capability_matrix()` omits promoted MT BayesRC and LV capabilities but is consumed by architecture audit. | Replace audit source with executable dispatch/schema assertions; remove or clearly quarantine stale matrix. |
| F09 | High | `MISSING_OR_INADEQUATE_TEST` | F01/F02 routes | existing tests do not independently validate the two incorrect transitions/outputs. | Analytical conditional tests described in F01/F02. |
| F10 | Moderate | `MISSING_OR_INADEQUATE_TEST` | SBayesRC final fields, CSR semantics, credible sets | tests check shape/finite values but not semantic equality or PIP domain. | Add draw-identity, indefinite-operator, and `[0,1]` tests. |
| F11 | High | `MISSING_OR_INADEQUATE_TEST` | installed-package validation/reproducibility | two frozen serialized MD5 gates differ after install; opt-in callr harness errors before testing. No posterior error is established, but release validation is not clean. | Compare named numerical leaves across load-all and installed package; repair callr harness context. |
| F12 | Low | `UNSUPPORTED_INTERNAL_CODE` | dormant BayesRC coupling/tempering/PX/particle/MCEM and legacy paths | source existence does not establish public support. | Keep excluded from dispatch/support claims or remove in a separate cleanup task. |
| F13 | High | `UNRESOLVED_SCIENTIFIC_DECISION` | MT `V_b` | choose and document the intended covariance prior/augmentation and whether `V_b` output is a posterior draw or estimator. | Approve target before F01 correction. |
| F14 | Moderate | `UNRESOLVED_SCIENTIFIC_DECISION` | sparse CSR | decide whether to require PSD, diagnose/repair it, or retain algebraic operator-relative outputs with explicit labels. | Simulation with controlled indefinite operators. |
| F15 | Moderate | `UNRESOLVED_SCIENTIFIC_DECISION` | SBayesRC/block | decide whether exact published/GCTB parity is a target or only a parameterized analogue. | End-to-end primary-reference fixture including grid, priors, residual policy, scale, and likelihood. |
| F16 | Moderate | `EQUIVALENT_PARAMETERIZATION` | SBayesRC | relative ladder equals the published absolute ladder only under `v_b=.01 v_g^SNP`; current `v_b` is learned. | Crosswalk is valid but not an exact-prior claim. |
| F17 | Moderate | `INTENTIONAL_NUMERICAL_APPROXIMATION` | sparse and block operators | sparsification, block boundaries, rank truncation, and float storage alter the fitted operator, not merely evaluation speed. | Existing reduction tests plus approximation sensitivity analyses. |
| F18 | Low | `IMPLEMENTATION_SPECIFIC_SAFEGUARD` | probability/variance/block policies | floors, exponent guards, eigen floors, and rescue rules restrict numerical states. Most are appropriate; the heuristic `V_b` clamps are part of F01, not a harmless safeguard. | Boundary tests and explicit diagnostics. |
| F19 | Informational | `NO_ACTION` | maintained core routes | BayesC/R conditionals, SBayesRC stick algebra, MAF-S, MT state probabilities, `Q/w`, draw-level `v_g`, and marginal-PIP semantics agree with theory. | Retain current analytical/reduction tests. |

## 16. Verified current support

`SUPPORTED` requires public dispatch, executable implementation, coherent
output handling, documentation, and permanent test ownership. A route with a
substantive output or covariance limitation is marked `PARTIAL` even if it
runs.

| Public model | Data level | ST/MT | Available operators | Important restriction | Verified status |
|---|---|---|---|---|---|
| `bayesc` | individual | ST | BED | no MAF-S or LV on BED | `SUPPORTED` |
| `sbayesc` | summary | ST | CSR; retained or dense block eigen | CSR may be indefinite; operator-relative variance semantics | `SUPPORTED` |
| `bayesr` | individual | ST | BED | no MAF-S on BED | `SUPPORTED` |
| `sbayesr` | summary | ST | CSR; retained or dense block eigen | retained residual policy is explicit; not complete GCTB parity | `SUPPORTED` |
| `bayesrc` | individual | ST | BED | published/default prior settings are not identical across implementations | `SUPPORTED` |
| `sbayesrc` | summary | ST | CSR; retained or dense block eigen | CSR annotation `final` fields are means; relative scale is not automatically the published absolute prior | `PARTIAL` |
| `bayesc` | individual | MT | BED | reported `V_b` is heuristic; common-individual likelihood | `PARTIAL` |
| `sbayesc` | summary | MT | CSR; reconstructed-dense block eigen | reported `V_b`; descriptive off-diagonal `cov_g`; no arbitrary overlap model | `PARTIAL` |
| `bayesr` | individual | MT | BED | full joint simplex is correct; reported `V_b` is heuristic | `PARTIAL` |
| `sbayesr` | summary | MT | CSR; reconstructed-dense block eigen | reported `V_b`; descriptive off-diagonal `cov_g` | `PARTIAL` |
| `bayesrc` | individual | MT | BED | `P x H` state model is implemented; reported `V_b` is heuristic | `PARTIAL` |
| `sbayesrc` | summary | MT | CSR; reconstructed-dense block eigen | `P x H`; reported `V_b`; descriptive off-diagonal `cov_g` | `PARTIAL` |
| BayesC-V / SBayesC-V | summary | ST | CSR; retained block eigen | learned log-variance only; no BED/MT | `SUPPORTED` within restriction |
| BayesR-V / SBayesR-V | summary | ST | CSR; retained block eigen | learned log-variance only; no BED/MT | `SUPPORTED` within restriction |

The learned-logistic probability-provider route is `PARTIAL`: it is public and
executable, but its default global-`pi` update changes the intended posterior
target. Internal experimental samplers and unsupported operator combinations
are not included in the supported table.

## 17. Models and extensions in development

| Model or extension | Goal | Present state and evidence | Next validation gate | Status |
|---|---|---|---|---|
| Combined probability and variance annotations | compose `P` and `Q` | limited old ST CSR BayesC composition exists; architecture plan is broader; learned-logistic `pi` is incorrect | approve probability baseline and repair F02 with analytical tests | `IN_DEVELOPMENT` |
| External marker-specific variance priors | fixed positive external `q_j` | provider exists on restricted ST CSR BayesC route | cross-family dispatch/schema tests and scale-identifiability checks | `PARTIAL` |
| MT variance-informed models | shared learned `q_j` with `V_b` | theory/plan only; public MT-LV rejected | first resolve MT `V_b`, then derive weighted covariance conditional | `PLANNED` |
| Annotation-dependent trait sharing | marker-specific `H` | explicitly outside current `P x H` architecture | approve identifiability and probability model | `RESEARCH` |
| Annotation-dependent effect covariance | annotation-dependent `V_b` | conceptual research extension only | choose covariance parameterization and positive-definiteness strategy | `RESEARCH` |
| GCTB SBayesRC parity | exact external-reference agreement | eigenspace and probability crosswalks qualified; priors/scales/residual details differ | independent end-to-end primary-reference fixture | `IN_DEVELOPMENT` |
| Formal multi-signal fine mapping | joint/effect-specific credible sets | current helpers are marginal-PIP and approximate LD-conditioned procedures | specify joint configuration/effect model and calibrated simulations | `PLANNED` |

## 18. Proposed correction sequence

These are separate future tasks; none was performed during this audit.

1. **Approve theory targets and parameterization decisions.** Resolve the MT
   `V_b` prior/augmentation, CSR PSD/reporting policy, target-population MT
   covariance semantics, and whether exact GCTB/SBayesRC parity is required.
   Likely affected: Methods only if the approved target changes, then current
   contracts. Results could change; migration notes are required.
2. **Correct the confirmed posterior-target error.** Replace the learned-
   logistic global-`pi` Beta update with a valid conditional update or fix
   global `pi`. Affected: learned-annotation C++ core, R validation/defaults,
   analytical tests, contracts. Existing learned-logistic results can change.
3. **Correct MT covariance output/state handling.** Remove overwritten
   heuristic draws or explicitly separate posterior `V_b` from a descriptive
   estimator. Affected: `mtblr.cpp`, both MT cores, raw schema/formatter/tests.
   Marker trajectories should be tested for preservation; covariance outputs
   and reproducibility hashes will change. Backward-compatible field names may
   be retained only if meaning is corrected and migration documented.
4. **Correct output definitions and schemas.** Retain true SBayesRC final draws
   or rename fields; qualify CSR algebraic summaries and unequal-sample MT
   bilinear forms; decide `pi_trace`/`pis` policy. Affected source aggregation,
   formatters, schema, tests, Notes, and migration documentation.
5. **Strengthen operator consistency.** Add PSD diagnostics/contract tests,
   full-rank reductions, and target-population covariance fixtures. Results may
   change only if a repair policy is adopted.
6. **Complete annotation/SBayesRC validation.** Add absolute-vs-relative grid,
   hyperprior, residual-policy, and global-scale reference tests. Exact parity,
   if chosen, can change posterior results and defaults.
7. **Close test and validation gaps.** Add independent conditional tests,
   repair the callr harness, replace opaque whole-object MD5 gates with named
   scientific invariants plus narrowly scoped reproducibility hashes, and add
   PIP-domain tests.
8. **Update developer/public contracts and then create the roadmap.** Update
   capability authority, block-eigen roxygen, schemas, Notes, and examples only
   after corrected behavior is executable. The public roadmap should use the
   verified support tables in this audit, not internal source existence.

## 19. Validation commands and outcomes

All commands below were run in the private checkpoint copy unless explicitly
marked as original-repository Git inspection.

| Validation | Command or scope | Outcome |
|---|---|---|
| compile and load | `Rscript -e "Rcpp::compileAttributes('.'); pkgload::load_all('.', compile=TRUE)"` | Passed. Generated/compiled files existed only in the private copy. Compiler emitted unused-symbol warnings only. |
| maintained fast tests | `Rscript -e "pkgload::load_all('.', compile=FALSE, quiet=TRUE); devtools::test('.', filter='blr-(unified-public-contract\|model-semantics\|operator-reductions\|unified-convergence\|selected-marker-diagnostics\|extended-(covariance\|probability\|parameter)-diagnostics)\|mtblr-(bayesr-model\|bayesrc-model\|raw-schema)\|stblr-(selection-s\|annotation-backends\|raw-schema)\|backend-consistency')"` | Passed; one warning from `diag(V)` in a synthetic covariance diagnostic. |
| generated-interface audit | `Rscript tools/audit/blr_generated_interfaces_audit.R` | Passed: 71 wrappers, 71 registrations, 31 exports. |
| architecture audit | `Rscript tools/audit/blr_architecture_audit.R .` | Passed 19/19, with the authority limitation in F08. |
| documentation audit | `Rscript tools/audit/blr_documentation_audit.R .` | Failed only for the known out-of-scope historical link `docs/dev/sbayesrc_s_em_phase5c.md -> docs/dev/delta, alpha`; no Methods failure. |
| focused scientific tests | `Rscript -e "pkgload::load_all('.', compile=FALSE, quiet=TRUE); for (f in c('test-alpha-hierarchy-conditionals.R','test-individual-bayesrc.R','test-sbayesrc-alpha-reference.R','test-sbayesrc-s-probit-reference.R','test-sbayesrc-s-proper-intercept-reference.R','test-stblr-maf-effect-s.R','test-stblr-low-rank-operator.R','test-block-eigen-gctb-residual-policy.R','test-logvar-math.R','test-logvar-bayesc.R','test-logvar-bayesr.R','test-mtblr-bayesr-model.R','test-mtblr-bayesrc-model.R','test-credible-sets.R','test-finemap-stblr-csr.R','test-extract-stblr-finemap-loci.R')) testthat::test_file(file.path('tests/testthat', f))"` | Passed. |
| full suite | with `SBLR_RUN_EXTENDED_REPRODUCIBILITY=true`: `Rscript -e "pkgload::load_all('.', compile=FALSE, quiet=TRUE); devtools::test('.')"` | All ordinary tests passed; 7 expected short-chain/synthetic warnings. One opt-in fresh-process callr test errored before execution because `testthat::test_path()` lacked a test context. |
| build | `Rscript tools/check/check_package.R .` build phase | `R CMD build` passed. |
| installed package check | `Rscript tools/check/check_package.R .` check phase | 5,333 expectations passed, 7 skipped, 1 warning, 2 failed; check status 1 ERROR plus 1 installed-size NOTE. The two failures were frozen serialized MD5 gates for retained block BayesC and BayesR (`0e54...` expected vs `2e62...`; `aa329...` expected vs `699c...`). Functional source-load tests for the same paths passed. |

Not run: large-data studies, research-only Study06 experiments, obsolete phase-
numbered benchmarks, and external GCTB large reference runs. These are not
required by the maintained fast validation contract and would not independently
resolve the confirmed findings.

## 20. Unresolved scientific decisions and audit integrity

The following are choices, not algebraic errors:

- the intended prior/update for MT `V_b`, including whether latent augmentation
  or a different identified covariance model is authoritative;
- whether sparse operators must be PSD, should be repaired, or should expose
  explicitly algebraic operator-relative summaries;
- whether cross-trait summary output should require a declared compatible
  target-population covariance or retain a separately named descriptive
  unequal-sample bilinear;
- whether the target is exact published SBayesRC/GCTB parity or a documented
  related parameterization;
- whether `fit$pis` should be introduced in a future intentional schema
  migration or `pi_trace` remains the sole canonical trace field.

This audit made no implementation correction. It did not modify R, C++, tests,
Methods, Notes, existing developer contracts, generated interfaces, package
metadata, or workflows. The only repository change produced by the audit is
this report, left unstaged for review. No file was staged, committed, pushed,
tagged, published, reset, restored, cleaned, or discarded; no branch was
created or switched; and no pull request was opened.
