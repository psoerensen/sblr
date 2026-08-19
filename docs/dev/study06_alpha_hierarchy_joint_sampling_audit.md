# Study 06 alpha-hierarchy joint-sampling audit

## Status and scope

This package-side development audit started from clean `master` at
`f2e3647920ed7e8b1ea9d47a6571b3753285682a` (`Evaluate exact pairwise
BayesRC allocation updates`). The read-only `sblrbench` evidence checkout was
clean at `de31f62e182d8540488d4135df4c58f052a515d9`. No package default,
public sampler argument, likelihood, prior, mixture, benchmark specification,
or qualification threshold was changed. Formal Study 06 qualification was not
rerun and the final benchmark remains unauthorized.

The immutable v2 truth passed both required identities:

* specification `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`;
* truth `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.

It has 2,000 people (1,400 training/600 validation), 1,500 markers in 15
blocks, component counts 1,329/84/50/37, realized heritability 0.50, and
`gamma = c(0, .01, .1, 1)`. The four annotation columns are an exact all-ones
intercept, enriched binary, continuous, and null annotations.

The diagnostic environment was R 4.4.1 UCRT, x86_64-w64-mingw32, Rtools GCC
13.2, and OpenMP (`-fopenmp`). Truth regeneration used the pinned isolated
`sblr` 0.2.0 install carrying `RemoteSha f2e364...`; source diagnostics then
loaded this working checkout. Key installed dependencies were Rcpp 1.1.1,
RcppArmadillo 15.2.3.1, posterior 1.6.1, and testthat 3.3.2. Windows R used its
configured reference BLAS (no external BLAS name was reported).

## Evidence entering the audit

The paired power-isolation result in `sblrbench` remains unchanged: ordinary
BayesR/SBayesR and true-alpha-fixed BayesRC/SBayesRC converged, whereas all
learned-alpha fits failed. Informative learned fits nevertheless improved
causal ranking. The committed result reported learned-informative AUPRC of
0.5948 (BED) and 0.5407 (block eigen), versus baseline 0.3975 and 0.3270.
These are descriptive under non-converged learned chains, not qualification
claims.

## Current implementation and mathematical contract

All scalar routes use `st_bayesrc_update_annotation_prior()` in
`src/st_bayesrc_annotation_prior.h`. Packed BED calls it from
`run_bed_bayesrc_chain()` after a complete marker sweep and variance updates;
CSR and retained block eigen call the same owner from
`run_csr_sbayesrc_core()`. The operator-specific code supplies current marker
allocations but does not redefine the annotation conditional.

For component `c_i in {0,...,K-1}`, the indicator for stick `j` (zero based)
is `y_ij = I(c_i > j)`. Stick 1 uses all markers. Stick `j > 1` is eligible
only when the preceding indicator is one, equivalently `c_i > j-1`. Thus the
observed truth has eligible/continuation/stopping counts 1500/171/1329,
171/87/84, and 87/37/50. Eligible rows, outcomes, and annotation rows are
constructed together and their dimensions agree.

Given latent Albert--Chib variables `z`, current residualized latent response,
and coefficient `k`, the scalar update is

```text
precision = x_k'x_k + prior_precision
variance  = 1 / precision
mean      = variance * (x_k' residual_without_k +
                        prior_precision * prior_mean)
```

The intercept uses its committed proper stick-specific normal prior. The first
annotation column was verified to be exactly one and is treated only as the
intercept. It is excluded from annotation-variance learning. Each of the three
remaining coefficients has `N(0, sigmaSqAlpha_j)` prior.

With `q = 3`, current parameters `a` and `b` produce

```text
sigmaSqAlpha_j | alpha[-1,j]
  = (sum(alpha[-1,j]^2) + b) / chi_square(q + a).
```

Equivalently the prior is inverse-gamma with shape `a/2` and scale `b/2`, or
scaled inverse-chi-square with `nu0 = a` and `scale0 = b/a`.

| Setting | a | b | nu0 | scale0 | prior mean | prior variance | survival-tail power |
|---|---:|---:|---:|---:|---:|---:|---:|
| current/teaching | 2 | 2 | 2 | 1 | infinite | infinite | 1 |
| production-equivalent | 4 | 4 | 4 | 1 | 2 | infinite | 2 |

The non-intercept count, not the intercept, enters the degrees of freedom.
Empty sticks draw from the proper prior; an unavailable non-intercept column
on a non-empty eligible subset also receives its prior conditional. Complete
separation remains finite under the proper intercept prior.

Conditional continuation probabilities are `Phi(A alpha_j)`. Sequential
component probabilities are `(1-p1)`, `p1(1-p2)`, ... and the final product.
Each component probability is floored at `1e-12`, followed by row
normalization. Deterministic ordinary states were far from the floor; it is a
numerical guard, not the source of this audit's mixing result.

## Component-trace defect and correction

The one-marker Study 06 smoke reproduced the native recursive termination.
`stblr_bed()` resolved character/numeric trace markers against every rsid in
the full BED Glist, but native packed-BED chains index the fitted `cls` subset.
For Study 06, even a valid fitted marker could therefore arrive as an index
greater than 1,500 and be dereferenced inside an OpenMP worker.

The R boundary now resolves trace markers against the canonical fitted marker
data produced by `.make_bed_marker_data()`. Both scalar BED BayesR and BayesRC
native entry points additionally validate every selected index before parallel
execution. This changes storage validation only, not any sampler transition.
The exact 20 x 4 x 1 Study 06 smoke passed, as did each full
9000 x 4 x 1500 dynamic trace. Tiny regression tests cover fitted-subset
identity, one marker, several/all tiny-model markers, multiple chains, invalid
indices, disabled tracing, RNG-neutral output, and CSR/block-eigen trace
contracts already exercised by the operator tests.

## Independent conditional validation

`research/sbayesrc/tools/study06_alpha_hierarchy_reference.R` independently constructs stick
eligibility, residualized Gaussian conditionals, the inverse-chi-square
conditional, and the stick-to-component transform. Native results agreed with:

* all eligible, continuation, and stopping counts exactly;
* alpha means and variances to `1e-14`;
* the variance conditional mean and 0.1/0.5/0.9 quantiles within preregistered
  Monte Carlo tolerances (50,000 draws);
* component probabilities and row normalization to `1e-15`;
* existing empty, one-information-direction, separated, and intercept-only
  posterior tests.

No conditional-formula, intercept-column, trace-retention, probability, or
BED/block route-adapter correctness defect was found.

## Frozen-allocation experiment

Truth allocations were fixed; phenotype, effects, and variance state were not
updated. Four seeds 860101/860202/860303/860404 ran 12,000 iterations with
3,000 burn-in and 9,000 retained draws. Every nonconstant monitored quantity
passed R-hat <= 1.01, bulk/tail ESS >= 400, and relative MCSE <= 0.05.

| Condition | max R-hat | min bulk ESS | min tail ESS | max rel. MCSE | Result |
|---|---:|---:|---:|---:|---|
| F1 current a=2,b=2 | 1.0015 | 2759.8 | 5017.6 | .0191 | pass |
| F2 production a=4,b=4 | 1.0022 | 2962.7 | 6235.4 | .0184 | pass |
| F3 fixed sigma=1 | 1.0015 | 2952.0 | 6201.7 | .0184 | pass |
| F4 alpha fixed to truth | maximum marker-prior error 0 | | | | pass |

This rejects H1 as a primary explanation and rejects H4: neither sampled
variance nor the probit conditional fails when stick outcomes are fixed.

## Dynamic-allocation ablations

All fits used fit seed 701020, chain seeds
701121/701222/701323/701424, 9,000 iterations, the existing 3,000 + 6,000
selected window, the unchanged phenotype/operator/priors other than the stated
variance ablation, and complete component traces.

| Condition | Route | Failed/monitored* | max R-hat | min bulk | min tail | max rel. MCSE | Runtime s |
|---|---|---:|---:|---:|---:|---:|---:|
| D1 fixed sigma=1 | BED | 29/32 | 1.0995 | 34.5 | 109.4 | .1723 | 472.7 |
| D2 a=4,b=4 | BED | 29/35 | 1.0854 | 38.0 | 75.1 | .1630 | 421.0 |
| D1 fixed sigma=1 | retained 0.995 | 29/32 | 1.0703 | 47.0 | 82.6 | .1516 | 288.2 |
| D2 a=4,b=4 | retained 0.995 | 32/35 | 1.1175 | 28.4 | 44.5 | .1869 | 240.6 |

`*` Constant fixed-variance traces are excluded rather than counted as failed
nonconstant estimands. Full quantity-level tables are retained locally.

Per-stick alpha/variance results were:

| Condition | stick 1 failed | stick 2 failed | stick 3 failed | sigmaSqAlpha result |
|---|---:|---:|---:|---|
| D1 BED | 4/4 | 4/4 | 4/4 | fixed, not an estimand |
| D2 BED | 4/5 | 4/5 | 4/5 | all three pass individually |
| D1 block | 4/4 | 4/4 | 4/4 | fixed, not an estimand |
| D2 block | 5/5 | 5/5 | 4/5 | sticks 1/2 fail; stick 3 passes |

The complete alpha/variance convergence audit is below. Fixed sigma rows are
reported as `fixed`; they are not failed random estimands.

| Condition | Stick | Parameter | R-hat | bulk/tail ESS | rel. MCSE | Result |
|---|---:|---|---:|---:|---:|---|
| D1 BED | 1 | intercept | 1.0995 | 35.0/109.4 | .169 | fail |
| D1 BED | 1 | enriched | 1.0553 | 92.6/290.7 | .104 | fail |
| D1 BED | 1 | continuous | 1.0523 | 92.3/251.5 | .107 | fail |
| D1 BED | 1 | null | 1.0752 | 45.8/117.9 | .159 | fail |
| D1 BED | 2 | intercept | 1.0208 | 180.5/220.8 | .073 | fail |
| D1 BED | 2 | enriched | 1.0201 | 392.5/870.0 | .051 | fail |
| D1 BED | 2 | continuous | 1.0089 | 374.8/572.0 | .052 | fail |
| D1 BED | 2 | null | 1.0841 | 34.5/624.6 | .172 | fail |
| D1 BED | 3 | intercept | 1.0070 | 338.8/416.1 | .055 | fail |
| D1 BED | 3 | enriched | 1.0141 | 265.5/617.4 | .061 | fail |
| D1 BED | 3 | continuous | 1.0120 | 346.5/688.2 | .054 | fail |
| D1 BED | 3 | null | 1.0130 | 404.0/1008.1 | .050 | fail |
| D1 BED | all | sigmaSqAlpha | fixed | fixed | fixed | not sampled |
| D2 BED | 1 | intercept | 1.0854 | 38.0/128.8 | .163 | fail |
| D2 BED | 1 | enriched | 1.0363 | 91.1/381.4 | .103 | fail |
| D2 BED | 1 | continuous | 1.0428 | 86.7/355.7 | .109 | fail |
| D2 BED | 1 | null | 1.0218 | 141.3/212.0 | .084 | fail |
| D2 BED | 1 | sigmaSqAlpha | 1.0053 | 1014.5/7014.8 | .021 | pass |
| D2 BED | 2 | intercept | 1.0249 | 181.8/435.3 | .074 | fail |
| D2 BED | 2 | enriched | 1.0119 | 199.7/197.9 | .076 | fail |
| D2 BED | 2 | continuous | 1.0182 | 218.6/193.8 | .074 | fail |
| D2 BED | 2 | null | 1.0107 | 217.2/351.4 | .068 | fail |
| D2 BED | 2 | sigmaSqAlpha | 1.0034 | 560.4/786.0 | .041 | pass |
| D2 BED | 3 | intercept | 1.0164 | 285.3/270.2 | .060 | fail |
| D2 BED | 3 | enriched | 1.0261 | 167.2/109.0 | .080 | fail |
| D2 BED | 3 | continuous | 1.0180 | 210.4/187.9 | .079 | fail |
| D2 BED | 3 | null | 1.0193 | 287.2/309.2 | .062 | fail |
| D2 BED | 3 | sigmaSqAlpha | 1.0079 | 401.1/420.1 | .047 | pass |
| D1 block | 1 | intercept | 1.0287 | 76.8/131.9 | .114 | fail |
| D1 block | 1 | enriched | 1.0249 | 97.0/114.9 | .104 | fail |
| D1 block | 1 | continuous | 1.0703 | 47.0/96.7 | .152 | fail |
| D1 block | 1 | null | 1.0536 | 85.3/140.8 | .109 | fail |
| D1 block | 2 | intercept | 1.0444 | 124.2/303.3 | .090 | fail |
| D1 block | 2 | enriched | 1.0205 | 183.1/365.4 | .075 | fail |
| D1 block | 2 | continuous | 1.0293 | 107.9/295.1 | .101 | fail |
| D1 block | 2 | null | 1.0548 | 121.0/246.2 | .091 | fail |
| D1 block | 3 | intercept | 1.0165 | 177.2/194.6 | .075 | fail |
| D1 block | 3 | enriched | 1.0173 | 292.7/642.0 | .058 | fail |
| D1 block | 3 | continuous | 1.0115 | 309.6/498.6 | .057 | fail |
| D1 block | 3 | null | 1.0190 | 137.9/310.0 | .087 | fail |
| D1 block | all | sigmaSqAlpha | fixed | fixed | fixed | not sampled |
| D2 block | 1 | intercept | 1.0616 | 49.1/66.7 | .147 | fail |
| D2 block | 1 | enriched | 1.0737 | 42.6/45.1 | .172 | fail |
| D2 block | 1 | continuous | 1.0657 | 46.5/66.0 | .153 | fail |
| D2 block | 1 | null | 1.1175 | 28.4/85.9 | .187 | fail |
| D2 block | 1 | sigmaSqAlpha | 1.0308 | 105.4/168.9 | .092 | fail |
| D2 block | 2 | intercept | 1.0465 | 84.8/109.2 | .107 | fail |
| D2 block | 2 | enriched | 1.0442 | 81.0/117.0 | .113 | fail |
| D2 block | 2 | continuous | 1.0489 | 79.7/98.7 | .117 | fail |
| D2 block | 2 | null | 1.0342 | 135.3/228.9 | .082 | fail |
| D2 block | 2 | sigmaSqAlpha | 1.0128 | 279.3/371.3 | .065 | fail |
| D2 block | 3 | intercept | 1.0246 | 158.8/314.8 | .079 | fail |
| D2 block | 3 | enriched | 1.0295 | 155.2/359.6 | .081 | fail |
| D2 block | 3 | continuous | 1.0180 | 123.4/201.4 | .097 | fail |
| D2 block | 3 | null | 1.0761 | 50.7/321.0 | .138 | fail |
| D2 block | 3 | sigmaSqAlpha | 1.0065 | 880.5/2581.5 | .030 | pass |

Thus later sticks are weak, but this is not only a late-stick problem: first
stick and null-versus-non-null occupancy fail in every dynamic condition.
Under D2, BED sigma ranges by stick were 0.211--43.68,
0.150--64.36, and 0.155--195.47; block ranges were 0.156--96.03,
0.143--41.08, and 0.121--27.95. The production prior makes variance tails
lighter but does not reconcile allocation regimes.

Mean dynamic occupancies (components 1/2/3; active) were 17.1/12.1/29.4;
58.7 (D1 BED), 15.0/12.0/29.3; 56.3 (D2 BED), 32.2/52.3/35.0;
119.5 (D1 block), and 49.9/50.9/35.5; 136.3 (D2 block). Active-count
ranges were 19--252, 20--190, 35--392, and 31--415, respectively.
Occupancy and expected-active diagnostics failed in every condition.
Selected-window mean expected-active counts by chain were
53.88/58.90/58.08/63.68 (D1 BED), 53.52/56.81/54.18/60.45
(D2 BED), 116.10/133.98/115.89/110.07 (D1 block), and
124.04/144.48/161.17/114.23 (D2 block). Corresponding realized occupancy
means closely tracked them. These chain-location differences confirm that the
failure is not caused by a single rare outlying draw.

Variance components remained finite. Mean BED genetic/residual variance and
heritability were 0.827/1.155/.417 in both D1 and D2. Block means were
1.005/0.976/.507 (D1) and 1.004/0.977/.506 (D2). The persistent route offset
is therefore still separate from hierarchy mixing.

## Posterior geometry

The stronger sampled prior produced moderate alpha-scale coupling, not a
single dominant funnel: maximum absolute chain-specific correlation between
`sigmaSqAlpha` and the non-intercept alpha norm was .595 (BED) and .590
(block). Sigma versus expected active count was at most .214 (BED) and .478
(block). Recorded alpha/occupancy correlations reached .623 and .803.
Fixed-variance fits still had severe occupancy and alpha failures, proving that
the observed hierarchy--occupancy feedback survives removal of the variance
dimension. A centred variance funnel is a secondary contributor in some
chains, not the primary bottleneck.

Core variance/active-count dependence was smaller for BED (maximum absolute
chain correlation .239 for effect variance and .169 across the other variance
summaries) and moderate for block eigen (up to .248 for effect variance, .330
for genetic variance, .462 for residual variance, and .453 for heritability).
This does not erase the separate BED/block calibration offset. The retained
operator rebuilt projected residuals 364 times across logical chains in each
condition; maximum absolute drift was `3.375078e-14`, with no integrity
failure.

## Causal ranking and prediction

| Condition | AUPRC | AUROC | recall 10/25/50/100 | precision 10/25/50/100 | FDR5 selected/true | FDR10 selected/true | prediction cor. |
|---|---:|---:|---|---|---|---|---:|
| D1 BED | .5913 | .8469 | .058/.140/.240/.404 | 1/.96/.82/.69 | 19/19 | 23/22 | .6546 |
| D2 BED | .5956 | .8526 | .058/.140/.246/.404 | 1/.96/.84/.69 | 19/19 | 23/22 | .6548 |
| D1 block | .5448 | .8330 | .058/.135/.216/.386 | 1/.92/.74/.66 | 24/22 | 32/27 | .6319 |
| D2 block | .5578 | .8385 | .058/.135/.228/.398 | 1/.92/.78/.68 | 24/23 | 33/27 | .6320 |

The ranking benefit remains descriptively strong, but none of these learned
hierarchies passes the retained convergence contract. Ranking power and joint
hierarchical reliability remain distinct conclusions.

## Decision

Decision C is primary: allocation feedback is the demonstrated bottleneck.
The exact same hierarchy converges rapidly with truth allocations frozen, but
fails with allocations dynamic under both fixed variance and the
production-equivalent prior, on both BED and retained block eigen. Decision F
is secondary because sampled variance shows moderate coupling and broad tails,
but neither fixing it nor regularizing it resolves the feedback. Decision E is
not adequate because first-stick inference also fails. No evidence supports D.

The next package task should diagnose posterior-preserving allocation/hierarchy
transition strategies (for example update frequency, blocked or partially
collapsed allocation, or interweaving) against this exact truth. It must be a
separate method-design review. This audit does not authorize pairwise or
collapsed allocation changes, default prior changes, requalification, or the
final benchmark.

## Validation

All changed R files and the machine decision JSON parsed. A clean
`R CMD INSTALL --preclean` into an isolated ignored library succeeded with GCC
13.2 and OpenMP. Focused hierarchy/component-trace/operator tests passed 693
expectations. The complete testthat suite passed 3,991 expectations with zero
failures/errors, one pre-existing covariance warning, and one opt-in fresh-
process reproducibility skip. Packaged `R CMD check --no-manual` completed with
zero errors and warnings and one installed-size NOTE. The direct directory
check was not authoritative because R 4.4.1 did not derive Author/Maintainer
from `Authors@R`; checking the standard source tarball resolved that metadata
step. Temporary tarballs, check directories, installed libraries, DLLs, and
object files were removed after validation.

Compact CSV evidence remains ignored under
`results/local/study06_alpha_hierarchy_audit/`. Raw fit objects and histories
were not tracked.
