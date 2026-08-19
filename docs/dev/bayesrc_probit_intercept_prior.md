# Proper probit-stick intercept prior

## Decision record

This note preregisters the intercept-scale decision before any new Study 06
pilot is examined. Study 06 is a diagnostic target, not a tuning set.

The previous BayesRC/SBayesRC implementation used a flat prior for the first
coefficient of every probit stick. A later stick can contain only continuation
outcomes, only stopping outcomes, or no eligible markers. The likelihood then
has an unbounded or absent intercept direction. The preserved Study 06 BED
histories exhibit exactly this configuration, so longer runs or a different
Gibbs order cannot identify the same flat-prior target.

## Model and centring

For stick `j`, `s_ij = 1` means that marker `i` continues beyond component
`j`. The package transform is

```text
Pr(component = j) = product_{l < j}(p_l) (1 - p_j)
Pr(component = K) = product_{l < K}(p_l)
p_j = Phi(A_i alpha_j).
```

For resolved global component probabilities `pi_0, ..., pi_K`, the baseline
continuation probability is

```text
p_j = sum_{k > j} pi_k / sum_{k >= j} pi_k,
```

and the default intercept mean is `mu_0j = qnorm(p_j)`. This mean is part of
the model. `alpha_init` remains only an MCMC starting state and never changes
the prior.

The new default is independently proper:

```text
alpha_0j ~ Normal(mu_0j, tau_0j^2).
```

Non-intercept coefficients retain the existing zero-centred hierarchical
normal prior and the existing `sigmaSqAlpha` update.

## Preregistered scale selection

The candidate grid is `0.5, 1.0, 1.5, 2.5`. The deterministic calibration
uses four component architectures (balanced, sparse equal-active, moderately
sparse, and skewed-active), none obtained from Study 06. A stick probability
is called practically degenerate below `1e-8` or above `1 - 1e-8`.

The selection rule is: choose the largest candidate SD (the least informative
candidate) for which the normal prior puts at most 1% probability on the
practically-degenerate region for every representative stick. The calibration
also reports 95% prior-predictive probability widths. Mixed-outcome synthetic
checks require that 50 or more balanced observations overwhelm the prior;
complete-separation checks require finite posterior moments without coefficient
caps. Run `Rscript research/sbayesrc/tools/calibrate_bayesrc_intercept_prior.R` to reproduce the
table. The selected value is recorded below only after running that script.

The maximum degenerate mass over the four architectures was approximately
`2.3e-7`, `5.84e-3`, `4.64e-2`, and `1.57e-1` for SDs 0.5, 1.0, 1.5, and 2.5,
respectively. The preregistered rule therefore selects default SD **1.0**.

## Full conditional

With Albert--Chib latent responses `z`, design `X`, and
`alpha ~ Normal(m0, P0^-1)`, the Gaussian conditional is

```text
P = X'X + P0
m = P^-1 (X'z + P0 m0)
alpha | z, ... ~ Normal(m, P^-1).
```

The coefficient-wise kernel therefore adds `1/tau_0j^2` to the intercept
likelihood diagonal and adds `mu_0j/tau_0j^2` to its right-hand side. It does
not include the intercept in the `sigmaSqAlpha` conditional.

## Empty and separated sticks

An empty stick contributes no likelihood. One Gibbs transition samples its
intercept from the proper intercept prior, samples each non-intercept
coefficient from its conditional hierarchical prior, and then updates
`sigmaSqAlpha` from its unchanged conditional. All-one and all-zero sticks use
the ordinary latent-variable update and remain finite because the intercept
prior is proper.

The historical flat prior is retained only through explicit
`intercept_flat = TRUE`. It emits a warning and fails for an empty or completely
separated stick, because those posterior targets are improper. Historical
`intercept_flat = FALSE` coupled the intercept to `sigmaSqAlpha`; that ambiguous
model is intentionally rejected in favour of an explicit prior.

## API and native contract

`annotation_intercept_prior` accepts a list with `distribution = "normal"`,
`mean = "initial_mixture"` or finite scalar/stick-specific means, and positive
finite scalar/stick-specific `sd`. One R resolver records means, SDs,
variances, precisions, component probabilities, stick probabilities, and the
legacy status. It passes one 3-by-stick native matrix containing type, mean,
and precision. `alpha_init` is passed separately.

The shared owner is `src/st_bayesrc_annotation_prior.h`; scalar BED, scalar
CSR, retained block eigen, multitrait BED, multitrait CSR, and multitrait block
eigen all call it. Operator approximation affects the marker likelihood. It
does not alter this annotation-prior identification contract.

## References and implementation comparison

Albert and Chib (1993, JASA 88:669--679,
doi:10.1080/01621459.1993.10476321) provide the latent-normal augmentation.
Zheng et al. (2024, Nature Genetics 56:767--777,
doi:10.1038/s41588-024-01704-y) define SBayesRC's annotation-informed mixture.
The manuscript repository (`zhilizheng/SBayesRC`, `main` commit
`2881f9a6f57374b80b4d19a09dbf9387939a2f46`) and GCTB (`jianzeng/GCTB`,
`master` commit `cc7fa7d765c83a89c6375946cf77fe50ba1a317e`) were inspected as
implementation references. Their empty-stick heuristics are not used as the
prior model here.

## Limitations

A proper prior identifies separation but does not guarantee rapid mixing of
all annotation, mixture, and variance parameters. Study 06 package-development
pilots are required after deterministic and synthetic tests; they are not
qualification evidence.

## Development validation (2026-08-04)

Deterministic tests cover mixture-to-stick round trips, scalar and
stick-specific priors, prior/initialization separation, prior-only empty-stick
updates, legacy failures, and reproducibility. Intercept-only Gibbs chains for
balanced, imbalanced, near-separated, all-continuation, and no-continuation
data agreed with numerical-grid posterior means and variances. Existing
operator-reduction, seed, thread, trace-retention, and truncated-normal tests
also passed in the full package suite.

Non-qualification Study 06 replicate-1 BED pilots used 300 iterations in each
of the original four chains. Both scenarios completed with finite annotation
and variance traces without caps or clipping. Informative-chain intercepts
ranged from -3.20 to 2.75; the largest absolute non-intercept coefficient was
8.04 and `sigmaSqAlpha` ranged from 0.114 to 103.7. Uninformative-chain
intercepts ranged from -3.21 to 3.23; the largest absolute non-intercept
coefficient was 7.54 and `sigmaSqAlpha` ranged from 0.089 to 46.99. These short
runs demonstrate finite identification, not satisfactory long-run mixing.

One-chain, 100-iteration retained-0.995 block-eigen pilots also completed for
both scenarios. Informative residual variance ranged from 0.769 to 1.147 and
maximum projected-residual drift was `2.22e-15`; uninformative residual
variance ranged from 0.713 to 1.138 and drift was `1.33e-15`. Genetic variance
remained positive. All checked sibling evidence files were byte-identical.

The remaining blocker is a longer multi-chain mixing assessment of the proper
posterior, especially the non-intercept coefficients and `sigmaSqAlpha` on
late, sparsely occupied sticks. The package change is therefore not by itself
a Study 06 qualification result.

The subsequent focused review is recorded in
`docs/dev/bayesrc_annotation_mixing_review.md`. It corrected one prior-only
conditional for annotation columns unavailable on the current eligible subset,
but did not change the prior or introduce a blocked/non-centred kernel.
