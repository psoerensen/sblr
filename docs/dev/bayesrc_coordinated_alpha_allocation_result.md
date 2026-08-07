# Coordinated BayesRC alpha--allocation transition result

## Decision

**AA-R5 -- no valid scalable coordinated transition was derived.**

The exact beta-inclusive construction in
`bayesrc_coordinated_alpha_allocation_transition.md` is posterior preserving
and passes a tiny-model oracle. It is not a viable production transition for
the demonstrated Study 06 bottleneck. With a fixed, scientifically meaningful
alpha proposal, its acceptance penalty grows with the number of markers left
outside the exactly refreshed subset. Shrinking the proposal by
`1 / sqrt(m)` preserves acceptance but makes the alpha movement vanish. Making
the refreshed subset global instead costs `K^|S|` component configurations.

No sampler transition, public argument, default, native route, or returned-fit
contract was changed. In particular, ordinary BED, CSR, and block-eigen
BayesRC/SBayesRC retain their exact historical RNG trajectories.

## Provenance and scope

- package source: `sblr` 0.2.0;
- source commit: `0c89234273389e14112ba0e08ef9d11d3e1819dc`;
- decision date: 2026-08-07;
- sibling `sblrbench` was inspected read-only at
  `fbe80603ff6fa09e0a611a56d09130cb4b2cbc8c`;
- no Study 06 qualification, large scientific fit, BED fit, CSR fit, or
  block-eigen fit was run in this task.

The validation used only an independent R development oracle. Compact results
are written to the ignored local directory
`results/local/bayesrc_coordinated_transition/`.

## Candidate decisions

### A: alpha plus allocation labels with beta fixed

Rejected on support, before tuning. The null component is a point mass at
`beta = 0`, whereas active components have continuous Gaussian effect priors.
A nonzero beta cannot be assigned to the null component, and an exactly zero
beta has probability zero under an active Gaussian component. A beta-fixed
move can therefore change labels only among active components. It cannot move
the realized active count and cannot address the observed bottleneck.

### B: collapsed alpha

Rejected as non-scalable. The Gaussian alpha prior times the probit stick
likelihood has no closed-form integral. Albert--Chib augmentation converts the
integral to a Gaussian orthant probability with dimension equal to the stick's
eligible marker count. Tiny quadrature is possible; an exact Study 06-scale
transition is not.

### C: alpha plus allocations plus beta

Selected for the exact feasibility oracle. For a state-independent subset
`S`, all `K^|S|` component configurations were enumerated. Effects were
integrated out when choosing a configuration and then sampled jointly from
their exact Gaussian conditional. With symmetric alpha proposal `q`, the log
Metropolis--Hastings ratio is

\[
\begin{aligned}
\log r={}&\log p(\alpha_j'\mid\sigma_{\alpha,j}^2)
-\log p(\alpha_j\mid\sigma_{\alpha,j}^2)\\
&+\sum_{i\notin S}
 [\log\Pr(c_i\mid\alpha')-\log\Pr(c_i\mid\alpha)]\\
&+\log Z_S(\alpha')-\log Z_S(\alpha),
\end{aligned}
\]

where `Z_S` is the sum of the exact collapsed likelihood/effect-prior/
allocation weights over the subset configurations. This cancellation was
verified against an independently evaluated full target and forward/reverse
proposal density.

## Tiny exact oracle

The registered fixture has eight independent marker likelihood terms, two
components (null and active), an intercept plus one binary annotation, a
proper `N(Phi^-1(0.15), 1)` intercept prior, a `N(0, 1)` annotation prior,
fixed residual variance 1, and active effect variance 0.7. The exact reference
enumerates all 256 allocations and integrates the two alpha dimensions on a
181 by 181 grid. The maximum omitted univariate prior-tail bound was
`4.70e-4`.

Four chains were run for each kernel, using seeds 12011, 12022, 12033, and
12044. Each chain used 30,000 iterations, 5,000 burn-in iterations, and
thinning by five. Monte Carlo tolerances were fixed as four times the
between-chain Monte Carlo standard error before interpreting the output.

| Kernel | max absolute alpha-mean error | max PIP error | max beta-mean error | active-count distribution L1 |
|---|---:|---:|---:|---:|
| ordinary Gibbs | 0.0131 | 0.00693 | 0.00391 | 0.0200 |
| coordinated only | 0.00819 | 0.00438 | 0.00411 | 0.0126 |
| ordinary + coordinated | 0.0168 | 0.00550 | 0.00200 | 0.0145 |

All alpha and PIP errors were below the preregistered Monte Carlo tolerances.
The exact reference alpha mean was `(-1.44075, -0.15632)` and all three
kernels recovered it. The coordinated-only and combined chains attempted
120,000 moves each, accepted 81.9%, and made 23,284 and 23,011 accepted
nonzero active-count changes, respectively. The largest possible and observed
jump was two because the refreshed subset contained two markers.

The independent detailed-balance comparison covered the full target,
collapsed configuration proposal, conditional beta density, and the proper
alpha prior. Forward and reverse log flows agreed at tolerance `1e-10`.
One-dimensional numerical integration independently reproduced the collapsed
component weights at tolerance `1e-10`. All-null, all-active, and rare-state
support was finite, and deterministic diagnostic helpers consumed no RNG.

## Scalability gate

The same two-marker exact refresh and alpha random walk were evaluated on
fixed synthetic marker arrays. These are feasibility calculations, not
smaller-model posterior samples.

| markers | alpha scale | median log MH | mean acceptance | median acceptance | mean alpha jump |
|---:|:---|---:|---:|---:|---:|
| 8 | fixed | -0.120 | 0.815 | 0.887 | 0.310 |
| 50 | fixed | -0.524 | 0.552 | 0.592 | 0.314 |
| 100 | fixed | -1.378 | 0.357 | 0.252 | 0.318 |
| 250 | fixed | -3.191 | 0.241 | 0.0411 | 0.315 |
| 500 | fixed | -6.726 | 0.137 | 0.00120 | 0.308 |
| 1,500 | fixed | -23.773 | 0.0760 | `4.80e-11` | 0.311 |
| 37,991 | fixed | -585.009 | 0.00574 | `8.60e-255` | 0.312 |
| 1,500 | `1/sqrt(m)` | -0.165 | 0.788 | 0.848 | 0.0235 |
| 37,991 | `1/sqrt(m)` | -0.117 | 0.828 | 0.889 | 0.00446 |

The mean acceptance at large marker counts is misleading because it is
dominated by a small upper tail; the median proposal is effectively never
accepted. The locally scaled proposal remains accepted but becomes an
ordinary infinitesimal random walk. It cannot coordinate the tens-of-markers
occupancy changes implicated by Study 06. Increasing the subset size enough
to absorb those changes is not feasible: with the Study 06 four-component
mixture, sizes 8, 10, and 20 require 65,536, 1,048,576, and about 1.1 trillion
configuration evaluations per proposal.

## Production and route consequences

The scalability gate failed before production integration. Consequently:

- no `alpha_allocation_transition` argument or experimental default was added;
- no native C++ or Rcpp binding changed;
- BED, CSR, and block-eigen routes require no route-specific equivalence test
  because none invokes the development reference kernel;
- the GCTB-compatible block residual policy and all existing likelihood,
  prior, annotation, allocation, chain, OpenMP, and output contracts are
  unchanged;
- ordinary-mode RNG identity follows structurally: the production call graph
  and every sampler instruction are byte-for-byte unchanged.

A 1,500-marker Study 06 comparison was deliberately not launched. The
predeclared requirement was to proceed only after deriving a transition that
was both exact and computationally capable of crossing occupancy/alpha
regimes. Candidate C is exact but fails that gate.

## Recommendation and limitations

Do not promote a small-subset alpha/allocation/beta MH move. It duplicates
local behavior already supplied more efficiently by ordinary Gibbs and the
existing exact pair update, while its global alpha proposal pays an extensive
outside-marker penalty.

The next method task, if pursued, should derive a genuinely global proposal
whose allocation and beta state can be marginalized or transported without
`K^m` enumeration. Candidates include a rigorously collapsed global move or a
deterministic transport with a computable Jacobian and reverse density. That
is new method development and requires its own exact tiny-model design gate;
it should not begin with `sblrbench` requalification.

Limitations are intentional: the oracle fixes `sigmaSqAlpha`, uses a diagonal
marker likelihood, and uses two mixture components. Those simplifications are
favorable to Candidate C. Adding correlated markers, more mixture components,
or sampled annotation variances increases rather than removes its enumeration
and coupling burden.

## Validation record

- focused reference tests: 20 passed, 0 failed, 0 warned, 0 skipped;
- complete `devtools::test()`: 4,755 passed, 0 failed, 1 known environment
  warning, and 1 expected opt-in reproducibility skip;
- source package build: passed;
- `R CMD check --no-manual --as-cran`: installation, loading, compiled-code,
  examples, R/Rd checks, and static checks passed, but the overall check ended
  with the pre-existing installed-package test error in
  `test-coupling-tempering-offline-ratios.R`, which sources the intentionally
  unshipped `tools/study06_partial_exchange_feasibility.R`; status was 1 error
  and 5 notes;
- `git diff --check`: passed;
- modified native sampler or Rcpp files: none;
- generated tarball, check directory, object files, and DLL: removed.

The compact ignored evidence has SHA-256
`D97CDD4814AD3F6872E17579E0944180ED99003E9F29FCDB5FDFEFC3CF432F44`
for `oracle_result.json` and
`A0CB20C8A23D2563854E390114538806CDADF0C8DBE0AB9C05B7F14B95B32152`
for `oracle_chains.rds`.
