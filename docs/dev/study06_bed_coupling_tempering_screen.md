# Study 06 packed-BED coupling-tempering screen

## Status and scope

This is a package-development proof of concept, not Study 06 qualification
evidence. It used packed BED only, ran 3,000 recorded iterations per ensemble,
and did not run CSR, block eigen, the formal qualification pipeline, or the
final benchmark. The starting package commit was
`8908267a68a46267fcccb910850b4f6380bfa978` and the read-only `sblrbench`
commit was `de31f62e182d8540488d4135df4c58f052a515d9`.

The immutable Study 06 informative identity was:

- specification hash: `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`;
- truth hash: `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`;
- 2,000 individuals, 1,400 training and 600 validation;
- 1,500 markers and realized components 1,329 / 84 / 50 / 37;
- realized heritability 0.50 and `gamma = c(0, 0.01, 0.1, 1)`.

## Target family

Let (c_i) be marker (i)'s component, (A_i=(1,A_{i,-0})), and let
(mu_k) be the registered baseline probit-stick intercept. At coupling level
(lambda\in[0,1]), the continuation predictor is

\[
\eta_{ik}^{(\lambda)}=(1-\lambda)\mu_k+
\lambda A_i\alpha_k.
\]

The existing stick-breaking map converts the finite continuation probabilities
to a normalized component-probability vector. The existing probability floor
and renormalization are applied identically at every level. At `lambda = 1`,
this is exactly the current BayesRC allocation prior. At `lambda = 0`, the
allocation prior is the fixed registered BayesR-like mixture and is independent
of alpha. The proper intercept prior, non-intercept hierarchy, effect prior,
likelihood, and variance priors are unchanged.

For the Albert--Chib step, the latent-normal mean is an offset
((1-\lambda)\mu_k) plus design (lambda A_i). Consequently the existing
Gaussian coefficient update remains valid after replacing each design entry by
(lambda A_{iq}). At zero coupling the likelihood precision for alpha is zero,
so alpha is updated from its proper hierarchical prior; `sigmaSqAlpha` receives
its usual conditional update. Later-stick eligibility continues to be determined
only by the current component allocation.

## Replica exchange

Each replica owns marker effects and components, the maintained phenotype
residual, marker/effect and residual variances, derived genetic variance and
heritability state, alpha, `sigmaSqAlpha`, marker probabilities, logical replica
identity, and its RNG stream. Latent probit variables are not persistent and are
regenerated at the next hierarchy update.

For adjacent levels (a,b), and persistent states (S_a,S_b), the exact log
exchange ratio is

\[
\log r = \log p_a(c_b\mid\alpha_b)+\log p_b(c_a\mid\alpha_a)
-\log p_a(c_a\mid\alpha_a)-\log p_b(c_b\mid\alpha_b).
\]

The phenotype likelihood, marker-effect priors, variance priors, alpha prior,
and `sigmaSqAlpha` prior cancel because they are identical across levels. The
calculation uses the complete marker allocation probability in log space. After
an accepted swap, marker probabilities are recomputed at the destination level.
The alternating adjacent-pair schedule is fixed and state independent.

## Implementation boundary

The proof of concept is reached only through the development option
`sblr.development.bed_coupling_tempering`, which appends an unambiguous native
control to the already resolved intercept-prior matrix. Its absence leaves the
ordinary path and public signature unchanged. The ladder is fixed at
`c(0, 0.5, 1)`, exchange is attempted every five complete 1-allocation / 1-hierarchy
cycles, and only the level-one state is retained. Compact coupling diagnostics
are carried through the existing per-chain diagnostic object only for tempered
runs.

## Tiny near-exact validation

The independent reference used 8 markers, 12 samples, one binary non-intercept
annotation, two components, one stick, fixed effect variance 0.20, and fixed
residual variance 1. It enumerated all 256 component allocations, integrated
marker effects analytically, and integrated the two alpha coefficients by
independent Gaussian quadrature. Tolerances were declared before sampling.

| Check | Error | Declared tolerance | Result |
|---|---:|---:|---|
| PIP, quadrature order 31 vs 41 | 1.34e-9 | 2e-4 | pass |
| Active count, quadrature order 31 vs 41 | 6.00e-9 | 5e-4 | pass |
| PIP, tempered target vs order-41 reference | 0.00266 | 0.035 | pass |
| Active count, tempered target vs reference | 0.00654 | 0.15 | pass |
| Alpha mean, tempered target vs reference | 0.00731 | 0.18 | pass |

Four 50,000-iteration tiny runs (10,000 burn-in) were used. Mean exchange rates
were 0.510 for 0--0.5 and 0.647 for 0.5--1. Algebraic reverse-ratio and numerical
detailed-balance tests passed. Endpoint probability tests passed, and an explicit
disabled control reproduced the ordinary native object exactly apart from no
timing field being present. Repeated tempered runs reproduced all scientific
draws and diagnostics exactly; only wall-clock timing values differ.

## Short Study 06 screen

Four ensembles used registered seeds 701121, 701222, 701323, and 701424. Each
contained levels 0, 0.5, and 1; used 3,000 recorded cycles, 1,000 burn-in cycles,
and 2,000 retained target draws; and attempted a fixed alternating adjacent swap
every five cycles.

### Exchange failure

Every ensemble made 300 attempts for each adjacent pair. All 2,400 attempts were
rejected and no replica identity ever left its initial level. Across ensembles:

- 0--0.5 maximum exchange probabilities ranged from 2.73e-48 to 3.21e-53;
- 0.5--1 maximum exchange probabilities ranged from 2.59e-56 to 3.76e-78;
- median log ratios ranged from -331.5 to -373.2 for 0--0.5;
- median log ratios ranged from -441.0 to -480.3 for 0.5--1;
- complete round trips: zero;
- longest period without a successful exchange: all 3,000 cycles.

These are finite exact log ratios. The result is target separation, not a
floating-point or acceptance-guard failure.

### Active-count behavior

| Ensemble | lambda 0 mean (range) | lambda 0.5 mean (range) | lambda 1 mean (range) | target lag-50 ACF |
|---:|---:|---:|---:|---:|
| 1 | 11.0 (7--17) | 58.6 (25--119) | 59.3 (23--164) | 0.532 |
| 2 | 11.0 (7--16) | 65.5 (24--144) | 50.7 (21--128) | 0.063 |
| 3 | 10.9 (7--17) | 69.2 (25--152) | 46.7 (20--105) | 0.100 |
| 4 | 11.0 (7--17) | 56.6 (24--119) | 57.6 (19--161) | 0.483 |

The zero-coupling state was confined to the 0--49 band. Target chains did visit
overlapping regions, but no visit was caused by exchange. The target active-count
R-hat was 1.074, bulk ESS 74.8, tail ESS 73.6, and relative MCSE 0.132. Expected
active count was similar (R-hat 1.074, bulk ESS 75.1). Component-1 occupancy was
especially unstable (R-hat 1.185, bulk ESS 20.7). Target lag-50 autocorrelation
was heterogeneous rather than consistently improved. The committed BED S1 and
H20 chain-1 comparators were 0.390 and 0.147, respectively.

### Annotation and variance summaries

The short target screen did not converge. Alpha R-hat ranged from 1.014 to 1.444;
the minimum alpha bulk ESS was 8.08. `sigmaSqAlpha` R-hat was 1.020, 1.319, and
1.021 by stick, with stick-2 reaching 872.5. Effect variance also remained slow
(R-hat 1.028, bulk ESS 133). Genetic variance, residual variance, and heritability
were stable: their R-hats were 1.001, 1.002, and 1.002, with mean heritability
0.416.

### Scientific summaries

The short target-only posterior retained the known annotation-ranking benefit:

- PIP AUPRC 0.590 and AUROC 0.848;
- recall / precision at 10: 0.0585 / 1.00;
- at 25: 0.140 / 0.96;
- at 50: 0.240 / 0.82;
- at 100: 0.404 / 0.69;
- Bayesian FDR 5%: 19 selected, 19 true;
- Bayesian FDR 10%: 23 selected, 22 true;
- posterior-effect correlation with truth 0.913;
- validation genetic-value correlation 0.917;
- phenotype prediction correlation 0.654.

These are close to committed BED S1/H20 results (AUPRC about 0.595, AUROC
0.852--0.853, validation correlation 0.918, prediction correlation 0.655).
Scientific agreement does not rescue an exchange mechanism that never moved.

## Runtime and decision

Wall time was 317.4 seconds for 36,000 replica marker sweeps. This is 0.0265
seconds per recorded target draw across ensembles. Native accumulated transition
time was 798.4 worker-seconds and exchange evaluation used 24.9 worker-seconds.
The ignored local evidence occupies the fit/checkpoint directory; raw chains are
not tracked.

Decision **T4 (exchange mechanism fails)** is selected. The three-level ladder
and this coupling path are unusable as configured: there were no accepted swaps,
no round trips, and therefore no possible persistent target-state change caused
by exchange. A full-length run is not recommended. Any denser-ladder proposal
requires a separate preregistered design task; it must first address the extensive
allocation-prior log-density gap at 1,500 markers. No package default or public
model contract should change from this screen.

## Reproducibility and limitations

The implementation and tiny reference are in
`research/sbayesrc/tools/coupling_tempering_tiny_reference.R`; the Study 06 harness is
`research/sbayesrc/tools/study06_bed_coupling_tempering_screen.R`. Large local evidence is ignored.
This screen is intentionally too short for qualification, tests one fixed ladder,
and says nothing about CSR, block eigen, substantial eigen truncation, or the
separate BED/block heritability calibration issue.
