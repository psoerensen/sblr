# Study 06 scalar SBayesRC stabilization

> **Status: RESEARCH/EXPERIMENTAL RECORD.** This document preserves historical
> Study06 diagnoses and corrections at its recorded commits. It is not current
> general architecture authority; see [`README.md`](README.md).

## Status

Study 06 is **not ready for four-entry requalification**. The CSR failures are
classified as unsupported hard-sparse likelihood states rather than residual
drift, and retained SBayesRC plus the annotation tail sampler received focused
corrections. The preserved BED histories, however, demonstrate complete
separation under the configured flat stick intercept. Resolving that requires a
focused prior/model decision and a new package SHA; it cannot be corrected by a
posterior-preserving mixing change alone.

## Provenance

- `sblr`: starting HEAD `02e8c74baa906e83c4a08d42a9cc6339b4e81072`,
  *Align BLR prior variance calibration*, version 0.2.0.
- `sblrbench`: read-only HEAD
  `964bf188fc990d5fbbad0ea00592b288c4442723`.
- Study 06 recorded package source `02088ea8`; result evidence `41ed452`;
  integration history also records `71f627`.
- Study 05 and Study 06 use `sblr` SHA `02e8c74...`; qgdata SHA
  `6cca5819e711d326cfb2614d7e9d9f34942612cd`.
- Primary external sources and commits are recorded in
  `sbayesrc_reference_crosswalk.md`.

The authoritative Study 06 identities were revalidated: 2,000 samples,
37,991 markers, 1,400 training samples, sample hashes
`9d079...`/`f6a556...`, marker hash `c39ca...`, annotation hash `4a0e...`,
informative phenotype hash `437a...`, uninformative phenotype hash `77035...`,
and the recorded sparse LD prefix with 8,888,984 upper-triangle entries,
`r2 >= 0.001`, and a 1,000-marker window.

## Study 05 evidence

Study 05 already demonstrated the same mechanism at lower severity. Exact CSR
closely reproduced full-sweep BED BayesR. The tested hard-sparse matrix was
positive definite (minimum eigenvalue 2822.033 in cross-product units), yet its
relative Frobenius error was 0.4421812 and it changed corrected scores,
conditional mixture probabilities, quadratic forms, component allocation, and
posterior recovery. Exact-minus-sparse quadratic errors were 39.44685 for the
truth, 32.20268 for the BED state, 32.17449 for the exact-CSR state, and
62.69658 for the sparse-CSR state. Retaining 99.5% of within-block positive
spectral mass improved some checks but could not restore omitted weak or
cross-block LD.

## Authoritative CSR reproduction and residual algebra

Both original four-chain qualification calls were reproduced before changing
sampler behavior. Informative replicate 1 failed logical chain 1 (one-based
chain 2); uninformative replicate 1 failed logical chain 0 (one-based chain 1).
Exact-seed diagnostic replays then localized the first invalid transitions:

| Scenario | Seed | Iteration | Sparse maintained scale | Sparse rebuilt scale | Sparse quadratic scale | Exact counterfactual scale |
|---|---:|---:|---:|---:|---:|---:|
| Informative | 701242 | 6073 | -169.254890930703 | -169.25489093070345 | -169.25489093070573 | 1715.927618303926 |
| Uninformative | 801141 | 3983 | -83.32739554486146 | -83.32739554486214 | -83.32739554486328 | 1715.242011175886 |

The maintained, rebuilt, and independently evaluated sparse identities agree
to numerical precision. Maximum residual drift was
`3.7801e-12` (relative `4.4438e-16`) for informative and `1.2790e-12`
(relative `1.5342e-16`) for uninformative. Scaling, alignment, and finite-state
checks passed. Residual-state drift and an invalid algebraic identity are not
the cause.

At the informative state, `b'R_sparse b=1986.860277678945` and
`b'R_exact b=3872.042786913567`, giving `delta_q=1885.182509234623`.
At the uninformative state the corresponding values are
`1898.753532133455`, `3697.322938854203`, and `1798.569406720748`.
Holding every non-operator quantity fixed therefore changes each invalid sparse
scale into a strongly positive exact scale.

Unlike the Study 05 sparse operator, the exact Study 06 hard-sparse matrix is
also indefinite. Sparse Cholesky factorization of its 37,991 by 37,991
operator (8,888,984 stored upper entries, implicit unit diagonal) failed with
`leading principal minor of order 881 is not positive`. This is an integrity
failure independent of residual bookkeeping. The observed failure-state
quadratics remain positive, so indefiniteness does not replace the direct
counterfactual: loss of quadratic fidelity explains the negative scales, while
indefiniteness can additionally distort sampler trajectories. No PSD projection
or operator substitution is performed.

## Omitted-LD and blow-up diagnosis

The active-marker-only decomposition did not materialize genome-wide LD:

| Scenario | Active markers | Diagonal contribution | Retained pair contribution | Omitted pair contribution | Omitted positive | Omitted negative |
|---|---:|---:|---:|---:|---:|---:|
| Informative | 2418 | 2003.198933 | -16.338655 | 1885.182509 | 24816.047362 | -22930.864853 |
| Uninformative | 569 | 1900.018214 | -1.264682 | 1798.569407 | 7369.297506 | -5570.728100 |

No active high-LD (`|r| >= 0.8`) opposite-signed pair explained either
failure. Many omitted below-threshold and outside-window terms have large
cancelling positive and negative totals, with a large positive remainder. This
is loss of quadratic fidelity from hard sparsification.

The informative trajectory also resembles GCTB's documented collective effect
blow-up: residual variance fell to 0.02383, hard-sparse heritability reached
0.9835, marker variance rose to 0.0390, and active count expanded to 2,418.
Maximum absolute effect remained only 0.424, so the failure is collective rather
than one extreme SNP pair. The uninformative state had heritability 0.9038 and
569 active markers. Hard-sparse LD both permits this state and understates its
quadratic penalty.

## Retained block eigen

The retained route implements the projected identity
`yy_projected - 2 b'w_projected + b'Qb` with retained positive eigenvectors and
one global projected residual variance. It does not silently use block-specific
fallbacks. During this audit, a real control-plumbing defect was corrected:
scalar SBayesRC previously documented `low_rank_residual_rebuild_every` but did
not pass it to native code. Periodic and final reduced-residual rebuilds now
record rebuild counts and maximum absolute drift, consistent with scalar
SBayesR. Deterministic full-positive-rank reduction tests compare the route with
its configured block target.

A labelled non-qualification control used the informative Study 06 data,
annotations, priors, initialization, summary statistics, seed 701242, 0.995
retention, and 50 iterations. It completed with residual variance in
`[0.7368855, 1.080391]`, minimum genetic variance `0.3626268`, 51 projected
residual rebuilds, and maximum drift `1.554312e-15`. Protected sibling inputs
were byte-identical afterward. This establishes numerical validity at the
tested control point; it is not qualification evidence.

This route remains an approximation: it omits cross-block LD and discarded
spectral mass, and its global residual-variance contract differs from current
GCTB. It is the canonical scalable route, not a universal replacement for exact
CSR or BED.

## BED annotation mixing

The preserved 9,000-by-4 histories were analyzed without rerunning them. The
qualification failure lists remain those frozen by Study 06: 17 of 23
quantities for informative and 21 of 23 for uninformative. Full-history
diagnostics confirm broad chain-location differences and very slow movement.
Examples:

- informative component-2 stick intercept: chain means 9.20--48.96, maxima
  40.58--96.82, lag-100 autocorrelation 0.907--0.962;
- uninformative component-1 stick intercept: means 15.92--28.26, maxima
  47.23--72.66, lag-100 autocorrelation 0.915--0.968;
- uninformative component-2 stick intercept: means 20.24--45.77 and a maximum
  of 110.09;
- informative worst full-history R-hat 1.902 with bulk ESS 5.64;
  uninformative worst R-hat 1.724 with bulk ESS 6.15.

Final component states show the mechanism. Informative chains 3 and 4 contain
only components 0 and 3; every uninformative chain contains only components 0
and 3. Component 2 is empty in all eight final states, and component 1 is empty
in six. Later eligible stick sets are nonempty but contain only successes, so
the probit intercept is completely separated. Under `intercept_flat=TRUE`, its
conditional posterior has an unbounded positive tail. This creates the extreme
intercepts, strong dependence on component occupancy/effect variance, and
derived-prior instability. It is not an initialization, annotation ordering,
or empty-eligible-set defect.

The shared latent sampler also had a separate numerical defect: inverse-CDF
probabilities clipped at `1e-12` could return a draw violating its truncation in
extreme tails. An exact normal/exponential rejection sampler now preserves
truncation and the intended posterior. That correction prevents invalid latent
draws but cannot make a flat-intercept separated posterior proper.

No manuscript/GCTB empty-stick fallback was copied. Setting a coefficient to
`-10` is not a draw from the configured posterior, and the Study 06 case has
nonempty eligible sets. A blocked alpha update would preserve the same separated
target and cannot resolve the blocker. A proper intercept prior or another
identified model constraint is the required next method-development decision;
it would change the Study 06 prior and therefore is outside this stabilization.

## Implemented contracts and diagnostics

- Opt-in native failure capture reports chain/iteration, all residual
  identities, quadratics, variance state, component counts, finite checks,
  residual drift, vectors, and preceding trajectories. Normal execution pays
  only an environment-flag check.
- Invalid scales still fail. There is no clamp, absolute value, zero
  replacement, previous-value fallback, operator switch, SNP removal, or PSD
  projection.
- CSR fits now classify recorded complete/full construction,
  recorded hard-sparse construction, and unclassified external CSR separately.
  Construction metadata is not mislabeled as source-fidelity evidence.
- BED is labeled an individual-level reference; retained low rank is labeled a
  canonical scalable summary approximation; reconstructed dense block eigen is
  retained explicitly for regression/reproducibility.
- The development tools reproduce exact Study 06 inputs, protect sibling files
  with hashes, write outside qualification directories, and avoid dense
  genome-wide LD.
- Deterministic tests cover extreme two-sided latent tails, operator
  reductions, retained residual rebuilding, metadata classification, and the
  contextual no-clamping invalid-scale path.

After the posterior-preserving tail correction, the exact authoritative chain
seeds were rerun with the unchanged 9,000-iteration ceiling. Informative seed
701242 failed at iteration 1,047 (sparse scale `-183.148944`, exact scale
`1968.405958`, `delta_q=2151.554902`, drift `9.81e-13`). Uninformative seed
801141 passed its original failure point but failed at iteration 5,913 (sparse
scale `-250.215723`, exact scale `1938.738074`, `delta_q=2188.953797`, drift
`1.85e-12`). Changed iterations reflect a corrected Markov trajectory; both
reruns confirm that the hard-sparse route remains unsupported for these states.

## Remaining gate

Validation on R 4.4.1/Rtools44 GCC 13.2.0 completed native compilation,
`devtools::load_all(quiet=TRUE)`, 627 focused assertions, and the full
`devtools::test()` suite: 3,908 passed, zero failed, one pre-existing covariance
warning, and one opt-in extended-reproducibility skip.

The CSR root cause is resolved as an unsupported hard-sparse operator/input
contract, not repaired by changing its posterior. Retained identities and drift
tests and the numerical annotation correction are package-testable. BED mixing
remains blocked by flat-intercept complete separation, so Study 06 must not be
declared complete or rerun for qualification until that model/prior decision is
made, implemented under a new package SHA, and all four qualification entries
are rerun unchanged.
