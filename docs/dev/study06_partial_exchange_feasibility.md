# Study 06 partial-exchange feasibility audit

## Decision

This offline audit selects **F6: retained state is insufficient**. The completed
packed-BED screen retained the full scientific history only for the `lambda = 1`
slot. It did not retain the allocations, alpha, `sigmaSqAlpha`, effects,
residuals, or variance state of the `lambda = 0` and `lambda = 0.5` replicas at
exchange attempts. Therefore neither the native complete-state exchange ratio
nor any proposed exact partial exchange can be reconstructed independently.

No MCMC was run, no sampler source was changed, and the existing T4 decision was
not altered.

## Provenance and preservation

The source branch remained `master` at
`8908267a68a46267fcccb910850b4f6380bfa978`. The read-only `sblrbench` repository
remained clean at `de31f62e182d8540488d4135df4c58f052a515d9`.

The Study 06 identities were revalidated from the retained provenance:

- specification hash: `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`;
- truth hash: `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.

Before this audit the complete tracked unstaged diff was saved to
`results/local/study06_bed_coupling_tempering_screen/pre_partial_exchange_audit.patch`.
Its SHA-256 is
`b425a5f20c757bc71d097740ea30c7c0b33cc00dd31312e8436adf2e3d5d5647`.
The ignored `preserved_file_hashes.csv` records SHA-256 hashes for all eight
modified tracked implementation files, all five untracked tempering files, and
the nine retained evidence files used here.

## Retained-state inventory

For each of four ensembles the formatted checkpoint contains:

| Quantity | Scope | Dimensions |
|---|---|---:|
| component allocations | `lambda = 1` only | 3000 x 1500 |
| marker effects | `lambda = 1` only | 3000 x 1500 |
| active indicators | `lambda = 1` only | 3000 x 1500 |
| alpha | `lambda = 1` only | 3000 x 12 |
| `sigmaSqAlpha` | `lambda = 1` only | 3000 x 3 |
| replica identity | all three slots | 3000 x 3 |
| active count | all three slots | 3000 x 3 |
| expected active count | all three slots | 3000 x 3 |
| swap record | each attempted adjacent pair | 600 x 9 |

The swap record contains iteration, adjacent-pair index, accepted flag,
acceptance probability, complete log ratio, pre-swap identities, and target
active count before and after the attempt.

The checkpoint does **not** contain per-iteration lower-replica components,
effects, residuals, alpha, `sigmaSqAlpha`, marker probabilities, or RNG state.
Final target summaries do not repair this loss because all swaps were rejected
and the missing replicas evolved independently.

## Complete exchange ratio

For adjacent levels $a,b$, the implemented complete-state ratio is

\[
R_c=\ell_a(c_b;\alpha_b)+\ell_b(c_a;\alpha_a)
-\ell_a(c_a;\alpha_a)-\ell_b(c_b;\alpha_b),
\]

where $\ell_\lambda(c;\alpha)=\log p_\lambda(c\mid\alpha)$.
Independent recomputation requires both component vectors and both alpha matrices
at every attempt. Only the target vector and alpha were retained. Consequently
the native/offline error, marker contributions, stick contributions, and
floating-point mismatch count are **not estimable**. Treating active count as an
allocation vector or reconstructing alpha from expected active count would not be
an exact calculation.

This is a state-retention failure for the requested audit, not evidence against
the already validated tiny-model exchange algebra.

## Exact partial-exchange derivations

Let $h(\alpha,\sigma)$ denote the proper hierarchy density, including normal
intercept priors, normal non-intercept coefficients conditional on
`sigmaSqAlpha`, and the scaled-inverse-chi-square variance prior. Let
$\alpha_a^{B\leftarrow b}$ replace block $B$ of alpha in replica $a$ with
the corresponding block from replica $b$.

### P1: all alpha and `sigmaSqAlpha`

\[
R_{P1}=\ell_a(c_a;\alpha_b)+\ell_b(c_b;\alpha_a)
-\ell_a(c_a;\alpha_a)-\ell_b(c_b;\alpha_b).
\]

The two $h(\alpha,\sigma)$ terms travel together and cancel pairwise because
the hierarchy prior is common across coupling levels. Exact evaluation needs
both lower-replica allocations, alpha, and variances; these are missing.

### P2: all alpha only

P2 adds

\[
\log p(\alpha_b\mid\sigma_a)+\log p(\alpha_a\mid\sigma_b)
-\log p(\alpha_a\mid\sigma_a)-\log p(\alpha_b\mid\sigma_b)
\]

to the P1 allocation term. Proper intercept-prior products cancel, while the
non-intercept normal densities under destination variances do not. Both lower
alpha and variances are missing.

### P3 and P7: complete stick hierarchy

For a selected stick $j$, replace alpha column $j$ and its variance in each
replica. The ratio is the allocation term computed with hybrid alpha matrices.
The stick-specific alpha/variance prior travels as one block and cancels. P3 uses
stick 1; P7 separately uses sticks 2 and 3. Required lower alpha, variance, and
allocation states are missing.

### P4: first-stick alpha only

The allocation term uses hybrid first-stick alpha. Add cross-evaluation of the
three non-intercept coefficients under the destination first-stick variances.
The intercept prior cancels. Required lower state is missing.

### P5: first-stick intercept only

Only the allocation term changes. The identical proper intercept-prior products
cancel under exchange. Both first-stick intercepts and allocations are required;
the lower values are missing.

### P6: first-stick non-intercept effects only

Use hybrid enriched, continuous, and null-annotation coefficients in the
allocation term and add their conditional-normal prior cross-ratio under retained
first-stick variances. Required lower state is missing.

### P8: complete non-annotation state

If $x=(c,\beta,r,v_b,v_e,\ldots)$ is swapped as a complete internally
consistent block, phenotype likelihood, effect prior, and common variance priors
cancel globally. The remaining ratio evaluates each moved allocation under the
destination coupling level and the destination replica's retained alpha. There
is no Jacobian for a permutation. Exact calculation requires both replicas'
complete marker, residual, and variance states; only the target state exists.

### P9: allocation only

An allocation-only permutation while effects remain fixed generally crosses the
spike-and-slab support: a nonzero effect assigned to the null component has zero
target density. The maintained residual remains a function of beta, not the
component label, but the state is usually outside the component-conditioned
effect-prior support. This is not a useful exchange block. In any case the lower
allocations and effects are absent.

## What the retained scalar evidence does show

The all-slot active and expected-active counts can be aligned exactly to the 600
recorded attempts per ensemble.

| Pair | Median active-count difference | Median expected-count difference | Median log ratio | Mean log ratio | Range |
|---|---:|---:|---:|---:|---:|
| 0--0.5 | 48.0 | 57.6 | -348.2 | -613.6 | -20615.6 to -109.5 |
| 0.5--1 | -7.0 | -8.1 | -454.6 | -461.5 | -889.1 to -128.0 |

For 0--0.5, the absolute active-count difference correlated only -0.097 with the
log ratio. For 0.5--1 the correlation was -0.377. Thus aggregate sparsity
mismatch contributes, especially at the lower bridge, but cannot by itself
explain the separation: the 0.5 and 1 slots had similar mean counts while their
median exchange penalty was even more negative. Distinguishing marker identity,
component allocation, and alpha-stick contributions requires the missing vectors.

This evidence supports only the narrow statement that **global sparsity is not
the sole cause**. It cannot identify first-stick versus later-stick alpha,
intercepts versus annotation effects, `sigmaSqAlpha`, causal/enriched marker
groups, or diffuse markerwise accumulation.

## Requested decompositions and bridge audits

The following are unavailable without approximation:

- marker, stick, component, null/active, causal, enrichment, and continuous-bin
  decomposition;
- positive/negative marker concentration and 50/80/90/95% concentration counts;
- P1--P8 hypothetical acceptance distributions;
- direction-specific partial ratios;
- sparsity-matched BayesR, target-mean, or oracle endpoints;
- mean-preserving predictor, continuation-probability, or component-probability
  paths;
- stratified marker-count scaling;
- active-set Jaccard, causal/enriched overlap, marker-prior correlation, and alpha
  predictor correlation.

Candidate baseline proportions can be defined from external summaries, but an
exchange overlap statistic still requires the missing source allocation and
alpha states. Reporting a one-state self-evaluation would not estimate exchange
acceptance and was intentionally not done.

## Minimal state for a future offline audit

A later **short diagnostic run**, if separately authorized, would need to retain
at each exchange attempt and for every slot:

1. component vector and alpha matrix;
2. all three `sigmaSqAlpha` values;
3. marker probabilities or enough information to reproduce them exactly;
4. for P8, marker effects, maintained residual, effect variance, residual
   variance, and every other persistent variance state;
5. immutable marker/annotation identities and causal/enrichment labels can remain
   external because they are already hashed.

Only exchange-attempt checkpoints are necessary; full per-iteration lower-replica
histories are not. Storing log allocation contributions by marker and stick at
attempt time would make decomposition cheaper, but raw states are preferable for
independent verification.

## Recommendation

No exchange proposal deserves implementation from the current evidence. The
correct next task is a narrowly designed **state-retention review** specifying a
compact exchange-attempt checkpoint schema and storage estimate. That review
must precede any further short MCMC and must not presume that partial exchange or
a matched endpoint will work.

The existing T4 conclusion remains unchanged. F5 cannot be selected because the
partial and matched-endpoint ratios were not observable; F6 is the only valid
decision.

## Offline outputs

The reproducible analysis is `research/sbayesrc/tools/study06_partial_exchange_feasibility.R`.
Ignored outputs under `partial_exchange_feasibility/` include the checkpoint
inventory, retained count compatibility, count summaries/correlations, proposal
feasibility table, SHA-256 manifest, and local machine-readable decision.
