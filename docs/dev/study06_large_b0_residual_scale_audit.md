# Study 06 large B0 residual-scale audit

## Decision

**B0-C4: unresolved.** The iteration-0 failure is genuine and is not numerical
cancellation. No tolerance, clamp, rank substitution, or posterior-changing
fallback was added.

The frozen experiment identities were:

```text
specification b001bc36a5531e5e6b342286a253fc1fd34dad4265359d89d2feaa026d4533df
truth         e94a511540f600e61ef47b52947836f19a15388f5e8ce795c179929956817507
markers       37,991
blocks        76 (75 x 500, 1 x 491)
retained rank 37,991
```

## Reproduction

The exact full-positive B0 state fails after its first marker/allocation sweep,
before the residual-variance draw:

| Term | Value |
|---|---:|
| `yy` | 10,223.891131 |
| initial residual variance | 1.022594 |
| base effect variance | 0.005909 |
| occupied markers | 1,945 |
| effect squared norm | 2.675477 |
| transformed-score squared norm | 82,949.69 |
| maintained reduced-residual squared norm | 68,502.41 |
| score product | 13,929.360314 |
| block fitted quadratic | 13,411.437722 |
| block SSE | -4,223.391775 |
| prior contribution | 2.045187 |
| residual scale | -4,221.346588 |

The expanded native error reports these inputs without changing the execution
path. An exact one-iteration diagnostic with `updateE = FALSE` captured the
state immediately before the failed draw. That diagnostic makes no residual
RNG call, so its state agrees with the failing trajectory through the point of
failure.

## Independent oracle

An independent R oracle reconstructed every retained factor, `w_b`,
`r_b = w_b - Q_b beta_b`, and block contribution. It reproduced the native SSE
within `4.04e-10`. Direct prediction from the same frozen BED genotypes gave:

| Oracle term | Value |
|---|---:|
| direct score product | 13,929.36 |
| direct full-genotype quadratic | 30,740.09 |
| direct SSE | 13,105.26 |
| direct residual scale | 13,107.31 |
| quadratic omitted by block factorization | 17,328.66 |

The score products agree within `2.20e-6`; the discrepancy is the omitted
cross-block quadratic. Adding it maps the invalid block SSE to the positive
direct SSE.

The native runtime diagonal spans 4,714.349 to 5,292.417. All 37,991 positive
modes were retained. The required diagnostic 0.995 run also failed materially
(`SSE = -3,950.014513`), and setting `adjE = 0` did not change the first failure.
Neither is a scientific substitute.

## Interpretation

The retained operator represents within-block cross-products. Retaining every
positive within-block mode does not restore cross-block cross-products. For
this 76-block same-sample panel, concatenated block projections are not one
orthogonal subspace of the phenotype. Consequently the documented global
projected residual identity

```text
yy - sum(w_b'w_b) + sum(r_b'r_b)
```

is not a nonnegative residual norm for this state. A scale-aware tolerance is
inapplicable because the deficit is thousands of units, not rounding error.

Possible changes—using the direct BED quadratic, adopting a reduced-coordinate
or GCTB-style block variance conditional, adding cross-block factors, or
changing blocks/rank—define a different likelihood or residual-variance
contract. None can be introduced while claiming the frozen experiment and
posterior target are unchanged. A separately reviewed contract decision is
required before B0/B1/B2 can run.

The machine-readable local audit is
`results/local/06_annotation_models/large_feasibility/continuation/b0_residual_scale_audit.json`;
block terms are beside it in `b0_block_terms.csv`.

