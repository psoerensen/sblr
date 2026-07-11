# Block-Eigen LD Operator

This note describes the experimental block-eigen LD backend infrastructure for
CSR summary-statistic models. Phase 1 adds reusable C++ infrastructure only; it
does not route production fitting through the new operator.

## Purpose

The block-eigen operator is intended to support eigenvalue filtering and
projection workflows used by SBayesRC-style methods while keeping the
implementation general for CSR models. It is an LD operator backend, not
SBayesRC-only code, because ordinary BayesC, BayesR, annotation-informed BayesC,
and SBayesRC all need the same sampler-facing LD operations.

## Operator Contract

CSR samplers should eventually depend on an LD operator abstraction with three
methods:

```text
diag()
apply_offdiag()
rebuild()
```

`diag()` returns the self terms for the active operator. Existing sparse CSR LD
uses the original `xx`. A hard-truncated block-eigen operator uses
`diag(tilde_A)`, which can differ from the original `xx`.

`apply_offdiag(i, diff, r)` applies the off-diagonal residual update for a
single marker effect change. `rebuild(wy, b, r)` rebuilds residuals as
`r = wy - A b`.

## Filtering And Storage

Block construction starts from exact standardized BED blocks. Each block forms
the exact PSD crossproduct

```text
A = X_S' X_S
```

in `X'X` units, then converts to correlation scale:

```text
C = D^{-1/2} A D^{-1/2}
```

Eigenvalue filtering and shrinkage thresholds are applied in correlation space.
The filtered result is then mapped back to `X'X` units and stored by the
operator. This keeps sampler residual semantics consistent with existing CSR
code.

The block-eigen operator stores each dense block as a packed float upper
triangle. Block membership and local marker indices are stored separately, and
the operator keeps a precomputed `diag_` vector for O(1) diagonal access.

## Projection Rules

For hard truncation, eigenvalues below `max(tau, 0.01)` in correlation space are
dropped, while at least the largest eigenvector direction is retained. The
operator stores

```text
tilde_A = D^{1/2} V_k diag(mu_k) V_k' D^{1/2}
```

and `wy` is projected in place for every trait using the same retained
subspace:

```text
W_block <- (W_block * Lk) * Rk
Lk = D^{-1/2} V_k
Rk = V_k' D^{1/2}
```

For fixed ridge and Ledoit-Wolf shrinkage, the operator remains full-rank:

```text
tilde_A = (1 - a) A + a D
```

No `wy` projection is performed for those modes.

## Future Integration Order

The intended integration order is:

1. ordinary CSR BayesC
2. CSR BayesR
3. CSR SBayesRC
4. CSR prior/group/learned annotation BayesC

## Phase 2A Status

- ordinary CSR BayesC now uses `CsrOperator` for the existing sparse-CSR path
- block-eigen is still infrastructure only and not yet exposed
- no public R API changes
