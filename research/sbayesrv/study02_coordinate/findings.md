# Study 02-coordinate findings

Status: `SBayesRV_PREDICTION_SIGNAL_OBSERVED_AT_STUDY02_COORDINATE`

Separate limitation: `ORDINARY_SBayesR_SPARSE_LD_BOUNDARY_AT_Q500_H050`

## Provenance and truth

This one replicate reused qgdata `human_independent`, 37,991 post-QC markers,
the validated 3,500/1,500 Study 02 split, training-only allele frequencies,
and the existing CSR LD operator. The shared annotation surface had proxy
correlation 0.697555, true-$q$ range 0.198613--8.764821, median 1.002801, and
geometric mean exactly one to numerical precision. These results are
descriptive research output, not accepted Study 02 evidence.

## Held-out and marker-level results

Errors use lower-is-better orientation. Values in parentheses are
`SBayesRV - ordinary SBayesR`.

| stage and metric | ordinary SBayesR | SBayesRV | paired difference |
|---|---:|---:|---:|
| `q50_h030` prediction correlation | 0.898962 | 0.903341 | +0.004380 |
| `q50_h030` prediction NMSE | 0.192086 | 0.184373 | -0.007713 |
| `q50_h030` effect RMSE | 0.00150353 | 0.00146464 | -0.00003889 |
| `q50_h030` PIP Brier score | 0.00117717 | 0.00115194 | -0.00002523 |
| `q50_h030` average precision | 0.297408 | 0.293791 | -0.003617 |
| `q500_h030` prediction correlation | 0.837002 | 0.842821 | +0.005819 |
| `q500_h030` prediction NMSE | 0.300387 | 0.290225 | -0.010162 |
| `q500_h030` effect RMSE | 0.00184355 | 0.00180252 | -0.00004104 |
| `q500_h030` PIP Brier score | 0.01273607 | 0.01271464 | -0.00002144 |
| `q500_h030` average precision | 0.052104 | 0.050888 | -0.001216 |

Phenotype prediction correlation was 0.491792 versus 0.489186 at `q50_h030`
and 0.487935 versus 0.491145 at `q500_h030`. It is secondary because the
primary prediction target was held-out genetic value.

| stage | method | top 50 | top 100 | top 500 | top 1000 |
|---|---|---:|---:|---:|---:|
| `q50_h030` | ordinary | 0.30 | 0.34 | 0.44 | 0.54 |
| `q50_h030` | SBayesRV | 0.32 | 0.34 | 0.46 | 0.50 |
| `q500_h030` | ordinary | 0.038 | 0.046 | 0.070 | 0.084 |
| `q500_h030` | SBayesRV | 0.038 | 0.044 | 0.064 | 0.078 |

Thus prediction and effect recovery improved modestly in both completed
stages, while causal prioritisation did not improve consistently.

## Theta and induced-q recovery

The true theta values in binary, continuous, proxy, and noise order were
1.386294, 0.693147, 0, and 0.

| stage | posterior theta means | log-q correlation | log-q RMSE | causal log-q correlation | causal log-q RMSE |
|---|---|---:|---:|---:|---:|
| `q50_h030` | 1.0094, 0.5906, 0.2803, -0.3043 | 0.9078 | 0.4152 | 0.9271 | 0.3667 |
| `q500_h030` | 1.5222, 0.6521, -0.6362, 0.3785 | 0.6640 | 0.7702 | 0.6712 | 0.7932 |

For `q50_h030`, maximum theta R-hat was 1.0304 and minimum bulk ESS was
129.9, giving encouraging functional recovery. For `q500_h030`, maximum R-hat
was 1.2007, minimum bulk ESS 14.6, and minimum tail ESS 45.7. Mixing and
separation of informative, correlated-proxy, and noise effects were therefore
inadequate at that stage; its theta estimates must not be treated as stable.

## Ordinary SBayesR operator boundary at `q500_h050`

The truth identities passed, but ordinary SBayesR failed at iteration zero
before SBayesRV ran. Diagnostics were $y^\mathsf{T}y=7130.54$,
$b^\mathsf{T}X^\mathsf{T}y=9402.56$, and
$b^\mathsf{T}C_{\mathrm{sparse}}b=9072.37$, producing SSE $-2602.22$ with
1,006 active markers. A nonnegative SSE required a quadratic of at least
11674.59, so the missing accumulated quadratic was at least 2602.22.

Thousands of individually weak correlations omitted by hard entrywise CSR
sparsification accumulated in this polygenic state. This is an operator-
approximation boundary, not evidence against gsim, SBayesRV, or its theta
transition. The full-panel Schur residual is not retained as a blocking
diagnostic for hard-thresholded LD when $M>N$.
