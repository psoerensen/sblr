# SBayesRV research checkpoint

SBayesRV is the research name for the variance-modulated BayesR family:

$$
\beta_j\mid c_j=k>0 \sim
N\!\left(0,v_b\gamma_k q_j\right),
\qquad q_j=\exp(X_j\theta).
$$

The null component remains exactly zero and the global BayesR component
probabilities remain annotation independent. This differs from SBayesRC,
where annotations change component-membership probabilities. The maintained
package implementation is still publicly named BayesR-LV/SBayesR-LV; this
research checkpoint does not rename or extend the public interface.

## Qualification status

Gate 1 passed. The self-contained theory, independent R oracle, conditional
and collapsed gradients, reductions, production/native crosswalk, and the
deterministic qualification suite agree. Run that bounded qualification from
the repository root with:

```powershell
Rscript research/sbayesrv/analysis.R
```

The earlier small low-LD pilot found `NO_ORACLE_Q_GAIN_IN_SMALL_PILOT`. It had
only 40 causal markers informing four correlated annotation coefficients and
did not provide a regime in which oracle $q$ materially improved inference.
That underpowered result is not a general negative conclusion about SBayesRV.

## Study 02-coordinate checkpoint

The retained one-replicate experiment reused the exact validated Study 02
coordinate: qgdata `human_independent`, 37,991 post-QC markers, the fixed
3,500/1,500 split, training-only allele frequencies, and the existing CSR LD
operator. No LD was rebuilt.

At both the 50-causal and 500-causal $h^2=0.30$ stages, SBayesRV modestly
improved held-out genetic-value correlation and NMSE and reduced effect RMSE
relative to ordinary SBayesR. Causal prioritisation did not improve
consistently. Theta and induced-$q$ recovery were encouraging with 50 causal
markers; with 500, theta mixing and separation of informative, proxy, and
noise coefficients were inadequate. These are descriptive observations from
one replicate, not benchmark evidence. The checkpoint status is
`SBayesRV_PREDICTION_SIGNAL_OBSERVED_AT_STUDY02_COORDINATE`.

The 500-causal, $h^2=0.50$ stage stopped in ordinary SBayesR before SBayesRV
ran. Hard entrywise CSR sparsification underestimated the accumulated effect-
state quadratic $b^\mathsf{T}C_{\mathrm{sparse}}b$ across 1,006 active
markers. This is recorded separately as
`ORDINARY_SBayesR_SPARSE_LD_BOUNDARY_AT_Q500_H050`; it is an operator-
approximation boundary, not evidence against gsim, SBayesRV, or the SBayesRV
transition. The previously proposed full-panel Schur-residual condition is not
an authoritative gate for hard-thresholded LD when $M>N$.

See [`study02_coordinate/findings.md`](study02_coordinate/findings.md) for the
compact numerical checkpoint. Further operator work belongs in
[`../operator_approximation/`](../operator_approximation/), outside SBayesRV.

## Active files

- `theory.md`: version-1 model, likelihood boundary, reductions, and native
  implementation crosswalk.
- `prototype.R`: independent conditional and collapsed reference equations.
- `qualification.R`: compact deterministic qualification.
- `analysis.R`: deterministic qualification entry point only.
- `study02_coordinate/`: the frozen design, rerunnable research code, and
  tracked findings for the one-replicate coordinate.

## Next research boundary

Any next SBayesRV experiment should be prespecified and replicated, use an LD
operator qualified for its intended polygenicity, and assess prediction,
effect recovery, prioritisation, theta/$q$ recovery, and convergence without
changing the production interface. Operator approximation must be qualified
as a separate problem before scientific interpretation.
