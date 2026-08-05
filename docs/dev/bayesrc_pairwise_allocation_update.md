# Pairwise BayesRC allocation design review

## Decision

An exact two-marker component/effect conditional was derived and validated,
but the proposed fixed-pair transition did not improve the demonstrated Study
06 occupancy bottleneck. It is therefore **not a production sampler option**.
The package retains only a development reference kernel and deterministic tests.
Scalar BED, CSR, and retained block-eigen trajectories and public interfaces
remain unchanged.

The next method-development step is a separately reviewed larger collapsed
allocation move. No Study 06 benchmark or qualification contract is changed by
this work.

## Existing LD swap

The current SBayesRC LD swap selects an active marker and a null high-LD
friend, transfers the effect magnitude and component label, and applies an MH
correction for the state-dependent proposal, component prior, effect prior, and
residual-likelihood change. Consequently it preserves total active-marker count
and every component-specific occupancy exactly. It can improve localization
mixing but cannot directly reconcile chains in different occupancy regimes.

## Exact pair conditional

For a pair (S=\{i,j\}), let (G_S) be the fitted operator submatrix and let
(s_S) be the corrected score after both current effects are removed. For a
candidate component pair (c=(c_i,c_j)), restrict to its active subset (A).
The diagonal prior covariance (V_A) includes the current common marker
variance, component multiplier, and any marker-specific effect scale. With
residual variance (v_e),

\[
P_A=G_A/v_e+V_A^{-1},\qquad h_A=s_A/v_e,
\]

\[
C_A=P_A^{-1},\qquad m_A=C_Ah_A.
\]

After integrating the active effects, the component-pair log weight is

\[
\log w(c)=\log\pi_i(c_i)+\log\pi_j(c_j)
-\tfrac12\{\log|V_A|+\log|P_A|\}
+\tfrac12h_A^\mathsf{T}P_A^{-1}h_A.
\]

The null-null determinant and quadratic contributions are zero. Conditional on
the sampled component pair, active effects follow (N(m_A,C_A)). Thus, for a
pair chosen independently of current components and effects, this is an
ordinary random-scan Gibbs transition and preserves the existing posterior.
There is no acceptance probability.

The existing one-marker conditional is recovered when the second marker is
fixed null. At zero cross-product, the integrated component probabilities
factor into the two one-marker conditionals.

## Deterministic validation

An independent dense-Gaussian reference verifies all component-pair weights,
conditional means, and covariances for:

- zero, moderate positive, strong positive, and negative cross-products;
- unequal operator diagonals;
- unequal annotation component probabilities;
- unequal marker-specific effect scales;
- marker-exchange symmetry and null-dominant priors.

The reference kernel rejects a pair submatrix that is not positive
semidefinite. Tests also verify zero-cross-product factorization. These checks
establish the mathematics of the candidate transition; they do not establish
practical mixing value.

## Study 06 development evidence

All runs used uninformative replicate 1 BED identities and the original four
chain seeds (`801121`, `801222`, `801323`, `801424`). They were development
diagnostics, not qualification runs, and the sibling evidence files were hash
checked unchanged.

The proposed fixed graph examined the next 20 marker indices using the actual
selected standardized-genotype cross-product:

| Minimum r-squared | Stored edges | Stage A result |
|---:|---:|:---|
| 0.50 | 0 | Pair transition inactive; 1,000-iteration result exactly reproduced the committed baseline diagnosis |
| 0.10 | 0 | Pair transition inactive |
| 0.01 | 140 | 191 completed pair Gibbs moves across four 1,000-iteration chains; only one changed component state and occupancy |

For the r-squared 0.01 ablation, maximum R-hat was 2.254, minimum bulk ESS was
5.14, and the final component-1 occupancies were 12, 3, 0, and 15. The committed
single-marker Stage A reference had maximum R-hat about 2.122, minimum bulk ESS
about 5.22, and component-1 occupancies 51, 14, 11, and 9. Implied prior
summaries still failed to converge. The candidate did not improve effective
exploration and sometimes produced broader late-stick excursions.

The diagnostic r-squared change was not used to change a default. It was an
explicit post-failure ablation to determine whether the empty high-LD graph or
the pair transition itself limited performance.

## Why the move was not retained

The useful Study 06 graph was tiny relative to 37,991 markers. At the tested
frequency it almost always redrew the same state. Lower-correlation arbitrary
pairs would approach two independent single-marker updates and provide little
reason to expect a new occupancy transition. Increasing frequency or selecting
pairs from the current active set would either add substantial cost or require
an exact state-dependent proposal calculation; neither had demonstrated value.

The evidence therefore supports the larger-block outcome: occupancy movement
requires a transition that can jointly reorganize more than two competing
marker states, likely with collapsed effect integration and an explicitly
validated fixed or corrected block-selection law. That design is outside this
task.

## Scope and limitations

The retained C++ code is a development reference conditional exposed only to
package tests. It is not called by BED, exact CSR, hard-sparse CSR, retained
block eigen, or reconstructed block eigen. Multitrait BayesRC is unchanged.
The result does not alter the proper probit-stick prior, LD-operator roles,
Study 06 thresholds, or any benchmark conclusion.
