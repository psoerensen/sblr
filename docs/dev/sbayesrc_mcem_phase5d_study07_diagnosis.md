# SBayesRC-EM Phase 5D: Study 07 dispersed-start diagnosis

## Scope and frozen provenance

Phase 5D is a diagnostic audit of the frozen Study 07
`informative_annotations`, replicate 1,
`v2_identifiable_qualification` case. It does not alter either EM method,
its priors, damping, genomic updates, convergence tolerances, or inner Monte
Carlo effort. The benchmark specification identity is
`241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`.
The input truth, truth export, learned block checkpoint, human Glist, and
training LD identities are those recorded in the frozen Study 07 report.

Study 07 remains historical evidence in `sblrbench`. The Phase 5D loader is a
read-only reconstruction of its input loader: it verifies all frozen hashes,
uses the already-existing block CSR files, and writes evidence only below
`sblr/results/local/sbayesrc_mcem/phase5D`.

## Responsibility-conditioned objective diagnostic

For fixed Monte Carlo responsibilities, stick `s` has expected eligibility
`e_is`, expected success `u_is`, and continuation probability
`q_is = Phi(A_i alpha_s)`. The recorded annotation surrogate is

```text
Q_annotation(alpha | r) =
  sum_is { u_is log(q_is) + (e_is-u_is) log(1-q_is) }.
```

The Gaussian coefficient prior contributes

```text
log_prior_alpha(alpha) =
  -1/2 sum_sj precision_js (alpha_js - mean_js)^2,
```

up to constants independent of `alpha`. The diagnostic total is
`Q_total = Q_annotation + log_prior_alpha`. For each outer iteration the
implementation records this decomposition at the current, undamped M-step
target, and damped-updated alpha. The undamped target agrees with the existing
minimized M-step objective by sign to numerical tolerance.

Current/target/updated values within one outer iteration use the same fixed
responsibilities and can therefore diagnose the conditional optimization.
Values from different outer iterations use different Monte Carlo
responsibilities. They are not evaluations of one common observed-data
marginal likelihood and must not be labelled `observed_loglik`.

The optional diagnostics perform deterministic arithmetic only. Their default
is off, and responsibility checkpoints are retained only for explicitly
requested outer iterations.

## Continuous SBayesRC-EM results

### Frozen failure reproduction

The four 200-outer reruns reproduce the original iteration-50 states exactly:
the maximum alpha and pre-existing history differences against every frozen
Study 07 checkpoint are both zero. At iteration 50 the reproduced maximum
between-start differences are 0.8365 for alpha, 0.4629 for component priors,
and 0.4118 for active priors; the minimum active-prior correlation is 0.7139.

### Outer 50 to 200

The unchanged outer trajectories show sustained contraction:

| outer | max alpha difference | max component-prior difference | max active-prior difference | min active-prior correlation | expected-active range |
|---:|---:|---:|---:|---:|---:|
| 50 | 0.8365 | 0.4629 | 0.4118 | 0.7139 | 67.49 |
| 100 | 0.7572 | 0.3541 | 0.3541 | 0.8922 | 50.40 |
| 150 | 0.6859 | 0.2961 | 0.2486 | 0.9624 | 31.10 |
| 200 | 0.3611 | 0.2394 | 0.1918 | 0.9876 | 17.49 |

No start met the unchanged convergence criterion by iteration 200. Final beta
means were nevertheless extremely stable (maximum difference 0.00636; minimum
correlation 0.99969). Final genomic PIPs were more variable (maximum difference
0.234; minimum correlation 0.9863), but substantially more stable than the
annotation-prior solution.

The within-iteration M-step target always improves the fixed-responsibility
surrogate over the current alpha. The trajectories continue to move materially
at iteration 50 and still move at iteration 200. Conditional surrogate values
cannot be ranked across starts as observed-data objectives because each start
has its own responsibilities.

### Fixed-state E-step Monte Carlo variation

Five independent current-effort E-steps were run at both baseline-like and
mixed-like frozen states. Across streams, mean pairwise RB responsibility RMSE
was approximately 0.00046. The induced maximum alpha M-step difference was
0.00670 for the baseline state and 0.00715 for the mixed state; expected-active
ranges were 0.418 and 0.480 markers. Thus current Monte Carlo variation exceeds
the strict 0.001 alpha convergence tolerance, but is far smaller than the
remaining 0.361 between-start alpha separation at outer 200.

In the optional baseline-state effort ladder, 4x effort reduced mean pairwise
RB RMSE to 0.000193, maximum alpha variation to 0.00192, and expected-active
range to 0.117. A three-stream 2x result was noisier than 1x and is treated as a
small-replicate fluctuation rather than a monotone estimate. The 4x result is
consistent with a Monte Carlo noise floor decreasing with effort.

### Frozen-responsibility M-step

For one fixed responsibility matrix, optimization from baseline, mixed,
negative, and positive alpha starts agreed to maximum alpha difference
`2.65e-6`, objective range `2.10e-10`, and maximum induced component-prior
difference `1.30e-6`. No continuous M-step start-dependence was found.

### Continuous diagnosis

The evidence supports a mixed mechanism: slow outer contraction explains the
large dispersed-start separation, while current E-step Monte Carlo precision
prevents the strict convergence criterion from settling near a common fixed
point. It does not support a fixed-responsibility optimizer defect or a frozen
Study 07 provenance/invocation mismatch.

## SBayesRC-S-EM model-space audit

The long selection runs also reproduce every iteration-50 frozen checkpoint
exactly (alpha, EB PIP, delta, and pre-existing history differences are zero).
They do not approach a common annotation model by iteration 200. The four delta
MAP states remain `100`, `101`, `111`, and `110`; maximum EB-PIP difference is
1.0. Maximum component-prior difference changes only from 0.701 at iteration
50 to 0.691 at iteration 200, and minimum active-prior correlation changes from
0.528 to 0.684. Final genomic beta remains substantially more stable than the
annotation solution (minimum beta correlation 0.988), whereas final genomic
PIP correlation is 0.832.

All eight shared-delta configurations were enumerated at each retained
responsibility checkpoint using the Phase 5C responsibility-conditioned
Laplace target. On exactly matched responsibilities, production MC3 and exact
enumeration agree closely. Across outer 50/100/150/200, the largest absolute
MC3-versus-enumerated annotation PIP differences are 0.00819, 0.00419, 0.00493,
and 0.00254, respectively. At iteration 200 the enumerated best models are
exactly the four production delta states and have conditional probabilities
0.979, 0.984, 1.000, and 0.999.

Thus model-space exploration is not the source of the Study 07 selection
failure. Different outer starts feed strongly different responsibility
matrices into an otherwise correctly explored eight-model target. Selection
instability is classified as outer responsibility feedback, not MC3 error.

## Production implications

No implementation bug has been established and production defaults remain
unchanged. The primary continuous classification is `MCEM5D-C6`: slow outer
contraction plus an E-step noise floor at the strict convergence tolerance.
Because dispersed starts have not reached essentially the same solution by
outer 200, the formal Phase 5D outcome is `MCEM5D-R2` rather than the stronger
slow-convergence-resolved label. Selection instability comes from outer
responsibility feedback, not model-space exploration.

A separate Phase 5E should evaluate a preregistered outer-run/Monte-Carlo
precision contract (including continuation/checkpoint support and effort
scheduling) against this same frozen case. It must not alter priors or the
statistical target. The diagnostic additions here do not modify the scientific
transition or RNG stream when disabled.
