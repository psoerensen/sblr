# Frozen Study 02-coordinate design

The experiment uses `human_independent`, the exact historical 37,991-marker
QC panel, split seed 3101 (3,500 training and 1,500 validation individuals),
training-subset allele frequencies, and the validated CSR operator with zero
base-pair window, 1,000-variant window, $r^2=0.001$, block size 1,024, and one
construction thread. Phenotype-specific score statistics are recomputed for
each stage; LD is resolved from exact cache identity
`8e990e54be49212386cb47b4f14023048d8ae142f211a9c3b549e031abb15565`
and is never rebuilt by this analysis.

The annotation seed is 20260821. The centred, no-intercept design contains a
binary informative annotation, a standardised continuous informative
annotation, a correlated proxy (target correlation 0.70), and noise. The
variance truth is

$$
\theta=(\log 4,\log 2,0,0),\qquad q_j=\exp(X_j\theta),
$$

with geometric mean one. Annotations affect active-effect variance only;
component membership is annotation independent.

| stage | causal markers | $h^2$ | architecture seed | training residual seed | validation residual seed |
|---|---:|---:|---:|---:|---:|
| `q50_h030` | 50 | 0.30 | 20260901 | 20260902 | 20260903 |
| `q500_h030` | 500 | 0.30 | 20260911 | 20260912 | 20260913 |
| `q500_h050` | 500 | 0.50 | 20260921 | 20260922 | 20260923 |

The 50-causal fitting prior is the exact Study 02 prior. The 500-causal prior
has active probability $500/37991$ split 0.60/0.30/0.10 across positive
components. Both methods use `nit=2000`, `nburn=250`, `nthin=1`, four chains,
four workers, fit seed 30104, chain seeds
130104/230104/330104/430104, and no LD-swap update.

Validation phenotypes are inaccessible to fitting and are evaluated only
after both fits for a stage pass canonical extraction and diagnostic checks.
The $h^2=0.50$ stage stopped in ordinary SBayesR before validation evaluation
or an SBayesRV fit. No tuning or rerun followed that failure.
