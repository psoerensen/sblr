# BLR model contracts

The public model name combines prior kernel and data level:

| Model | Prior kernel | Data level |
|---|---|---|
| `bayesc` | BayesC | individual |
| `sbayesc` | BayesC | summary statistics |
| `bayesr` | BayesR | individual |
| `sbayesr` | BayesR | summary statistics |
| `bayesrc` | annotation-informed BayesR | individual |
| `sbayesrc` | annotation-informed BayesR | summary statistics |

The `s` prefix means summary statistics. It never activates MAF scaling.
`selection_s` is an independent effect-scale policy using
`[2p(1-p)]^(S+1)`; `NULL` means absent and `-1` is an explicitly requested
unit-scale reduction.

Annotation policies are `global`, `fixed_marker`, `group`,
`learned_logistic`, and `annotation_probit_stick`. MT BayesRC uses
annotation-dependent component probabilities and global conditional trait
pattern probabilities. Annotation values do not alter covariance matrices,
component multipliers, operators, or residual likelihoods.
