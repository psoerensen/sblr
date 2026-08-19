# SBayesRV research

SBayesRV is the canonical research name for the family in which annotations or
supplied marker signals modify active-effect variance within non-null BayesR
components. It is the research continuation of the maintained implementation
currently labelled BayesR-LV/SBayesR-LV. SBayesRC is different: it changes
component-membership probabilities, whereas SBayesRV changes the relative
active-effect multiplier $q_j$ while keeping global BayesR component
probabilities separate.

This directory is research-only. The future identifier `"sbayesrv"` is
reserved, but no public method rename, alias, package export, or second
production sampler has occurred.

Maintained files:

- `theory.md`: frozen version-1 mathematical target, implementation crosswalk,
  reductions, identification, and likelihood boundary.
- `prototype.R`: independent base-R conditional and collapsed reference
  calculations plus a small elliptical-slice transition.
- `qualification.R`: compact deterministic Gate 1 fixtures and checks.
- `analysis.R`: the single Gate 1 entry point; any generated tables go only to
  ignored `output/`.
- `prior_study12/`: preserved historical input. Its SBayesR-LV naming, designs,
  results, and conclusions remain unchanged for provenance and are not promoted
  benchmark evidence.

The research gate sequence is:

```text
theory -> independent prototype -> deterministic qualification
       -> gsim evaluation -> production decision -> benchmark
```

Passing Gate 1 establishes internal mathematical agreement only. It does not
show improved prediction, prioritization, effect recovery, theta recovery, or
operating characteristics.
