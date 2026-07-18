# Phase 14A packed-BED BayesRC references

The `.rds` files in this directory are permanent normalized raw-and-formatted
references captured from commit `c446d1c` with R 4.4.1/Rtools44. They use the
two-marker Phase 11A SNP-major BED fixture, four ordered components
`0, 0.01, 0.1, 1`, null component zero, seed 141 or 143, ten total iterations,
and the current marker-by-annotation probit-stick model. Fixture regeneration is
manual maintenance tooling only; ordinary tests never overwrite these files.
