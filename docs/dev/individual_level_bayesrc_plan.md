# Individual-Level BayesRC Plan

This developer note tracks the staged internal work toward a possible
individual-level BayesRC backend. It does not define or expose a public model.

## Sequence 1 Status

- likelihood-independent probit stick-breaking annotation-prior utilities were
  extracted into `src/st_bayesrc_annotation_prior.h`
- the existing CSR SBayesRC sampler remains the only consumer
- no individual-level BED BayesRC backend was added
- there are no public R API changes and no new statistical model
- the existing SBayesRC component ordering, annotation updates, RNG draw order,
  and formatted output are preserved

## Sequence 2 Status

- added the internal native full-sweep BED BayesRC backend
  `stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc`
- the backend uses the packed-BED BayesR likelihood with the shared probit
  stick-breaking annotation prior
- the public `stblr_bed()` API is not extended in this sequence
- every MCMC iteration visits every selected marker; adaptive null scheduling
  is intentionally absent
- annotation rows are supplied in the final selected BED marker order defined
  by the concatenated `cls` marker map; the low-level caller is responsible for
  applying that same index to `A`
- `selection_s`, LD-swap, annotation-dependent effect variances, and
  multivariate covariance sampling remain deferred

