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

