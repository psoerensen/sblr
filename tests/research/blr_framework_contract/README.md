# Unified BLR Phase 0 contract fixtures

These base-R fixtures make the approved Phase 0 contracts independently
executable before production resolvers, native boundaries, or schema builders
exist. They are not package validators and are not called by a production
fitter.

The fixtures freeze:

- `blr_resolved_spec` version 1 envelope names;
- analysis/execution combinations and logical task-seed shapes;
- seed-contract version 1, including fixed reference vectors;
- retention-contract version 1 and its exact post-burn indices;
- immutable operator resources referenced by one or more likelihood providers;
- common-sample joint providers versus independent singleton summary providers;
- `blr_raw` schema version 2 names, fixed axes, and present-but-`NULL` fields;
- fixed-full and sampled-full residual-covariance policies.

Run from the repository root:

```powershell
Rscript tests/research/blr_framework_contract/test_contract_fixtures.R
```

The fixtures deliberately do not emulate current raw schema version 1. Current
behavior, Phase 1 implementation status, and migration rules are maintained in
the developer contracts. A passing fixture test remains an independent
contract oracle rather than evidence derived from production constructors.
Fixture compatibility identifiers never bypass probability, covariance, axis,
or provenance validation.
