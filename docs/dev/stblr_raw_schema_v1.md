# ST-BLR Raw Schema v1

## Purpose

The named raw schema replaces fragile positional native return slots with a
stable internal list contract. Different ST-BLR CSR models currently reuse the
same slot positions for different concepts, which makes formatter changes risky
as new diagnostics and model families are added.

Raw schema v1 is an internal development interface. Public `fit` object field
names remain unchanged.

## Top-Level Namespaces

Version 1 uses these top-level names:

```r
schema
meta
marker
trace
variance
pi
diagnostics
chains
prior
group
annotation
component
selection
```

Unused namespaces may be empty named lists or `NULL` for a backend that does
not use them.

## Canonical Dimensions

Marker-level matrices use marker rows and trait columns:

```text
m x nt
```

Trace matrices use MCMC iteration rows and trait columns:

```text
n_trace x nt
```

For ordinary CSR BayesC, `marker$dm` is `P(d = 1)` and `trace$pis` is the
sampled active-marker probability `pi1` for every saved iteration.

## First Migrated Backend

Phase 1 migrates only ordinary summary-statistic CSR BayesC:

```text
src/st_cpg_omp_csr.cpp
```

This backend returns `schema`, `meta`, `marker`, `trace`, `variance`, `pi`,
`diagnostics`, `chains`, and `selection`. The `prior`, `group`, `annotation`,
and `component` namespaces are present as empty lists.

The R formatter `.format_stblr_raw_v1()` consumes only this named schema and
maps it back to the existing user-facing fit fields such as `bm`, `dm`, `vbs`,
`vgs`, `ves`, `vle`, `vld`, `pi`, `pim`, `pis`, `chains`, `ld_swap`, and
sampled `selection_s` summaries.

## Non-Migrated Backends

The old positional formatter remains active for non-migrated backends,
including scheduled CSR BayesC, CSR BayesR, CSR SBayesRC, prior/group/learned
annotation CSR BayesC, BED backends, and individual-level scheduled backends.

Wrappers should explicitly detect:

```r
raw$schema$class == "stblr_raw"
as.integer(raw$schema$version) == 1L
```

before routing through `.format_stblr_raw_v1()`.

## Planned Migration Order

The intended follow-up order is:

1. Marker-prior CSR BayesC.
2. Group CSR BayesC.
3. Learned annotation CSR BayesC.
4. CSR BayesR.
5. CSR SBayesRC.
6. Scheduled and BED backends after the summary-statistic CSR schemas have
   stabilized.
