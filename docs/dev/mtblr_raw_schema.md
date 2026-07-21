# MT-BLR raw schema version 1

`mtblr_raw` is the internal named result boundary for canonical multivariate
BayesC. Its `schema` namespace contains `class = "mtblr_raw"` and
`version = 1L`.

The namespaces are `schema`, `meta`, `marker`, `trace`, `variance`, `pi`,
`model`, `diagnostics`, `data`, and `alignment`. Marker matrices are `m x nt`;
trace matrices are `(nit + nburn) x nt`; covariance matrices are `nt x nt`.
`pi$final` and `pi$mean` have `nmodels` named entries and
`model$patterns` is an `nmodels x nt` binary matrix. Diagnostics contain the
marker, B, G, E, and probability contribution counts. R-side validated marker,
trait, resource, and alignment metadata are attached after native execution.

The schema contains no CSR buffers and is validated by `.validate_mtblr_raw()`
before `.as_mtblr_fit()` formats the public result.
