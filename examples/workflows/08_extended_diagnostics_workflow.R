# Phase 21 extended convergence workflow. These small iteration counts show the
# interface and are not convergence recommendations. The objects `stats`,
# `ld_prefix`, `ld_metadata`, `y`, `Glist`, and `annotations` are prepared by
# the maintained unified-operator and annotation workflows.

fit_core <- mtblr_csr(
  stats, ld_prefix = ld_prefix, ld_metadata = ld_metadata,
  method = "sbayesr", mixture_var = c(0, .01, .1, 1),
  nchains = 2, ncores = 2, nit = 100, nburn = 50,
  convergence = "core", convergence_control = list(warn = FALSE))

fit_extended <- mtblr_bed(
  y, Glist, method = "bayesr", mixture_var = c(0, .01, .1, 1),
  nchains = 2, ncores = 2, nit = 100, nburn = 50,
  convergence = "extended",
  convergence_control = list(
    warn = FALSE,
    extended_groups = c("covariance", "probability"),
    selected_markers = c(2L, 1L),
    selected_marker_quantities = c("b", "d", "component"),
    keep_traces = TRUE,
    max_trace_gb = 1,
    allow_large_traces = FALSE))

fit_annotation <- mtblr_csr(
  stats, ld_prefix = ld_prefix, ld_metadata = ld_metadata,
  method = "sbayesrc", annotations = annotations,
  mixture_var = c(0, .01, .1, 1), nchains = 2, ncores = 2,
  nit = 100, nburn = 50, convergence = "extended",
  convergence_control = list(
    warn = FALSE, extended_groups = c("probability", "annotations")))

fit_extended$convergence$group_overview
fit_extended$convergence$summary
fit_extended$convergence_traces$quantities
fit_annotation$convergence$summary

# Diagnostics assess chain mixing. Passing thresholds does not prove model
# correctness or absolute convergence, and selected-marker rows are not
# automatic fine-mapping evidence.
