pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("tests/testthat/helper-mtblr-bed-contract.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-public.R", local = TRUE)

rows <- list()
index <- 0L
for (mode in c("diagonal", "full"))
  for (convergence in c("none", "auto", "core"))
    for (nchains in c(1L, 2L, 4L))
      for (ncores in c(1L, 2L, 4L))
        for (keep_chains in c(FALSE, TRUE))
          for (keep_traces in c(FALSE, TRUE)) {
            if (convergence == "none" && keep_traces) next
            case <- phase17p_case(nt = 2L)
            control <- list(warn = FALSE, keep_traces = keep_traces)
            args <- phase17p_public_args(
              case, mode, updates = TRUE, center = FALSE,
              nchains = nchains, ncores = ncores,
              keep_chains = keep_chains, convergence = convergence,
              convergence_control = control,
              nit = 8L, nburn = 2L, memory_warning_gb = Inf)
            elapsed <- system.time(
              fit <- suppressWarnings(do.call(mtblr_bed, args)))["elapsed"]
            index <- index + 1L
            rows[[index]] <- data.frame(
              mode = mode, convergence = convergence,
              nchains = nchains, requested_cores = ncores,
              used_workers = fit$chain_diagnostics$used_workers,
              keep_chains = keep_chains, keep_traces = keep_traces,
              native_route = fit$input$convergence_trace_route,
              native_dispatch_seconds = fit$chain_diagnostics$dispatch_seconds,
              total_public_seconds = unname(elapsed),
              ordinary_memory_gib =
                fit$memory_estimate$estimated_total_gib -
                fit$memory_estimate$convergence_estimated_total_gib,
              convergence_memory_gib =
                fit$memory_estimate$convergence_estimated_total_gib,
              total_memory_gib = fit$memory_estimate$estimated_total_gib,
              fit_bytes = as.numeric(object.size(fit)),
              convergence_bytes = as.numeric(object.size(fit$convergence)),
              trace_bytes = as.numeric(object.size(fit$convergence_traces)),
              warning_emitted = fit$input$convergence_warning_emitted)
            unlink(case$fixture$paths)
          }
print(do.call(rbind, rows), row.names = FALSE)
cat("Phase 17V public regression signals only; no pure-adapter or linear-speedup claim.\n")
