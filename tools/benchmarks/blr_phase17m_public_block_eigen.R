root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pkgload::load_all(root, compile = FALSE, quiet = TRUE)
source(file.path(root, "tests/testthat/helper-mtblr-block-eigen-fixtures.R"),
       local = TRUE)

measure <- function(label, case, sharing, filter, blocks) {
  nt <- length(case$stats$wy)
  shared <- sharing == "shared"
  filter_values <- rep(filter, length.out = nt)
  block_values <- if (is.list(blocks)) blocks else rep(list(blocks), nt)
  native_blocks <- lapply(block_values, function(x) as.integer(x - 1L))
  internal <- phase17l_case(
    nt = nt, filters = filter_values, shared = shared,
    blocks = native_blocks)
  build_elapsed <- system.time({
    inspections <- if (shared) {
      list(phase17l_inspect(internal$block$operator_descriptors[[1L]],
                            internal$block$wy))
    } else {
      lapply(seq_len(nt), function(t)
        phase17l_inspect(internal$block$operator_descriptors[[t]],
                         list(internal$block$wy[[t]])))
    }
  })[["elapsed"]]
  internal_elapsed <- system.time(
    do.call(sblr:::mtblr_block_eigen_internal, internal$block)
  )[["elapsed"]]
  gc()
  elapsed <- system.time({
    fit <- phase17m_call(
      case, operator_sharing = sharing, eigen_filter = filter,
      block_start = blocks, nit = 20L, nburn = 5L)
  })[["elapsed"]]
  data.frame(
    case = label, sharing = sharing,
    owner_count = fit$block_diagnostics$owner_count,
    filter = paste(filter, collapse = ","),
    total_public_seconds = elapsed,
    fit_bytes = as.numeric(object.size(fit)),
    operator_bytes = as.numeric(object.size(inspections)),
    public_preparation_seconds = max(0, elapsed - internal_elapsed),
    native_build_seconds = build_elapsed,
    mcmc_seconds = max(0, internal_elapsed - build_elapsed),
    internal_call_seconds = internal_elapsed,
    stringsAsFactors = FALSE)
}

small <- phase17m_public_case()
results <- rbind(
  measure("small_shared_hard", small, "shared", "hard_truncate",
          c(1L, 3L)),
  measure("small_shared_ridge", small, "shared", "ridge_fixed",
          c(1L, 3L)),
  measure("small_trait_specific", small, "trait_specific",
          c("hard_truncate", "ridge_lw"),
          list(c(1L, 3L), c(1L, 2L, 4L))))
print(results, row.names = FALSE)
cat("Phase 17M benchmark is a regression signal; preparation and MCMC",
    "components are residual timing estimates around the separable internal",
    "build measurement.\n")
