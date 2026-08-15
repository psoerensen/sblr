blr_with_legacy_execution <- function(callback) {
  original <- sblr:::resolve_blr_spec_from_wrapper
  testthat::with_mocked_bindings(
    callback(),
    resolve_blr_spec_from_wrapper = function(...) {
      args <- list(...)
      args[["execution_contract_version"]] <- 0L
      do.call(original, args)
    },
    .package = "sblr"
  )
}
