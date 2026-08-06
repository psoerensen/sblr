write_aggregate_bayesr_bed <- function(path, dosage) {
  dosage_to_code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    code <- unname(dosage_to_code[as.character(dosage[marker, ])])
    code <- c(code, rep(0L, (-length(code)) %% 4L))
    vapply(seq(1L, length(code), by = 4L), function(index)
      sum(code[index:(index + 3L)] * c(1L, 4L, 16L, 64L)), integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

aggregate_bayesr_fixture <- function() {
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 1, 0),
    c(2, 0, 1, 2, 0, 1, 2, 1),
    c(0, 0, 1, 1, 2, 2, 1, 0),
    c(2, 2, 1, 1, 0, 0, 1, 2))
  path <- tempfile(fileext = ".bed")
  write_aggregate_bayesr_bed(path, dosage)
  ids <- paste0("rs", seq_len(nrow(dosage)))
  list(path = path, ids = ids, y = matrix(
    c(-1, -.5, 0, .5, 1, 1.5, -.75, .25), ncol = 1L,
    dimnames = list(NULL, "T1")), Glist = list(
      n = ncol(dosage), ids = paste0("id", seq_len(ncol(dosage))),
      bedfiles = path, rsids = list(ids), rsidsLD = list(ids),
      chr = list(rep(1L, length(ids))), pos = list(seq_along(ids) * 100),
      af = list(rowMeans(dosage) / 2), maf = list(rowMeans(dosage) / 2)))
}

aggregate_block_fixture <- function() {
  fixture <- aggregate_bayesr_fixture()
  fixture$stats <- make_summary_stats(
    fixture$Glist, fixture$y, chr = 1L, rows = seq_len(nrow(fixture$y)),
    scale = TRUE, nthreads = 1L)
  fixture
}

expect_native_aggregate_oracle <- function(fit, marker_count, component_count) {
  quantities <- fit$convergence_traces$quantities
  values <- fit$convergence_traces$values
  selected <- which(quantities$group == "selected_component")
  components <- values[, , selected, drop = FALSE]
  count <- fit$component_count_trace
  for (draw in seq_len(dim(components)[1L])) for (chain in
      seq_len(dim(components)[2L])) {
    state <- as.integer(components[draw, chain, ])
    oracle <- tabulate(state + 1L, nbins = component_count)
    expect_identical(as.integer(count[draw, chain, ]), oracle)
    expect_identical(as.integer(fit$realized_active_count_trace[draw, chain, ]),
      sum(state > 0L))
    for (stick in seq_len(component_count - 1L)) {
      component_threshold <- stick - 1L
      eligible <- sum(state >= component_threshold)
      continued <- sum(state > component_threshold)
      expect_identical(as.integer(
        fit$stick_eligible_count_trace[draw, chain, stick]), eligible)
      expect_identical(as.integer(
        fit$stick_continue_count_trace[draw, chain, stick]), continued)
      expect_identical(as.integer(
        fit$stick_stop_count_trace[draw, chain, stick]), eligible - continued)
    }
    expect_identical(sum(oracle), marker_count)
  }
}
