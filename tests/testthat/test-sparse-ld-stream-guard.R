test_that("sparseLD_stream_CSR fails closed when both distance filters are disabled", {
  expect_error(
    sparseLD_stream_CSR(
      bed_files = "missing.bed",
      n = 10L,
      cls = list(1L),
      out_prefix = tempfile(),
      pos_bp = NULL,
      max_distance_bp = 0L,
      max_distance_variants = 0L
    ),
    "Both sparse LD distance filters are disabled"
  )

  expect_error(
    sparseLD_stream_CSR(
      bed_files = "missing.bed",
      n = 10L,
      cls = list(1L),
      out_prefix = tempfile(),
      pos_bp = list(integer()),
      max_distance_bp = 1000000L,
      max_distance_variants = 0L
    ),
    "Both sparse LD distance filters are disabled"
  )
})
