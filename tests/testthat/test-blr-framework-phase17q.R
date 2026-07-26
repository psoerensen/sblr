test_that("Phase 17Q task and seed contracts are deterministic", {
  expect_equal(phase17q_tasks(4, 3)$chain, 0:3)
  expect_equal(nrow(phase17q_tasks(4, 30)), 4L)
  expect_equal(phase17q_tasks(4, 2)$result_slot, 0:3)
  expect_equal(phase17q_seeds(11, 4), c(11, 9187, 18363, 27539))
  expect_equal(phase17q_seeds(-1, 2), c(2^32 - 1, 9175))
  expect_equal(phase17q_seeds(2^32 - 10, 2), c(2^32 - 10, 9166))
  expect_equal(phase17q_seeds(1, 3, c(9, -1, 7)), c(9, 2^32 - 1, 7))
})

test_that("Phase 17Q aggregation pools samples and preserves primary state", {
  a <- phase17q_chain(1L, 2L, matrix(c(1, 2), 1), matrix(c(.1, .2), 1),
                      matrix(c(1L, 0L), 1), 0)
  b <- phase17q_chain(2L, 6L, matrix(c(5, 6), 1), matrix(c(.5, .6), 1),
                      matrix(c(0L, 1L), 1), 6)
  out <- phase17q_aggregate(list(b, a), keep_chains = TRUE)
  expect_equal(out$bm, matrix(c(4, 5), 1))
  expect_equal(out$dm, matrix(c(.4, .5), 1))
  expect_equal(out$b, a$b)
  expect_identical(out$state, a$state)
  expect_false(any(out$state == .5))
  expect_equal(out$vbs, (a$vbs + b$vbs) / 2)
  expect_equal(out$covb, (a$covb_acc + b$covb_acc) / 8)
  expect_equal(names(out$chains), c("chain1", "chain2"))
  expect_false(any(c("r", "phenotype", "packed", "wy", "marker_order") %in%
                     names(out$chains[[1]])))
  expect_equal(out$bm_stability$sd, apply(simplify2array(list(
    a$bm_acc / 2, b$bm_acc / 6)), 1:2, sd))
})

test_that("Phase 17Q single-chain stability and failures are explicit", {
  a <- phase17q_chain(1L, 3L, matrix(c(2, 4), 1), matrix(c(.2, .4), 1),
                      matrix(c(1L, 1L), 1))
  out <- phase17q_aggregate(list(a))
  expect_equal(out$bm_stability$sd, matrix(0, 1, 2))
  expect_equal(out$bm_stability$min, out$bm)
  expect_equal(out$bm_stability$max, out$bm)
  expect_null(out$chains)
  first <- phase17q_chain(1L, 1L, matrix(1), matrix(1), matrix(1),
                          failed = TRUE, error = "first")
  second <- phase17q_chain(2L, 1L, matrix(1), matrix(1), matrix(1),
                           failed = TRUE, error = "second")
  expect_error(phase17q_aggregate(list(second, first)),
               "chain 1: first; chain 2: second", fixed = TRUE)
})

test_that("Phase 17Q memory scales shared and private objects separately", {
  one <- phase17q_memory(100, 200, 3, 8, 20, 1, 1)
  many <- phase17q_memory(100, 200, 3, 8, 20, 8, 4, TRUE)
  expect_equal(one$worker_count, 1L)
  expect_equal(many$worker_count, 4L)
  expect_equal(one$shared_bytes, many$shared_bytes)
  expect_equal(many$estimated_concurrent_bytes,
               many$shared_bytes + 4 * many$private_state_bytes_per_chain)
  expect_equal(many$estimated_retained_output_bytes,
               8 * many$retained_chain_bytes_per_chain)
  expect_match(many$estimate_kind, "not measured RSS")
})

test_that("Phase 17Q future controls and current public boundary are explicit", {
  controls <- phase17q_future_controls()
  expect_identical(controls$nchains, 1L)
  expect_identical(controls$ncores, 1L)
  expect_null(controls$chain_seeds)
  expect_false(controls$keep_chains)
  expect_identical(controls$openmp_unavailable, "warn once and run serial")
  expect_identical(controls$raw_schema, "mtblr_raw version 1")
  expect_identical(controls$primary_chain, 1L)
  expect_match(phase17q_future_internal_signature(), "int nchains")
  current <- names(formals(mtblr_bed))
  expect_false(any(c("nchains", "ncores", "chain_seeds", "keep_chains") %in% current))
})

test_that("Phase 17Q source architecture is audit-only", {
  root <- blr_repo_path()
  skip_if(is.null(root), "source checkout unavailable")
  public <- readLines(file.path(root, "R", "mtblr-bed.R"), warn = FALSE)
  native <- readLines(file.path(root, "src", "blr_mt_bed_core_impl.h"), warn = FALSE)
  exports <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  contract <- readLines(file.path(root, "docs", "dev",
                                  "blr_mt_bed_multichain_contract.md"), warn = FALSE)
  expect_false(any(grepl("nchains|ncores|chain_seeds|keep_chains", public)))
  expect_false(any(grepl("#pragma omp", native, fixed = TRUE)))
  expect_false(any(grepl("mtblr_bed_chains_internal", exports, fixed = TRUE)))
  required <- c("task_count = nchains", "schedule(static)",
                "primary_chain = 1", "warn once and run serial",
                "mtblr_raw version 1", "9176", "no partial")
  expect_true(all(vapply(required, function(x) any(grepl(x, contract, fixed = TRUE)),
                         logical(1))))
})
