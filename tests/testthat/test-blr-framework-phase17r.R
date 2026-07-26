test_that("Phase 17R one-chain execution is an exact serial reduction", {
  for (mode in c("diagonal", "full")) {
    for (updates in c(FALSE, TRUE)) {
      case <- phase17o_case(residual_covariance = mode, updates = updates)
      on.exit(phase17o_cleanup(case), add = TRUE)
      serial <- phase17o_call(case)
      chains <- phase17r_call(case)
      expect_identical(phase17r_numerical_raw(chains),
                       phase17r_numerical_raw(serial))
      expect_true(all(chains$marker$bm_sd == 0))
      expect_true(all(chains$marker$dm_sd == 0))
      expect_identical(chains$marker$bm_min, chains$marker$bm)
      expect_identical(chains$marker$bm_max, chains$marker$bm)
      expect_identical(chains$marker$dm_min, chains$marker$dm)
      expect_identical(chains$marker$dm_max, chains$marker$dm)
      expect_null(chains$chains)
      expect_silent(sblr:::.validate_mtblr_raw(chains))
    }
  }
})

test_that("Phase 17R resolves logical chain seeds independently of workers", {
  case <- phase17o_case(updates = FALSE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  serial <- phase17r_call(case, nchains = 4L, ncores = 1L)
  parallel <- phase17r_call(case, nchains = 4L, ncores = 3L)
  expected <- (as.double(case$seed) + 9176 * 0:3) %% 2^32
  expect_identical(serial$diagnostics$mt_bed$chain_seeds, expected)
  expect_identical(parallel$diagnostics$mt_bed$chain_seeds, expected)
  expect_identical(phase17r_without_timing(serial),
                   phase17r_without_timing(parallel))
  expect_identical(parallel$diagnostics$mt_bed$used_workers, 3L)

  explicit <- phase17r_call(case, nchains = 2L, ncores = 1L,
                            chain_seeds = c(-1L, 19L))
  expect_identical(explicit$diagnostics$mt_bed$chain_seeds,
                   c(2^32 - 1, 19))
})

test_that("Phase 17R serial multichain pooling and stability are deterministic", {
  case <- phase17o_case(nt = 3L, updates = TRUE, multiple_sets = TRUE,
                        nonzero = TRUE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  first <- phase17r_call(case, nchains = 2L, ncores = 1L,
                         keep_chains = TRUE)
  second <- phase17r_call(case, nchains = 2L, ncores = 1L,
                          keep_chains = TRUE)
  expect_identical(phase17r_without_timing(first),
                   phase17r_without_timing(second))
  bm <- simplify2array(lapply(first$chains, function(x) x$marker$bm))
  dm <- simplify2array(lapply(first$chains, function(x) x$marker$dm))
  expect_equal(first$marker$bm, unname(apply(bm, 1:2, mean)),
               tolerance = 1e-15)
  expect_equal(first$marker$dm, unname(apply(dm, 1:2, mean)),
               tolerance = 1e-15)
  expect_equal(first$marker$bm_sd, unname(apply(bm, 1:2, sd)),
               tolerance = 1e-15)
  expect_equal(first$marker$dm_sd, unname(apply(dm, 1:2, sd)),
               tolerance = 1e-15)
  expect_identical(first$marker$b, first$chains$chain1$marker$b)
  expect_identical(first$marker$state, first$chains$chain1$marker$state)
  expect_identical(first$variance$vb, first$chains$chain1$variance$vb)
  expect_identical(first$pi$final, first$chains$chain1$pi$final)
  expect_identical(first$trace$vgs,
                   (first$chains$chain1$trace$vgs +
                      first$chains$chain2$trace$vgs) / 2)
})

test_that("Phase 17R retained chains are compact and extended raw formats", {
  case <- phase17o_case(updates = TRUE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  raw <- phase17r_call(case, nchains = 2L, ncores = 2L,
                       keep_chains = TRUE)
  expect_identical(names(raw$chains), c("chain1", "chain2"))
  expect_identical(names(raw$chains[[1]]),
                   c("chain", "seed", "marker", "trace", "variance", "pi",
                     "diagnostics"))
  expect_identical(names(raw$chains[[1]]$marker),
                   c("bm", "dm", "b", "state"))
  forbidden <- c("r", "wy", "order", "models", "sets", "phenotype",
                 "sample_residual", "genetic_values", "packed_genotype",
                 "marker_maps", "bed_files")
  expect_false(any(forbidden %in% unlist(lapply(raw$chains, names))))
  expect_silent(sblr:::.validate_mtblr_raw(raw))

  raw$model$names <- paste0("M", seq_len(raw$meta$nmodels))
  ids <- paste0("rs", seq_len(raw$meta$m))
  traits <- colnames(case$Y)
  fit <- sblr:::.as_mtblr_fit(
    raw, ids, traits, data.frame(marker_id = ids),
    data.frame(trait_id = traits), list(), list())
  expect_s3_class(fit, "mtblr_fit")
  expect_identical(fit$nchains, 2L)
  expect_identical(names(fit$chains), c("chain1", "chain2"))
  expect_identical(dimnames(fit$bm_sd), list(ids, traits))
  expect_identical(dimnames(fit$chains$chain1$marker$bm), list(ids, traits))
})

test_that("Phase 17R reports all worker failures in deterministic chain order", {
  case <- phase17o_case(updates = FALSE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  case$state[[1]][1] <- 2L
  serial <- tryCatch(phase17r_call(case, 3L, 1L), error = conditionMessage)
  parallel <- tryCatch(phase17r_call(case, 3L, 2L), error = conditionMessage)
  expect_identical(serial, parallel)
  expect_match(serial, "mtblr_bed_chains_internal failed", fixed = TRUE)
  expect_lt(regexpr("chain 1:", serial, fixed = TRUE),
            regexpr("chain 2:", serial, fixed = TRUE))
  expect_lt(regexpr("chain 2:", serial, fixed = TRUE),
            regexpr("chain 3:", serial, fixed = TRUE))
})

test_that("Phase 17R RNG state is fit-local and fresh-process reproducible", {
  skip_if_not_installed("callr")
  case <- phase17o_case(residual_covariance = "full", updates = TRUE)
  on.exit(phase17o_cleanup(case), add = TRUE)
  args <- phase17r_args(case, nchains = 2L, ncores = 1L,
                        keep_chains = TRUE)
  expected <- do.call(sblr:::mtblr_bed_chains_internal, args)
  dense <- phase17o_dense_args(phase17o_case())
  invisible(do.call(sblr:::mtblr, dense))
  observed <- do.call(sblr:::mtblr_bed_chains_internal, args)
  expect_identical(phase17r_without_timing(observed),
                   phase17r_without_timing(expected))

  scalar <- phase17p_case()
  on.exit(phase17p_cleanup(scalar), add = TRUE)
  invisible(stblr_bed(
    scalar$Y[, 1], scalar$Glist, rows = scalar$rows,
    nit = 1L, nburn = 0L, nthin = 1L, ncores = 1L,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE))
  after_scalar <- do.call(sblr:::mtblr_bed_chains_internal, args)
  expect_identical(phase17r_without_timing(after_scalar),
                   phase17r_without_timing(expected))

  root <- if (exists("blr_has_source_tree", mode = "function") &&
              blr_has_source_tree()) blr_source_root else NULL
  fresh <- callr::r(function(args, root) {
    if (!is.null(root)) pkgload::load_all(root, compile = FALSE, quiet = TRUE)
    else library(sblr)
    do.call(getFromNamespace("mtblr_bed_chains_internal", "sblr"), args)
  }, list(args = args, root = root))
  expect_identical(phase17r_without_timing(fresh),
                   phase17r_without_timing(expected))
})

test_that("Phase 17R preserves the public route and source architecture", {
  expect_false(any(c("nchains", "ncores", "chain_seeds", "keep_chains") %in%
                     names(formals(mtblr_bed))))
  root <- blr_repo_path()
  skip_if(is.null(root), "source checkout unavailable")
  mt <- readLines(file.path(root, "src", "mtblr.cpp"), warn = FALSE)
  execution <- readLines(file.path(
    root, "src", "blr_mt_bed_chains_execution_impl.h"), warn = FALSE)
  public <- readLines(file.path(root, "R", "mtblr-bed.R"), warn = FALSE)
  expect_equal(sum(grepl("mtblr_bed_chains_internal(", mt, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("#pragma omp parallel for schedule(static)",
                         execution, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("run_mt_bed_bayesc_core(", execution,
                         fixed = TRUE)), 1L)
  expect_equal(sum(grepl("mtblr_bed_internal(", public, fixed = TRUE)), 1L)
  expect_false(any(grepl("mtblr_bed_chains_internal", public, fixed = TRUE)))
})
