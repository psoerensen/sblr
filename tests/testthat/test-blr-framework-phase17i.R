phase17i_src <- function(file) paste(readLines(
  blr_repo_path(file), warn = FALSE), collapse = "\n")

phase17i_float <- function(x) {
  con <- rawConnection(raw(), "w+b"); on.exit(close(con))
  writeBin(as.numeric(x), con, size = 4L, endian = .Platform$endian)
  seek(con, 0L); readBin(con, "numeric", n = length(x), size = 4L,
                         endian = .Platform$endian)
}

phase17i_write_operator <- function(prefix, correlation, diagonal) {
  m <- nrow(correlation)
  edges <- which(upper.tri(correlation) & correlation != 0, arr.ind = TRUE)
  edges <- edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
  row_counts <- tabulate(edges[, 1L], nbins = m)
  row_ptr <- c(0, cumsum(row_counts))
  # These tiny fixtures fit in 32 bits; write each uint64 as low/high words.
  writeBin(as.integer(c(rbind(row_ptr, rep(0, length(row_ptr))))),
           paste0(prefix, ".row_ptr.u64.bin"), size = 4L,
           endian = .Platform$endian)
  writeBin(as.integer(edges[, 2L] - 1L),
           paste0(prefix, ".col_idx.u32.0based.bin"),
           size = 4L, endian = .Platform$endian)
  writeBin(as.numeric(correlation[edges]), paste0(prefix, ".values.f32.bin"),
           size = 4L, endian = .Platform$endian)
  writeLines(c(sprintf("n_variants=%d", m), sprintf("nnz=%d", nrow(edges))),
             paste0(prefix, ".meta.txt"))

  indices <- vector("list", m); values <- vector("list", m)
  for (k in seq_len(nrow(edges))) {
    i <- edges[k, 1L]; j <- edges[k, 2L]
    x <- phase17i_float(correlation[i, j] * sqrt(diagonal[i] * diagonal[j]))
    indices[[i]] <- c(indices[[i]], j - 1L); values[[i]] <- c(values[[i]], x)
    indices[[j]] <- c(indices[[j]], i - 1L); values[[j]] <- c(values[[j]], x)
  }
  list(indices = indices, values = values)
}

phase17i_case <- function(nt = 2L, trait_specific = FALSE,
                          independent = FALSE, nonzero_b = FALSE,
                          multiple_sets = FALSE, updates = TRUE) {
  m <- 4L
  base <- matrix(c(1, .18, 0, .05, .18, 1, .12, 0,
                   0, .12, 1, .22, .05, 0, .22, 1), m, m, byrow = TRUE)
  matrices <- rep(list(base), nt)
  if (trait_specific) {
    for (t in seq_len(nt)) {
      matrices[[t]][upper.tri(matrices[[t]])] <-
        matrices[[t]][upper.tri(matrices[[t]])] * (1 - .12 * (t - 1L))
      matrices[[t]][lower.tri(matrices[[t]])] <-
        t(matrices[[t]])[lower.tri(matrices[[t]])]
    }
  }
  if (independent && nt >= 2L) {
    matrices[[2L]][1, 2] <- matrices[[2L]][2, 1] <- 0
    matrices[[2L]][1, 3] <- matrices[[2L]][3, 1] <- .09
  }
  ww <- lapply(seq_len(nt), function(t) rep(1 + if (trait_specific) .1*(t-1) else 0, m))
  prefix <- file.path(tempdir(), paste0("phase17i_", Sys.getpid(), "_", sample.int(1e8, 1)))
  prefixes <- paste0(prefix, "_", seq_len(if (trait_specific || independent) nt else 1L))
  operators <- vector("list", nt)
  if (length(prefixes) == 1L) {
    built <- phase17i_write_operator(prefixes, matrices[[1]], ww[[1]])
    operators <- rep(list(built), nt)
  } else {
    for (t in seq_len(nt)) operators[[t]] <-
      phase17i_write_operator(prefixes[t], matrices[[t]], ww[[t]])
  }
  union_indices <- lapply(seq_len(m), function(i)
    unique(unlist(lapply(operators, function(op) op$indices[[i]]))))
  XXindices <- lapply(seq_len(m), function(i) c(i - 1L, union_indices[[i]]))
  XXvalues <- lapply(seq_len(nt), function(t) lapply(seq_len(m), function(i) {
    positions <- match(union_indices[[i]], operators[[t]]$indices[[i]])
    c(ww[[t]][i], ifelse(is.na(positions), 0, operators[[t]]$values[[i]][positions]))
  }))
  wy <- lapply(seq_len(nt), function(t) c(1.2, -.6, .5, .9) + .08*(t-1))
  models <- as.matrix(expand.grid(rep(list(0:1), nt)))
  models <- lapply(seq_len(nrow(models)), function(i) as.integer(models[i, ]))
  pi <- c(.7, rep(.3/(length(models)-1), length(models)-1))
  b <- lapply(seq_len(nt), function(t) if (nonzero_b) c(.04*t, 0, -.02, 0) else rep(0, m))
  sets <- if (multiple_sets) list(c(0L, 2L), c(1L, 3L)) else list(0:(m-1L))
  mat_list <- function(x) split(x, rep(seq_len(ncol(x)), each = nrow(x)))
  common <- list(wy=wy, ww=ww, yy=rep(50, nt), b=b, sets=sets,
    B=diag(.15, nt), E=diag(.8, nt), ssb_prior=mat_list(diag(.05, nt)),
    sse_prior=mat_list(diag(.3, nt)), models=models, pi=pi, nub=4, nue=4,
    updateB=updates, updateE=updates, updatePi=updates, n=as.integer(40+seq_len(nt)),
    nit=8L, nburn=3L, nthin=2L, seed=17001L, method=4L)
  dense <- c(common, list(XXvalues=XXvalues, XXindices=XXindices))
  dense <- dense[c("wy","ww","yy","b","XXvalues","XXindices","sets","B","E",
                   "ssb_prior","sse_prior","models","pi","nub","nue","updateB",
                   "updateE","updatePi","n","nit","nburn","nthin","seed","method")]
  csr <- c(common, list(ld_prefixes=prefixes))
  csr <- csr[c("wy","ww","yy","b","ld_prefixes","sets","B","E","ssb_prior",
               "sse_prior","models","pi","nub","nue","updateB","updateE",
               "updatePi","n","nit","nburn","nthin","seed","method")]
  list(dense=dense, csr=csr, prefixes=prefixes)
}

phase17i_compare <- function(case, tolerance = 1e-12) {
  dense <- do.call(sblr:::mtblr, case$dense)
  csr <- do.call(sblr:::mtblr_csr_internal, case$csr)
  expect_equal(csr, dense, tolerance = tolerance)
  invisible(csr)
}

test_that("Phase 17I contracts are binding-neutral and trait-specific", {
  types <- phase17i_src("src/blr_mt_csr_types.h")
  access <- phase17i_src("src/blr_mt_ld_access.h")
  core <- phase17i_src("src/blr_mt_default_core_impl.h")
  expect_match(types, "struct MtSparseLdBundleView", fixed=TRUE)
  expect_match(types, "std::vector<sblr::core::SparseLdCsrView> trait_ld", fixed=TRUE)
  expect_match(types, "struct MtCsrDataView", fixed=TRUE)
  expect_false(grepl("Rcpp|SEXP|std::shared_ptr", types))
  expect_match(access, "std::vector<double>& residual", fixed=TRUE)
  expect_match(core, "run_mt_bayesc_core_impl", fixed=TRUE)
  expect_match(core, "run_mt_default_core", fixed=TRUE)
  expect_match(core, "run_mt_csr_core", fixed=TRUE)
  expect_identical(source_match_count("for ( int it =", core, fixed=TRUE), 1L)
})

test_that("shared CSR reductions match the dense numerical oracle", {
  phase17i_compare(phase17i_case(updates=FALSE))
  phase17i_compare(phase17i_case(updates=TRUE))
  phase17i_compare(phase17i_case(nonzero_b=TRUE, updates=TRUE))
  phase17i_compare(phase17i_case(multiple_sets=TRUE, updates=TRUE))
  phase17i_compare(phase17i_case(nt=3L, updates=TRUE))
})

test_that("trait-specific values and independent patterns match dense oracles", {
  phase17i_compare(phase17i_case(trait_specific=TRUE, updates=TRUE), 1e-12)
  phase17i_compare(phase17i_case(trait_specific=TRUE, independent=TRUE,
                                 updates=FALSE), 1e-12)
})

test_that("invalid shared and marker-domain inputs fail before sampling", {
  x <- phase17i_case(updates=FALSE)
  x$csr$ww[[2]][1] <- x$csr$ww[[2]][1] + .1
  expect_error(do.call(sblr:::mtblr_csr_internal, x$csr), "identical trait diagonals")
  x <- phase17i_case(trait_specific=TRUE, updates=FALSE)
  writeLines(c("n_variants=3", "nnz=0"), paste0(x$prefixes[2], ".meta.txt"))
  expect_error(do.call(sblr:::mtblr_csr_internal, x$csr), "marker dimension")
})

test_that("internal CSR scientific identities and trait residuals hold", {
  x <- phase17i_case(trait_specific=TRUE, independent=TRUE, updates=FALSE)
  fit <- do.call(sblr:::mtblr_csr_internal, x$csr)
  expect_true(all(is.finite(unlist(fit))))
  expect_true(all(unlist(fit[[6]]) %in% c(0, 1)))
  expect_true(all(unlist(fit[[2]]) >= 0 & unlist(fit[[2]]) <= 1))
  expect_equal(fit[[14]], unname(split(x$csr$B,
    rep(seq_len(ncol(x$csr$B)), each=nrow(x$csr$B)))))
  expect_equal(fit[[16]], unname(split(x$csr$E,
    rep(seq_len(ncol(x$csr$E)), each=nrow(x$csr$E)))))
  expect_equal(fit[[17]][[1]], x$csr$pi)
  for (trait in seq_along(x$dense$wy)) {
    residual <- x$dense$wy[[trait]]
    effects <- fit[[5]][[trait]]
    for (marker in seq_along(effects)) if (effects[marker] != 0) {
      for (j in seq_along(x$dense$XXindices[[marker]])) {
        index <- x$dense$XXindices[[marker]][j] + 1L
        residual[index] <- residual[index] -
          x$dense$XXvalues[[trait]][[marker]][j] * effects[marker]
      }
    }
    expect_equal(fit[[4]][[trait]], residual, tolerance=1e-12)
  }
})

test_that("one core, finalizer, adapter, and internal-only wrapper remain", {
  cpp <- phase17i_src("src/mtblr.cpp")
  adapter <- phase17i_src("src/blr_mt_default_legacy_adapter.h")
  rroute <- phase17i_src("R/interface_mtblr.R")
  namespace <- phase17i_src("NAMESPACE")
  expect_match(cpp, "run_mt_csr_core", fixed=TRUE)
  expect_match(cpp, "make_mt_default_legacy_result", fixed=TRUE)
  expect_match(adapter, "MtDefaultLegacyResult result(22)", fixed=TRUE)
  expect_false(grepl("mtblr_csr_internal", rroute, fixed=TRUE))
  expect_false(grepl("mtblr_csr_internal", namespace, fixed=TRUE))
  expect_false(grepl("mtblr_cpg_omp_csr", cpp, fixed=TRUE))
  expect_silent(yaml::read_yaml(blr_repo_path(".github", "workflows", "blr-framework.yml")))
  expect_silent(yaml::read_yaml(blr_repo_path(".github", "workflows", "blr-framework-extended.yml")))
})
