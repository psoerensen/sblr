phase17k_float32 <- function(x) {
  con <- rawConnection(raw(), "w+b")
  on.exit(close(con))
  writeBin(as.numeric(x), con, size = 4L, endian = .Platform$endian)
  seek(con, 0L)
  readBin(con, "numeric", n = length(x), size = 4L, endian = .Platform$endian)
}

phase17k_write_bed <- function(path, dosage) {
  code <- function(x) ifelse(is.na(x), 1L, c(`0` = 3L, `1` = 2L, `2` = 0L)[as.character(x)])
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    values <- as.integer(code(dosage[marker, ]))
    values <- c(values, rep(0L, (-length(values)) %% 4L))
    vapply(seq(1L, length(values), by = 4L), function(i) {
      sum(values[i:(i + 3L)] * c(1L, 4L, 16L, 64L))
    }, integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

phase17k_case <- function() {
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 1, 0),
    c(2, 1, 0, 2, 1, 0, NA, 2),
    c(0, 0, 1, 1, 2, 2, 1, 0),
    c(2, 2, 1, 1, 0, 0, 1, 2),
    c(0, 1, 0, 2, 2, 1, 1, NA)
  )
  path <- tempfile(fileext = ".bed")
  phase17k_write_bed(path, dosage)
  list(
    bed = path, dosage = dosage, af = c(0.25, 0.35, 0.4, 0.3, 0.45),
    blocks = c(0L, 1L, 3L), wy = rbind(c(1, -.5, .75, -.25, .4), c(.2, .7, -.1, .5, -.3)),
    effects = c(.1, 0, -.2, .05, .3)
  )
}

phase17k_oracle <- function(case, filter, tau = .01, eta = 0) {
  dmat <- function(x) diag(x, nrow = length(x), ncol = length(x))
  z <- case$dosage
  for (i in seq_len(nrow(z))) {
    p <- case$af[i]
    z[i, ] <- (z[i, ] - 2 * p) / sqrt(2 * p * (1 - p))
    z[i, is.na(z[i, ])] <- 0
  }
  z <- matrix(phase17k_float32(z), nrow(z), ncol(z))
  starts <- case$blocks + 1L
  ends <- c(starts[-1L] - 1L, nrow(z))
  packed <- vector("list", length(starts))
  diagonal <- numeric(nrow(z))
  transformed <- case$wy
  diagnostics <- vector("list", length(starts))
  for (g in seq_along(starts)) {
    ix <- starts[g]:ends[g]
    A <- tcrossprod(z[ix, , drop = FALSE])
    d <- diag(A)
    C <- A / outer(sqrt(d), sqrt(d))
    C <- .5 * (C + t(C))
    if (filter == "hard_truncate") {
      ee <- eigen(C, symmetric = TRUE)
      keep <- which(ee$values >= max(tau, .01))
      if (!length(keep)) keep <- which.max(ee$values)
      V <- ee$vectors[, keep, drop = FALSE]
      M <- dmat(sqrt(d)) %*% V %*% dmat(ee$values[keep]) %*% t(V) %*% dmat(sqrt(d))
      transformed[, ix] <- transformed[, ix, drop = FALSE] %*% dmat(1 / sqrt(d)) %*% V %*% t(V) %*% dmat(sqrt(d))
      diagnostics[[g]] <- c(n_kept = length(keep), mu_min = min(ee$values), shrink = 1 - length(keep) / length(ix))
    } else {
      if (filter == "ridge_fixed") {
        a <- min(1, max(0, eta / (1 + eta)))
      } else {
        C2 <- sum(C * C); d2 <- C2 - length(ix)
        q <- colSums(z[ix, , drop = FALSE]^2 / d)
        a <- if (d2 > 0) min(sum(q^2) - C2 / ncol(z), d2) / d2 else 0
        a <- min(1, max(0, a))
      }
      M <- (1 - a) * A + a * dmat(diag(A))
      diagnostics[[g]] <- c(n_kept = length(ix), mu_min = 0, shrink = a)
    }
    M <- matrix(phase17k_float32(M), nrow(M), ncol(M))
    packed[[g]] <- unlist(lapply(seq_len(nrow(M)), function(i) M[i, i:ncol(M)]), use.names = FALSE)
    diagonal[ix] <- diag(M)
  }
  list(packed = packed, diagonal = diagonal, wy = transformed, diagnostics = diagnostics)
}

phase17k_inspect <- function(case, filter = "hard_truncate", tau = .01, eta = 0, mutation = "") {
  sblr:::stblr_block_eigen_contract_internal(
    case$bed, ncol(case$dosage), list(seq_len(nrow(case$dosage))), NULL,
    case$af, case$blocks, case$wy, case$effects, filter, tau, eta, mutation
  )
}

test_that("canonical block-eigen inspection matches hard-truncation oracle", {
  case <- phase17k_case()
  for (tau in c(0, .01, .4, 100)) {
    actual <- phase17k_inspect(case, "hard_truncate", tau)
    expected <- phase17k_oracle(case, "hard_truncate", tau)
    expect_equal(actual$packed_upper_triangle, expected$packed, tolerance = 1e-10)
    expect_equal(as.numeric(actual$diagonal), expected$diagonal, tolerance = 1e-12)
    expect_equal(actual$transformed_wy, expected$wy, tolerance = 1e-10)
    expect_equal(as.numeric(actual$rebuilt_residual), actual$vector_rebuilt_residual, tolerance = 1e-12)
    expect_true(all(actual$diagnostics$n_kept >= 1L))
  }
})

test_that("canonical block-eigen inspection matches fixed-ridge oracle", {
  case <- phase17k_case()
  for (eta in c(0, 1, 1e12)) {
    actual <- phase17k_inspect(case, "ridge_fixed", eta = eta)
    expected <- phase17k_oracle(case, "ridge_fixed", eta = eta)
    expect_equal(actual$packed_upper_triangle, expected$packed, tolerance = 1e-12)
    expect_equal(as.numeric(actual$diagonal), expected$diagonal, tolerance = 1e-12)
    expect_identical(actual$transformed_wy, case$wy)
    expect_equal(actual$diagnostics$shrink, rep(eta / (1 + eta), 3), tolerance = 1e-12)
  }
})

test_that("canonical block-eigen inspection matches Ledoit-Wolf oracle", {
  case <- phase17k_case()
  actual <- phase17k_inspect(case, "ridge_lw")
  expected <- phase17k_oracle(case, "ridge_lw")
  expect_equal(actual$packed_upper_triangle, expected$packed, tolerance = 1e-10)
  expect_equal(as.numeric(actual$diagonal), expected$diagonal, tolerance = 1e-12)
  expect_identical(actual$transformed_wy, case$wy)
  expect_true(all(is.finite(actual$diagnostics$shrink)))
  expect_true(all(actual$diagnostics$shrink >= 0 & actual$diagnostics$shrink <= 1))
})

test_that("block-eigen builder rejects invalid domains and frequencies", {
  case <- phase17k_case()
  for (starts in list(integer(), 1L, c(0L, 0L), c(0L, 3L, 2L), c(0L, 5L))) {
    changed <- case; changed$blocks <- starts
    expect_error(phase17k_inspect(changed), "block_start|empty")
  }
  for (bad in c(NA_real_, NaN, Inf, 0, 1)) {
    changed <- case; changed$af[2] <- bad
    expect_error(phase17k_inspect(changed), "af values")
  }
  changed <- case; changed$af <- changed$af[-1]
  expect_error(phase17k_inspect(changed), "af length")
  changed <- case; changed$wy <- changed$wy[, -1, drop = FALSE]
  expect_error(phase17k_inspect(changed), "wy_mat")
})

test_that("central block-eigen validation rejects corrupted storage", {
  case <- phase17k_case()
  expect_error(phase17k_inspect(case, mutation = "mapping"), "mappings")
  expect_error(phase17k_inspect(case, mutation = "packed_length"), "packed triangle")
  expect_error(phase17k_inspect(case, mutation = "nonfinite"), "finite")
  expect_error(phase17k_inspect(case, mutation = "diagonal"), "diagonal")
})

test_that("block-eigen contract is binding-neutral and singular", {
  header <- blr_source_text("src/blr_block_eigen.h")
  operator <- blr_source_text("src/st_ld_operator.h")
  expect_match(header, "struct BlockEigenStorage", fixed = TRUE)
  expect_match(header, "struct BlockEigenView", fixed = TRUE)
  expect_match(header, "std::vector<double>& residual", fixed = TRUE)
  expect_false(grepl("Rcpp|SEXP", header))
  expect_match(operator, "using BlockEigenOperator = sblr::core::BlockEigenStorage", fixed = TRUE)
  expect_equal(sum(grepl("struct BlockEigenStorage", c(header, operator), fixed = TRUE)), 1)
})

test_that("all scalar block-eigen routes use the canonical builder and view-backed storage", {
  files <- vapply(c("src/st_cpg_omp_csr.cpp", "src/st_cpg_omp_csr_bayesr.cpp", "src/st_sbayesrc_omp_csr.cpp"), blr_source_text, character(1))
  expect_true(all(grepl("build_block_eigen", files, fixed = TRUE)))
  expect_true(all(grepl("BlockEigenOperator", files, fixed = TRUE)))
  expect_equal(sum(grepl("parse_block_eigen_filter_mode", files, fixed = TRUE)), 3)
  expect_equal(sum(grepl("block_eigen_diagnostics_to_data_frame", files, fixed = TRUE)), 3)
})

test_that("public CSR routing and retained MT eigen disposition remain unchanged", {
  r_text <- paste(vapply(c("R/sparse_ld_bed_helper.R", "R/stblr-csr-sbayesrc.R", "R/interface_mtblr.R"), blr_source_text, character(1)), collapse = "\n")
  mt_text <- blr_source_text("src/mtblr.cpp")
  expect_false("ld_backend" %in% names(formals(sblr::stblr_csr)))
  expect_false("ld_backend" %in% names(formals(sblr::stblr_csr_bayesr)))
  expect_false("ld_backend" %in% names(formals(sblr::stblr_csr_annot)))
  expect_false("mtblr_eigen" %in% getNamespaceExports("sblr"))
  expect_match(mt_text, "mtblr_eigen", fixed = TRUE)
  expect_false(grepl("run_mt_bayesc_core_impl.*mtblr_eigen", mt_text))
})
