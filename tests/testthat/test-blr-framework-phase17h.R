phase17h_src <- function(file) {
  paste(readLines(blr_repo_path("src", file), warn = FALSE),
        collapse = "\n")
}

test_that("one binding-neutral sparse-LD storage and view contract is active", {
  shared <- phase17h_src("blr_sparse_ld_csr.h")
  common <- phase17h_src("st_csr_common.h")
  types <- phase17h_src("blr_csr_bayesc_types.h")

  expect_match(shared, "struct SparseLdCsrStorage", fixed = TRUE)
  expect_match(shared, "struct SparseLdCsrView", fixed = TRUE)
  expect_match(common,
               "using STLDCSR = sblr::core::SparseLdCsrStorage", fixed = TRUE)
  expect_false(grepl("struct STLDCSR", common, fixed = TRUE))
  expect_false(grepl("Rcpp", shared, fixed = TRUE))
  expect_false(grepl("SEXP", shared, fixed = TRUE))
  expect_match(shared, "std::vector<std::uint64_t> ptr", fixed = TRUE)
  expect_match(shared, "std::vector<int> idx", fixed = TRUE)
  expect_match(shared, "std::vector<float> xij", fixed = TRUE)
  expect_match(shared, "const std::uint64_t* row_ptr", fixed = TRUE)
  expect_match(shared, "const int* column_index", fixed = TRUE)
  expect_match(shared, "const float* offdiag_xij", fixed = TRUE)
  expect_match(shared, "const arma::rowvec* diagonal", fixed = TRUE)
  expect_false(grepl("trait_count", shared, fixed = TRUE))
  expect_false(grepl("sample_size", shared, fixed = TRUE))
  expect_false(grepl("std::mt19937", shared, fixed = TRUE))
  expect_false(grepl("Chain", shared, fixed = TRUE))
  expect_match(types, "SparseLdCsrView ld", fixed = TRUE)
  expect_false(grepl("const float* values", types, fixed = TRUE))
})

test_that("ordinary canonical CSR BayesC composes and validates the shared view", {
  adapter <- phase17h_src("st_cpg_omp_csr.cpp")
  core <- phase17h_src("blr_csr_bayesc_core_impl.h")
  operator <- phase17h_src("st_ld_operator.h")

  expect_match(adapter, "input.data.ld = op.view()", fixed = TRUE)
  expect_match(operator, "SparseLdCsrView view() const", fixed = TRUE)
  expect_match(core, "validate_sparse_ld_csr_view(input.data.ld)", fixed = TRUE)
  expect_match(core, "op.ld.apply_offdiag", fixed = TRUE)
  expect_match(core, "input.data.ld.rebuild", fixed = TRUE)
})

test_that("the canonical builder preserves disk and numerical conventions", {
  builder <- phase17h_src("st_csr_common.h")
  expect_match(builder, ".row_ptr.u64.bin", fixed = TRUE)
  expect_match(builder, ".col_idx.u32.0based.bin", fixed = TRUE)
  expect_match(builder, ".values.f32.bin", fixed = TRUE)
  expect_match(builder, "if (j == i) continue", fixed = TRUE)
  expect_match(builder, "++degree[static_cast<std::size_t>(i)]", fixed = TRUE)
  expect_match(builder, "++degree[static_cast<std::size_t>(j)]", fixed = TRUE)
  expect_match(builder,
    "rij * std::sqrt(xx[static_cast<std::size_t>(i)] * xx[static_cast<std::size_t>(j)])",
    fixed = TRUE)
  expect_match(builder, "const float xij_f = static_cast<float>(xij)", fixed = TRUE)
})

test_that("the centralized borrowed-view validator covers malformed views", {
  shared <- phase17h_src("blr_sparse_ld_csr.h")
  required <- c(
    "marker_count == 0", "row_ptr == nullptr",
    "row_ptr_size != view.marker_count + 1", "view.row_ptr[0] != 0",
    "view.row_ptr[row] > view.row_ptr[row + 1]",
    "view.row_ptr[view.marker_count] != view.nonzero_count",
    "view.column_index == nullptr", "view.offdiag_xij == nullptr",
    "column < 0", "static_cast<std::size_t>(column) == row",
    "!std::isfinite(static_cast<double>(view.offdiag_xij[p]))",
    "view.diagonal == nullptr", "view.diagonal->n_elem != view.marker_count",
    "!std::isfinite(value) || value <= 0.0"
  )
  for (condition in required) {
    expect_match(shared, condition, fixed = TRUE, info = condition)
  }
})

test_that("future trait-specific sharing is representable without a bundle type", {
  shared <- phase17h_src("blr_sparse_ld_csr.h")
  research <- phase17h_src("mt_cpg_omp_csr.cpp")
  route <- phase17h_src("../R/interface_mtblr.R")

  expect_false(grepl("MtSparseLdBundleView", shared, fixed = TRUE))
  expect_false(grepl("trait_count", shared, fixed = TRUE))
  expect_false(grepl("trait_id", shared, fixed = TRUE))
  expect_match(research, "INTERNAL RESEARCH ONLY", fixed = TRUE)
  expect_match(research, "noncanonical", fixed = TRUE)
  expect_false(grepl("mtblr_cpg_omp_csr", route, fixed = TRUE))
})
