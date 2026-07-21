// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <RcppArmadillo.h>

#include "st_bed_decode.h"
#include "st_block_eigen.h"
#include "st_block_eigen_rcpp.h"

static double ledoit_wolf_lambda(
    const float* Z,
    int ma,
    int n,
    const std::vector<double>& xx_blk,
    const arma::mat& C) {
  const double p = static_cast<double>(ma);
  const double nf = static_cast<double>(n);
  const double C2 = arma::accu(C % C);
  const double d2 = C2 - p;

  if (!(d2 > 0.0)) return 0.0;

  double sum_q2 = 0.0;

  for (int k = 0; k < n; ++k) {
    double q = 0.0;

    for (int ii = 0; ii < ma; ++ii) {
      const double xxi = xx_blk[static_cast<std::size_t>(ii)];
      if (!(xxi > 0.0)) continue;

      const double z =
        static_cast<double>(Z[static_cast<std::size_t>(ii) * n + k]);
      q += z * z / xxi;
    }

    sum_q2 += q * q;
  }

  const double bbar = sum_q2 - C2 / nf;
  double a = std::min(bbar, d2) / d2;

  if (a < 0.0) a = 0.0;
  if (a > 1.0) a = 1.0;

  return a;
}

static void pack_upper(const arma::mat& M, EigenBlock& B) {
  const int s = B.size;
  B.upper_triangle.resize(static_cast<std::size_t>(s) * (s + 1) / 2);

  for (int i = 0; i < s; ++i) {
    for (int j = i; j < s; ++j) {
      B.upper_triangle[B.pidx(i, j)] =
        static_cast<float>(
          M(static_cast<arma::uword>(i), static_cast<arma::uword>(j))
        );
    }
  }
}

BlockEigenOperator build_block_eigen(
    const PackedBedMatrix& G,
    const std::vector<double>& af,
    const std::vector<int>& block_start,
    EigenFilterMode mode,
    double tau,
    double eta,
    arma::mat& wy_mat,
    int nthreads,
    std::vector<BlockEigenDiag>* diag_out) {
  const int m = G.m;
  const int n = G.n;
  const int nb = static_cast<int>(block_start.size());

  if (block_start.empty() || block_start[0] != 0) {
    throw std::runtime_error("build_block_eigen: block_start must start at 0.");
  }

  if (static_cast<int>(af.size()) != G.m) {
    throw std::runtime_error("build_block_eigen: af length != m.");
  }

  for (double frequency : af) {
    if (!std::isfinite(frequency) || frequency <= 0.0 || frequency >= 1.0) {
      throw std::runtime_error("build_block_eigen: af values must be finite and in (0, 1).");
    }
  }

  if (wy_mat.n_cols != static_cast<arma::uword>(G.m)) {
    throw std::runtime_error("build_block_eigen: wy_mat must be nt x m.");
  }

  for (int b = 1; b < nb; ++b) {
    if (block_start[b] <= block_start[b - 1] || block_start[b] > m) {
      throw std::runtime_error(
        "build_block_eigen: block_start must be ascending and within marker range."
      );
    }
  }

  const std::vector<double> xx =
    compute_xx_from_packed_standardized(G, af.data(), nthreads);
  const double mu_floor = 0.01;

  BlockEigenOperator op;
  op.marker_count = static_cast<std::size_t>(m);
  op.blocks.reserve(static_cast<std::size_t>(nb));
  op.block_of.assign(static_cast<std::size_t>(m), -1);
  op.local_of.assign(static_cast<std::size_t>(m), -1);
  op.diagonal.set_size(static_cast<arma::uword>(m));

  std::vector<float> Z;

  for (int b = 0; b < nb; ++b) {
    const int start = block_start[b];
    const int end = (b + 1 < nb) ? block_start[b + 1] : m;
    const int ma = end - start;

    if (ma <= 0) {
      throw std::runtime_error("build_block_eigen: empty LD block.");
    }

    Z.assign(static_cast<std::size_t>(ma) * n, 0.0f);
    decode_packed_block_float(G, start, ma, af.data(), Z.data(), nthreads);

    arma::mat A(
      static_cast<arma::uword>(ma),
      static_cast<arma::uword>(ma),
      arma::fill::zeros
    );

    for (int ia = 0; ia < ma; ++ia) {
      const float* za = Z.data() + static_cast<std::size_t>(ia) * n;

      for (int ib = ia; ib < ma; ++ib) {
        const float* zb = Z.data() + static_cast<std::size_t>(ib) * n;
        double acc = 0.0;

        for (int k = 0; k < n; ++k) {
          acc += static_cast<double>(za[k]) * static_cast<double>(zb[k]);
        }

        A(static_cast<arma::uword>(ia), static_cast<arma::uword>(ib)) = acc;
        A(static_cast<arma::uword>(ib), static_cast<arma::uword>(ia)) = acc;
      }
    }

    arma::vec d = A.diag();
    arma::vec dsqrt(static_cast<arma::uword>(ma));
    arma::vec dinvsqrt(static_cast<arma::uword>(ma));
    std::vector<double> xx_blk(static_cast<std::size_t>(ma));

    for (int i = 0; i < ma; ++i) {
      const double di = d(static_cast<arma::uword>(i));
      const double xxi = xx[static_cast<std::size_t>(start + i)];

      if (!std::isfinite(di) || di <= 0.0 || !std::isfinite(xxi) || xxi <= 0.0) {
        throw std::runtime_error(
          "build_block_eigen: block contains non-positive standardized diagonal."
        );
      }

      dsqrt(static_cast<arma::uword>(i)) = std::sqrt(di);
      dinvsqrt(static_cast<arma::uword>(i)) = 1.0 / dsqrt(static_cast<arma::uword>(i));
      xx_blk[static_cast<std::size_t>(i)] = di;
    }

    arma::mat C = A;
    C.each_col() %= dinvsqrt;
    C.each_row() %= dinvsqrt.t();
    C = 0.5 * (C + C.t());

    arma::mat tildeA;
    arma::mat Lk;
    arma::mat Rk;
    bool do_project = false;

    BlockEigenDiag dg;
    dg.start = start;
    dg.size = ma;
    dg.n_kept = ma;

    if (mode == EigenFilterMode::hard_truncate) {
      arma::vec mu;
      arma::mat V;

      if (!arma::eig_sym(mu, V, C)) {
        throw std::runtime_error("build_block_eigen: eig_sym failed.");
      }

      dg.mu_min = mu.min();

      const double threshold = std::max(tau, mu_floor);
      arma::uvec keep = arma::find(mu >= threshold);

      if (keep.n_elem == 0) {
        keep = arma::uvec(1);
        keep(0) = mu.index_max();
      }

      arma::mat Vk = V.cols(keep);
      arma::vec mk = mu.elem(keep);
      arma::mat Wk = Vk;
      Wk.each_col() %= dsqrt;

      tildeA = Wk * arma::diagmat(mk) * Wk.t();

      Lk = Vk;
      Lk.each_col() %= dinvsqrt;
      Rk = Vk.t();
      Rk.each_row() %= dsqrt.t();

      do_project = true;
      dg.n_kept = static_cast<int>(keep.n_elem);
      dg.shrink =
        1.0 - static_cast<double>(keep.n_elem) / static_cast<double>(ma);
    } else {
      double a = 0.0;

      if (mode == EigenFilterMode::ridge_lw) {
        a = ledoit_wolf_lambda(Z.data(), ma, n, xx_blk, C);
      } else if (mode == EigenFilterMode::ridge_fixed) {
        a = eta / (1.0 + eta);
      } else {
        throw std::runtime_error("build_block_eigen: unsupported filter mode.");
      }

      if (!std::isfinite(a)) {
        throw std::runtime_error("build_block_eigen: shrinkage weight is invalid.");
      }

      if (a < 0.0) a = 0.0;
      if (a > 1.0) a = 1.0;

      tildeA = (1.0 - a) * A;
      tildeA.diag() += a * d;

      dg.shrink = a;
    }

    EigenBlock eb;
    eb.start = start;
    eb.size = ma;
    pack_upper(tildeA, eb);

    op.blocks.push_back(std::move(eb));
    const int gid = static_cast<int>(op.blocks.size()) - 1;
    const EigenBlock& stored = op.blocks.back();

    for (int i = 0; i < ma; ++i) {
      const int gi = start + i;
      op.block_of[static_cast<std::size_t>(gi)] = gid;
      op.local_of[static_cast<std::size_t>(gi)] = i;
      op.diagonal(static_cast<arma::uword>(gi)) = stored.sym_at(i, i);
    }

    if (do_project && wy_mat.n_rows > 0) {
      const arma::uword c0 = static_cast<arma::uword>(start);
      const arma::uword c1 = static_cast<arma::uword>(end - 1);
      arma::mat W = wy_mat.cols(c0, c1);
      wy_mat.cols(c0, c1) = (W * Lk) * Rk;
    }

    if (diag_out != nullptr) {
      diag_out->push_back(dg);
    }
  }

  sblr::core::validate_block_eigen_storage(op);
  return op;
}

// Maintenance-only inspection of the canonical builder and borrowed view.
// [[Rcpp::export]]
Rcpp::List stblr_block_eigen_contract_internal(
    Rcpp::CharacterVector bed_files,
    int n_bed,
    Rcpp::List cls,
    Rcpp::Nullable<Rcpp::IntegerVector> rows,
    Rcpp::NumericVector af,
    Rcpp::IntegerVector block_start,
    Rcpp::NumericMatrix wy,
    Rcpp::NumericVector effects,
    std::string eigen_filter = "hard_truncate",
    double eigen_tau = 0.01,
    double eigen_eta = 0.0,
    std::string validation_mutation = "") {
  std::vector<std::string> files = Rcpp::as<std::vector<std::string>>(bed_files);
  std::vector<std::vector<int>> columns(static_cast<std::size_t>(cls.size()));
  for (int group = 0; group < cls.size(); ++group) {
    Rcpp::IntegerVector values = cls[group];
    columns[static_cast<std::size_t>(group)].reserve(static_cast<std::size_t>(values.size()));
    for (int value : values) {
      if (value == NA_INTEGER || value <= 0) throw std::runtime_error("cls must be positive and 1-based.");
      columns[static_cast<std::size_t>(group)].push_back(value);
    }
  }
  std::vector<int> selected_rows;
  if (rows.isNotNull()) {
    Rcpp::IntegerVector values(rows);
    selected_rows.reserve(static_cast<std::size_t>(values.size()));
    for (int value : values) {
      if (value == NA_INTEGER || value <= 0 || value > n_bed)
        throw std::runtime_error("rows must lie in [1, n_bed].");
      selected_rows.push_back(value - 1);
    }
  }
  PackedBedMatrix packed = read_bedfiles_to_packed_matrix(
    files, n_bed, selected_rows.empty() ? nullptr : selected_rows.data(),
    static_cast<int>(selected_rows.size()), columns
  );
  std::vector<double> frequencies = Rcpp::as<std::vector<double>>(af);
  std::vector<int> starts = Rcpp::as<std::vector<int>>(block_start);
  arma::mat transformed_wy(wy.begin(), wy.nrow(), wy.ncol(), true);
  std::vector<BlockEigenDiag> diagnostics;
  BlockEigenOperator storage = build_block_eigen(
    packed, frequencies, starts, parse_block_eigen_filter_mode(eigen_filter),
    eigen_tau, eigen_eta, transformed_wy, 1, &diagnostics
  );
  if (validation_mutation == "mapping") storage.local_of[0] = 1;
  else if (validation_mutation == "packed_length") storage.blocks[0].upper_triangle.pop_back();
  else if (validation_mutation == "nonfinite") storage.blocks[0].upper_triangle[0] = std::numeric_limits<float>::quiet_NaN();
  else if (validation_mutation == "diagonal") storage.diagonal(0) = 0.0;
  else if (!validation_mutation.empty()) throw std::runtime_error("unknown validation mutation.");
  const sblr::core::BlockEigenView view = storage.view();
  sblr::core::validate_block_eigen_view(view);
  if (effects.size() != static_cast<int>(storage.marker_count))
    throw std::runtime_error("effects length must equal marker count.");
  arma::rowvec effect_row(effects.begin(), effects.size(), true);
  arma::rowvec residual;
  view.rebuild(transformed_wy.row(0), effect_row, residual);
  std::vector<double> wy_vector(transformed_wy.n_cols);
  std::vector<double> effect_vector(static_cast<std::size_t>(effects.size()));
  for (arma::uword i = 0; i < transformed_wy.n_cols; ++i) wy_vector[i] = transformed_wy(0, i);
  for (int i = 0; i < effects.size(); ++i) effect_vector[static_cast<std::size_t>(i)] = effects[i];
  std::vector<double> vector_residual;
  view.rebuild(wy_vector, effect_vector, vector_residual);

  Rcpp::List triangles(storage.blocks.size());
  Rcpp::IntegerVector starts_out(storage.blocks.size()), sizes_out(storage.blocks.size());
  for (std::size_t block = 0; block < storage.blocks.size(); ++block) {
    const EigenBlock& item = storage.blocks[block];
    starts_out[block] = item.start;
    sizes_out[block] = item.size;
    triangles[block] = Rcpp::wrap(item.upper_triangle);
  }
  return Rcpp::List::create(
    Rcpp::Named("marker_count") = static_cast<int>(storage.marker_count),
    Rcpp::Named("block_start") = starts_out,
    Rcpp::Named("block_size") = sizes_out,
    Rcpp::Named("packed_upper_triangle") = triangles,
    Rcpp::Named("block_of") = Rcpp::wrap(storage.block_of),
    Rcpp::Named("local_of") = Rcpp::wrap(storage.local_of),
    Rcpp::Named("diagonal") = Rcpp::wrap(storage.diagonal),
    Rcpp::Named("transformed_wy") = Rcpp::wrap(transformed_wy),
    Rcpp::Named("rebuilt_residual") = Rcpp::wrap(residual),
    Rcpp::Named("vector_rebuilt_residual") = Rcpp::wrap(vector_residual),
    Rcpp::Named("diagnostics") = block_eigen_diagnostics_to_data_frame(diagnostics),
    Rcpp::Named("filter") = Rcpp::List::create(
      Rcpp::Named("mode") = eigen_filter,
      Rcpp::Named("tau") = eigen_tau,
      Rcpp::Named("eta") = eigen_eta,
      Rcpp::Named("mu_floor") = 0.01
    )
  );
}
