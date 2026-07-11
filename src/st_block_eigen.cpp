// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <RcppArmadillo.h>

#include "st_bed_decode.h"
#include "st_block_eigen.h"

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
  B.tri.resize(static_cast<std::size_t>(s) * (s + 1) / 2);

  for (int i = 0; i < s; ++i) {
    for (int j = i; j < s; ++j) {
      B.tri[B.pidx(i, j)] =
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
  op.blocks.reserve(static_cast<std::size_t>(nb));
  op.block_of.assign(static_cast<std::size_t>(m), -1);
  op.local_of.assign(static_cast<std::size_t>(m), -1);
  op.diag_.set_size(static_cast<arma::uword>(m));

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
      op.diag_(static_cast<arma::uword>(gi)) = stored.sym_at(i, i);
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

  return op;
}
