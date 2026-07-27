#ifndef SBLR_BLR_MT_BED_CORE_IMPL_H
#define SBLR_BLR_MT_BED_CORE_IMPL_H

#include "blr_mt_bed_access.h"
#include "blr_mt_covariance_rng.h"

#include <armadillo>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr {
namespace mt {

inline bool mt_bed_is_symmetric(const arma::mat& x, const double tol = 1e-12) {
  return x.n_rows == x.n_cols &&
    arma::approx_equal(x, x.t(), "absdiff", tol);
}

inline void validate_mt_bed_spd(
    const arma::mat& x,
    const std::size_t nt,
    const char* label) {
  if (x.n_rows != nt || x.n_cols != nt || !x.is_finite() ||
      !mt_bed_is_symmetric(x)) {
    throw std::invalid_argument(std::string(label) +
      " must be a finite symmetric nt by nt matrix");
  }
  arma::mat lower;
  if (!arma::chol(lower, x, "lower")) {
    throw std::invalid_argument(std::string(label) +
      " must be symmetric positive definite");
  }
}

inline bool mt_bed_is_exactly_diagonal(const arma::mat& x) {
  for (arma::uword col = 0; col < x.n_cols; ++col) {
    for (arma::uword row = 0; row < x.n_rows; ++row) {
      if (row != col && x(row, col) != 0.0) {
        return false;
      }
    }
  }
  return true;
}

inline arma::mat mt_bed_prior_matrix(
    const std::vector<std::vector<double>>& x,
    const std::size_t nt,
    const char* label) {
  if (x.size() != nt) {
    throw std::invalid_argument(std::string(label) +
      " must have nt rows");
  }
  arma::mat out(nt, nt);
  for (std::size_t row = 0; row < nt; ++row) {
    if (x[row].size() != nt) {
      throw std::invalid_argument(std::string(label) +
        " must be nt by nt");
    }
    for (std::size_t col = 0; col < nt; ++col) {
      out(row, col) = x[row][col];
    }
  }
  return out;
}

inline bool mt_bed_state_is_model(
    const std::vector<std::vector<int>>& state,
    const std::size_t marker,
    const std::vector<std::vector<int>>& models) {
  for (const auto& model : models) {
    bool same = true;
    for (std::size_t trait = 0; trait < model.size(); ++trait) {
      if (state[trait][marker] != model[trait]) {
        same = false;
        break;
      }
    }
    if (same) {
      return true;
    }
  }
  return false;
}

template <class PackedGenotype>
inline void validate_mt_bed_problem(
    const MtBedDataView<PackedGenotype>& data,
    const MtBedInitialState& initial,
    const std::vector<std::vector<int>>& sets,
    const std::vector<std::vector<double>>& ssb_prior,
    const std::vector<std::vector<double>>& sse_prior,
    const std::vector<std::vector<int>>& models,
    const MtBedExecutionSpec& execution,
    const double nub,
    const double nue) {
  validate_mt_bed_genotype_view(data.genotype);
  const std::size_t m = data.genotype.marker_count;
  const std::size_t n = data.genotype.sample_count;
  if (data.marker_maps.size() != m) {
    throw std::invalid_argument("marker-map count must equal marker count");
  }
  for (const auto& map : data.marker_maps) {
    for (double value : map.value) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument("marker-map values must be finite");
      }
    }
    if (!std::isfinite(map.xx) || map.xx <= 0.0) {
      throw std::invalid_argument("marker-map xx must be finite and positive");
    }
  }

  if (data.phenotype.n_rows != n || data.phenotype.n_cols == 0 ||
      !data.phenotype.is_finite()) {
    throw std::invalid_argument(
      "phenotype must be a finite sample_count by nt matrix");
  }
  const std::size_t nt = data.phenotype.n_cols;
  for (std::size_t trait = 0; trait < nt; ++trait) {
    long double sum = 0.0L;
    long double sumsq = 0.0L;
    for (std::size_t sample = 0; sample < n; ++sample) {
      const double value = data.phenotype(sample, trait);
      sum += static_cast<long double>(value);
      sumsq += static_cast<long double>(value) *
        static_cast<long double>(value);
    }
    const double mean = static_cast<double>(sum / n);
    const double rms = std::sqrt(static_cast<double>(sumsq / n));
    if (std::abs(mean) > 1e-10 * std::max(1.0, rms)) {
      throw std::invalid_argument(
        "phenotype columns must be centered before native execution");
    }
  }
  if (data.marker_wy.n_rows != m || data.marker_wy.n_cols != nt ||
      !data.marker_wy.is_finite()) {
    throw std::invalid_argument("marker_wy must be a finite m by nt matrix");
  }
  if (data.marker_order.size() != m) {
    throw std::invalid_argument("marker order must cover every marker");
  }
  std::vector<int> seen_order(m, 0);
  for (int marker : data.marker_order) {
    if (marker < 0 || static_cast<std::size_t>(marker) >= m ||
        seen_order[marker]++) {
      throw std::invalid_argument(
        "marker order must be a permutation of zero-based marker indices");
    }
  }

  if (models.empty()) {
    throw std::invalid_argument("models must be nonempty");
  }
  bool has_null = false;
  std::set<std::vector<int>> unique_models;
  for (const auto& model : models) {
    if (model.size() != nt) {
      throw std::invalid_argument("every model must have width nt");
    }
    bool is_null = true;
    for (int value : model) {
      if (value != 0 && value != 1) {
        throw std::invalid_argument("model entries must be binary");
      }
      is_null = is_null && value == 0;
    }
    has_null = has_null || is_null;
    if (!unique_models.insert(model).second) {
      throw std::invalid_argument("model patterns must be unique");
    }
  }
  if (!has_null) {
    throw std::invalid_argument("models must include the null pattern");
  }
  if (initial.pi.size() != models.size()) {
    throw std::invalid_argument(
      "model-probability length must equal model count");
  }
  double pi_total = 0.0;
  for (double value : initial.pi) {
    if (!std::isfinite(value) || value < 0.0) {
      throw std::invalid_argument(
        "model probabilities must be finite and nonnegative");
    }
    pi_total += value;
  }
  if (!std::isfinite(pi_total) || pi_total <= 0.0 ||
      std::abs(pi_total - 1.0) > 1e-12) {
    throw std::invalid_argument(
      "model probabilities must be normalized before core execution");
  }

  if (sets.empty()) {
    throw std::invalid_argument("sets must be nonempty");
  }
  std::vector<int> marker_membership(m, 0);
  for (const auto& set : sets) {
    if (set.empty()) {
      throw std::invalid_argument("every set must be nonempty");
    }
    std::set<int> within_set;
    for (int marker : set) {
      if (marker < 0 || static_cast<std::size_t>(marker) >= m) {
        throw std::invalid_argument("set marker index is out of range");
      }
      if (!within_set.insert(marker).second || marker_membership[marker]++) {
        throw std::invalid_argument(
          "sets must be a disjoint marker partition");
      }
    }
  }
  if (std::find(marker_membership.begin(), marker_membership.end(), 0) !=
      marker_membership.end()) {
    throw std::invalid_argument("sets must cover every marker");
  }

  validate_mt_bed_spd(initial.B, nt, "B");
  validate_mt_bed_spd(initial.E, nt, "E");
  const arma::mat ssb = mt_bed_prior_matrix(ssb_prior, nt, "ssb_prior");
  const arma::mat sse = mt_bed_prior_matrix(sse_prior, nt, "sse_prior");
  validate_mt_bed_spd(ssb, nt, "ssb_prior");
  validate_mt_bed_spd(sse, nt, "sse_prior");
  if (!std::isfinite(nub) || nub <= static_cast<double>(nt) - 1.0) {
    throw std::invalid_argument("nub must be finite and greater than nt - 1");
  }
  if (!std::isfinite(nue) || nue <= static_cast<double>(nt) - 1.0) {
    throw std::invalid_argument("nue must be finite and greater than nt - 1");
  }
  if (execution.residual_covariance != "full" &&
      execution.residual_covariance != "diagonal") {
    throw std::invalid_argument(
      "residual_covariance must be 'full' or 'diagonal'");
  }
  if (execution.residual_covariance == "diagonal" &&
      (!mt_bed_is_exactly_diagonal(initial.E) ||
       !mt_bed_is_exactly_diagonal(sse))) {
    throw std::invalid_argument(
      "diagonal residual-covariance mode requires diagonal E and sse_prior");
  }

  if (initial.beta.size() != nt || initial.b.size() != nt ||
      initial.state.size() != nt) {
    throw std::invalid_argument(
      "beta, b, and state must each have nt rows");
  }
  for (std::size_t trait = 0; trait < nt; ++trait) {
    if (initial.beta[trait].size() != m ||
        initial.b[trait].size() != m ||
        initial.state[trait].size() != m) {
      throw std::invalid_argument(
        "beta, b, and state must each be nt by m");
    }
    for (std::size_t marker = 0; marker < m; ++marker) {
      const double latent = initial.beta[trait][marker];
      const double effective = initial.b[trait][marker];
      const int active = initial.state[trait][marker];
      if (!std::isfinite(latent) || !std::isfinite(effective)) {
        throw std::invalid_argument(
          "initial marker effects must be finite");
      }
      if (active != 0 && active != 1) {
        throw std::invalid_argument("initial state must be binary");
      }
      if (active == 0 && effective != 0.0) {
        throw std::invalid_argument(
          "inactive effective marker effects must be exactly zero");
      }
      if (active == 1 &&
          std::abs(effective - latent) >
            1e-12 * std::max({1.0, std::abs(effective), std::abs(latent)})) {
        throw std::invalid_argument(
          "active effective and latent marker effects must agree");
      }
    }
  }
  for (std::size_t marker = 0; marker < m; ++marker) {
    if (!mt_bed_state_is_model(initial.state, marker, models)) {
      throw std::invalid_argument(
        "every initial marker-state column must equal a supplied model");
    }
  }

  if (execution.method != 4) {
    throw std::invalid_argument("mtblr_bed_internal supports method = 4 only");
  }
  if (execution.nit <= 0 || execution.nburn < 0 || execution.nthin <= 0) {
    throw std::invalid_argument(
      "nit and nthin must be positive and nburn must be nonnegative");
  }
}

inline arma::mat mt_bed_cholesky_with_safeguard(
    arma::mat& precision,
    MtBedCoreDiagnostics* diagnostics) {
  arma::mat lower;
  if (arma::chol(lower, precision, "lower")) {
    return lower;
  }
  double increment = 1e-8;
  for (std::size_t attempt = 0; attempt < 7; ++attempt) {
    precision.diag() += increment;
    if (diagnostics != nullptr) {
      ++diagnostics->marker_cholesky_jitter_attempts;
      diagnostics->marker_cholesky_max_increment = std::max(
        diagnostics->marker_cholesky_max_increment, increment);
    }
    if (arma::chol(lower, precision, "lower")) {
      return lower;
    }
    increment *= 10.0;
  }
  throw std::runtime_error(
    "marker posterior precision is not positive definite after jitter");
}

inline MtBedMarkerKernelResult mt_bed_marker_kernel(
    const arma::vec& score,
    const double xx,
    const arma::mat& B_inverse,
    const arma::mat& E_inverse,
    const std::vector<std::vector<int>>& models,
    const std::vector<double>& pi,
    MtBedCoreDiagnostics* diagnostics = nullptr) {
  const std::size_t nt = score.n_elem;
  if (!std::isfinite(xx) || xx <= 0.0 ||
      B_inverse.n_rows != nt || B_inverse.n_cols != nt ||
      E_inverse.n_rows != nt || E_inverse.n_cols != nt ||
      !score.is_finite() || !B_inverse.is_finite() ||
      !E_inverse.is_finite() || models.empty() ||
      models.size() != pi.size()) {
    throw std::invalid_argument("invalid MT BED marker-kernel input");
  }

  MtBedMarkerKernelResult out;
  out.models.resize(models.size());
  out.log_weight.assign(
    models.size(), -std::numeric_limits<double>::infinity());
  out.probability.assign(models.size(), 0.0);
  for (std::size_t k = 0; k < models.size(); ++k) {
    if (models[k].size() != nt) {
      throw std::invalid_argument("marker-kernel model width must equal nt");
    }
    arma::mat D(nt, nt, arma::fill::zeros);
    for (std::size_t trait = 0; trait < nt; ++trait) {
      if (models[k][trait] != 0 && models[k][trait] != 1) {
        throw std::invalid_argument("marker-kernel models must be binary");
      }
      D(trait, trait) = static_cast<double>(models[k][trait]);
    }
    MtBedMarkerKernelModel& model = out.models[k];
    model.precision = B_inverse + xx * D * E_inverse * D;
    model.rhs = D * E_inverse * score;
    model.lower = mt_bed_cholesky_with_safeguard(
      model.precision, diagnostics);
    model.mean = arma::solve(
      arma::trimatu(model.lower.t()),
      arma::solve(arma::trimatl(model.lower), model.rhs));
    arma::mat inverse_lower = arma::solve(
      arma::trimatl(model.lower),
      arma::eye<arma::mat>(nt, nt));
    model.covariance = inverse_lower.t() * inverse_lower;
    model.covariance = 0.5 * (model.covariance + model.covariance.t());
    if (std::isfinite(pi[k]) && pi[k] > 0.0) {
      const double log_det_half =
        arma::sum(arma::log(model.lower.diag()));
      out.log_weight[k] = std::log(pi[k]) - log_det_half +
        0.5 * arma::dot(model.rhs, model.mean);
    }
  }

  const double max_log = *std::max_element(
    out.log_weight.begin(), out.log_weight.end());
  if (!std::isfinite(max_log)) {
    throw std::runtime_error(
      "marker model probabilities cannot be normalized");
  }
  double total = 0.0;
  for (std::size_t k = 0; k < models.size(); ++k) {
    out.probability[k] = std::isfinite(out.log_weight[k]) ?
      std::exp(out.log_weight[k] - max_log) : 0.0;
    total += out.probability[k];
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::runtime_error(
      "marker model probabilities cannot be normalized");
  }
  for (double& value : out.probability) {
    value /= total;
  }
  return out;
}

inline std::size_t mt_bed_draw_model(
    const std::vector<double>& probability,
    std::mt19937& rng) {
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  const double draw = uniform(rng);
  double cumulative = 0.0;
  for (std::size_t k = 0; k < probability.size(); ++k) {
    cumulative += probability[k];
    if (draw < cumulative || k + 1 == probability.size()) {
      return k;
    }
  }
  throw std::runtime_error("failed to select a marker model");
}

inline arma::mat mt_bed_genetic_covariance(
    const arma::mat& phenotype,
    const arma::mat& residual) {
  const arma::mat genetic_value = phenotype - residual;
  arma::mat covariance =
    genetic_value.t() * genetic_value /
    static_cast<double>(phenotype.n_rows);
  return 0.5 * (covariance + covariance.t());
}

template <class PackedGenotype>
inline MtBedCoreResult run_mt_bed_bayesc_core(
    const MtBedDataView<PackedGenotype>& data,
    const MtBedInitialState& initial,
    const std::vector<std::vector<int>>& sets,
    const std::vector<std::vector<double>>& ssb_prior,
    const std::vector<std::vector<double>>& sse_prior,
    const std::vector<std::vector<int>>& models,
    const double nub,
    const double nue,
    const MtBedExecutionSpec& execution) {
  validate_mt_bed_problem(
    data, initial, sets, ssb_prior, sse_prior, models,
    execution, nub, nue);
  const std::size_t n = data.genotype.sample_count;
  const std::size_t m = data.genotype.marker_count;
  const std::size_t nt = data.phenotype.n_cols;
  const std::size_t nmodels = models.size();
  const int total_iterations = execution.nit + execution.nburn;

  MtBedCoreResult output;
  MtDefaultCoreResult& result = output.posterior;
  result.nt = nt;
  result.m = m;
  result.nmodels = nmodels;
  result.bm.assign(nt, std::vector<double>(m, 0.0));
  result.dm.assign(nt, std::vector<double>(m, 0.0));
  result.r.assign(nt, std::vector<double>(m, 0.0));
  result.b = initial.b;
  result.d = initial.state;
  result.order = data.marker_order;
  result.vbs.assign(nt, std::vector<double>(total_iterations, 0.0));
  result.vgs.assign(nt, std::vector<double>(total_iterations, 0.0));
  result.ves.assign(nt, std::vector<double>(total_iterations, 0.0));
  result.vle.assign(nt, std::vector<double>(total_iterations, 0.0));
  result.vld.assign(nt, std::vector<double>(total_iterations, 0.0));
  result.pis.assign(nmodels, 0.0);
  result.cvbm.assign(nt, std::vector<double>(nt, 0.0));
  result.cvgm.assign(nt, std::vector<double>(nt, 0.0));
  result.cvem.assign(nt, std::vector<double>(nt, 0.0));
  result.B = initial.B;
  result.E = initial.E;
  result.pi = initial.pi;
  result.pistrait.assign(nt, std::vector<double>(4, 0.0));
  result.pismarker.assign(2, 0.0);

  std::vector<std::vector<double>> beta = initial.beta;
  arma::mat residual = data.phenotype;
  arma::vec marker_workspace(n, arma::fill::zeros);
  for (std::size_t marker = 0; marker < m; ++marker) {
    decode_mt_bed_marker(
      data.genotype, data.marker_maps[marker], marker, marker_workspace);
    for (std::size_t trait = 0; trait < nt; ++trait) {
      if (result.b[trait][marker] != 0.0) {
        residual.col(trait) -= marker_workspace *
          result.b[trait][marker];
      }
    }
  }

  const arma::mat ssb =
    mt_bed_prior_matrix(ssb_prior, nt, "ssb_prior");
  const arma::mat sse =
    mt_bed_prior_matrix(sse_prior, nt, "sse_prior");
  std::mt19937 rng(static_cast<std::mt19937::result_type>(execution.seed));
  arma::mat B_inverse = arma::inv(result.B);
  arma::mat E_inverse = arma::inv(result.E);
  arma::mat genetic_covariance =
    mt_bed_genetic_covariance(data.phenotype, residual);

  for (int iteration = 0; iteration < total_iterations; ++iteration) {
    std::vector<double> cmodel(nmodels, 1.0);
    for (const auto& set : sets) {
      if (execution.updateB) {
        sampleBset(
          static_cast<int>(nt), static_cast<int>(m),
          static_cast<int>(nub), result.B, result.d, result.b,
          ssb_prior, set, rng);
        sampleB_latent(
          static_cast<int>(nt), static_cast<int>(m),
          static_cast<int>(nub), result.B, beta, ssb_prior, rng);
      }
      B_inverse = arma::inv(result.B);

      for (int marker_index : data.marker_order) {
        if (std::find(set.begin(), set.end(), marker_index) == set.end()) {
          continue;
        }
        const std::size_t marker = static_cast<std::size_t>(marker_index);
        decode_mt_bed_marker(
          data.genotype, data.marker_maps[marker], marker, marker_workspace);
        arma::vec score(nt, arma::fill::zeros);
        for (std::size_t trait = 0; trait < nt; ++trait) {
          long double dot = 0.0L;
          for (std::size_t sample = 0; sample < n; ++sample) {
            dot += static_cast<long double>(marker_workspace[sample]) *
              static_cast<long double>(residual(sample, trait));
          }
          score[trait] = static_cast<double>(dot) +
            data.marker_maps[marker].xx * result.b[trait][marker];
        }
        MtBedMarkerKernelResult kernel = mt_bed_marker_kernel(
          score, data.marker_maps[marker].xx, B_inverse, E_inverse,
          models, result.pi, &output.diagnostics);
        const std::size_t selected =
          mt_bed_draw_model(kernel.probability, rng);
        cmodel[selected] += 1.0;

        arma::vec standard_normal(nt);
        std::normal_distribution<double> normal(0.0, 1.0);
        for (std::size_t trait = 0; trait < nt; ++trait) {
          standard_normal[trait] = normal(rng);
        }
        const MtBedMarkerKernelModel& selected_model =
          kernel.models[selected];
        arma::vec beta_new = selected_model.mean +
          arma::solve(
            arma::trimatu(selected_model.lower.t()), standard_normal);
        arma::vec delta(nt, arma::fill::zeros);
        for (std::size_t trait = 0; trait < nt; ++trait) {
          const double effective =
            static_cast<double>(models[selected][trait]) * beta_new[trait];
          delta[trait] = effective - result.b[trait][marker];
          beta[trait][marker] = beta_new[trait];
          result.b[trait][marker] = effective;
          result.d[trait][marker] = models[selected][trait];
        }
        for (std::size_t trait = 0; trait < nt; ++trait) {
          if (delta[trait] != 0.0) {
            residual.col(trait) -= marker_workspace * delta[trait];
          }
        }
      }

    }

    if (iteration >= execution.nburn &&
        (iteration - execution.nburn) % execution.nthin == 0) {
      for (std::size_t marker = 0; marker < m; ++marker) {
        for (std::size_t trait = 0; trait < nt; ++trait) {
          if (result.d[trait][marker] > 0) {
            result.bm[trait][marker] += result.b[trait][marker];
            result.dm[trait][marker] += 1.0;
          }
        }
      }
      ++result.marker_retained_count;
    }

    if (execution.updatePi) {
      samplePi(cmodel, result.pi, rng);
    }
    if (execution.updateB) {
      sampleB(
        static_cast<int>(nt), static_cast<int>(m),
        static_cast<int>(nub), result.B, result.d, result.b,
        ssb_prior, rng);
      for (std::size_t trait = 0; trait < nt; ++trait) {
        result.vbs[trait][iteration] = result.B(trait, trait);
      }
    }

    genetic_covariance =
      mt_bed_genetic_covariance(data.phenotype, residual);
    for (std::size_t trait = 0; trait < nt; ++trait) {
      result.vgs[trait][iteration] = genetic_covariance(trait, trait);
      long double diagonal_contribution = 0.0L;
      for (std::size_t marker = 0; marker < m; ++marker) {
        const long double effect =
          static_cast<long double>(result.b[trait][marker]);
        diagonal_contribution +=
          static_cast<long double>(data.marker_maps[marker].xx) *
          effect * effect;
      }
      result.vle[trait][iteration] = static_cast<double>(
        diagonal_contribution / static_cast<long double>(n));
      result.vld[trait][iteration] =
        result.vgs[trait][iteration] - result.vle[trait][iteration];
    }

    if (execution.updateE) {
      if (execution.residual_covariance == "diagonal") {
        result.E.zeros(nt, nt);
        for (std::size_t trait = 0; trait < nt; ++trait) {
          long double sse_trait = 0.0L;
          for (std::size_t sample = 0; sample < n; ++sample) {
            const double value = residual(sample, trait);
            sse_trait += static_cast<long double>(value) *
              static_cast<long double>(value);
          }
          const double scale =
            static_cast<double>(sse_trait) + nue * sse(trait, trait);
          std::chi_squared_distribution<double> chisq(nue + n);
          const double draw = std::max(chisq(rng), 1e-12);
          result.E(trait, trait) = scale / draw;
        }
        ++output.diagnostics.diagonal_e_updates;
      } else {
        arma::mat crossproduct = residual.t() * residual;
        crossproduct = 0.5 * (crossproduct + crossproduct.t());
        const arma::mat posterior_scale = nue * sse + crossproduct;
        result.E = draw_inverse_wishart(
          nue + static_cast<double>(n), posterior_scale, rng);
        ++output.diagnostics.full_e_updates;
      }
      E_inverse = arma::inv(result.E);
      for (std::size_t trait = 0; trait < nt; ++trait) {
        result.ves[trait][iteration] = result.E(trait, trait);
      }
    }

    if (iteration >= execution.nburn) {
      if (execution.updateB) {
        for (std::size_t row = 0; row < nt; ++row) {
          for (std::size_t col = 0; col < nt; ++col) {
            result.cvbm[row][col] += result.B(row, col);
          }
        }
        ++result.covb_retained_count;
      }
      for (std::size_t row = 0; row < nt; ++row) {
        for (std::size_t col = 0; col < nt; ++col) {
          result.cvgm[row][col] += genetic_covariance(row, col);
        }
      }
      ++result.covg_retained_count;
      if (execution.updateE) {
        for (std::size_t row = 0; row < nt; ++row) {
          for (std::size_t col = 0; col < nt; ++col) {
            result.cvem[row][col] += result.E(row, col);
          }
        }
        ++result.cove_retained_count;
      }
      if (execution.updatePi) {
        for (std::size_t k = 0; k < nmodels; ++k) {
          result.pis[k] += result.pi[k];
        }
        ++result.pi_retained_count;
      }
    }
  }

  if (result.marker_retained_count <= 0.0) {
    throw std::runtime_error("mtblr_bed: no retained marker-summary samples");
  }
  result.G = genetic_covariance;
  for (std::size_t marker = 0; marker < m; ++marker) {
    decode_mt_bed_marker(
      data.genotype, data.marker_maps[marker], marker, marker_workspace);
    for (std::size_t trait = 0; trait < nt; ++trait) {
      long double dot = 0.0L;
      for (std::size_t sample = 0; sample < n; ++sample) {
        dot += static_cast<long double>(marker_workspace[sample]) *
          static_cast<long double>(residual(sample, trait));
      }
      result.r[trait][marker] = static_cast<double>(dot);
    }
  }
  return output;
}

}  // namespace mt
}  // namespace sblr

#endif
