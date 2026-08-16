// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "blr_phase3_execution.h"
#include "blr_phase4a_cheng_mt.h"
#include "packed_bed.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Phase4aPreparedBed {
  std::unique_ptr<PackedBedMatrix> owner;
  std::vector<sblr::mt::MtBedMarkerMap> marker_maps;
  arma::mat phenotype;

  sblr::core::BedPackedGenotypeView<PackedBedMatrix> view() const {
    return sblr::core::BedPackedGenotypeView<PackedBedMatrix>{
      *owner, owner->data,
      static_cast<std::size_t>(owner->m) * owner->stride,
      static_cast<std::size_t>(owner->m),
      static_cast<std::size_t>(owner->n), owner->nbytes, owner->stride
    };
  }
};

sblr::phase4a::ActivityPatterns phase4a_activity_patterns(
    Rcpp::IntegerMatrix value, std::size_t trait_count) {
  sblr::phase4a::ActivityPatterns out(
    static_cast<std::size_t>(value.nrow()),
    sblr::phase4a::ActivityPattern(trait_count, 0));
  if (value.ncol() != static_cast<int>(trait_count)) {
    throw std::invalid_argument(
      "Activity-pattern columns must match the phenotype traits.");
  }
  for (int state = 0; state < value.nrow(); ++state) {
    for (int trait = 0; trait < value.ncol(); ++trait) {
      const int bit = value(state, trait);
      if (bit == NA_INTEGER) {
        throw std::invalid_argument("Activity patterns cannot be missing.");
      }
      out[static_cast<std::size_t>(state)][static_cast<std::size_t>(trait)] = bit;
    }
  }
  sblr::phase4a::validate_activity_patterns(out, trait_count);
  return out;
}

sblr::phase4a::ActivityPatterns phase4a_canonical_patterns(
    std::size_t trait_count) {
  if (trait_count < 2u ||
      trait_count > sblr::phase4a::maximum_trait_count) {
    throw std::invalid_argument(
      "Complete Cheng activity-pattern enumeration supports T in [2, 12].");
  }
  const std::size_t pattern_count = std::size_t{1u} << trait_count;
  sblr::phase4a::ActivityPatterns out(
    pattern_count, sblr::phase4a::ActivityPattern(trait_count, 0));
  for (std::size_t state = 0u; state < pattern_count; ++state) {
    for (std::size_t trait = 0u; trait < trait_count; ++trait) {
      out[state][trait] = static_cast<int>((state >> trait) & 1u);
    }
  }
  return out;
}

Rcpp::CharacterVector phase4a_pattern_ids(
    const sblr::phase4a::ActivityPatterns& patterns) {
  Rcpp::CharacterVector out(static_cast<R_xlen_t>(patterns.size()));
  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    std::string id;
    for (std::size_t trait = 0u; trait < patterns[state].size(); ++trait) {
      if (trait) id += "_";
      id += patterns[state][trait] ? "1" : "0";
    }
    out[static_cast<R_xlen_t>(state)] = id;
  }
  return out;
}

std::size_t phase5a_exact_size(double value, const char* component) {
  constexpr double exact_integer_limit = 9007199254740992.0;
  if (!std::isfinite(value) || value < 0.0 || value != std::floor(value) ||
      value > exact_integer_limit ||
      static_cast<long double>(value) >
        static_cast<long double>(std::numeric_limits<std::size_t>::max())) {
    throw std::invalid_argument(
      std::string("Packed BED allocation component '") + component +
      "' is not an exactly representable native size.");
  }
  return static_cast<std::size_t>(value);
}

void phase5a_guard_packed_bed_allocation(
    std::size_t selected_sample_count,
    std::size_t source_sample_count,
    std::size_t marker_count,
    bool selected_rows_used,
    double memory_limit_bytes) {
  if (source_sample_count < selected_sample_count) {
    throw std::invalid_argument(
      "Packed BED source sample count is smaller than the selected count.");
  }
  if (!selected_rows_used && source_sample_count != selected_sample_count) {
    throw std::invalid_argument(
      "Packed BED all-sample preparation requires equal source and selected counts.");
  }
  if (std::isnan(memory_limit_bytes) || memory_limit_bytes < 0.0 ||
      (std::isinf(memory_limit_bytes) && memory_limit_bytes < 0.0)) {
    throw std::invalid_argument(
      "The native packed-BED memory limit must be nonnegative or positive Inf.");
  }
  const std::size_t owner_bytes = packed_bed_owner_bytes(
    selected_sample_count, marker_count);
  const std::size_t source_row_bytes = selected_rows_used ?
    packed_bed_row_bytes(source_sample_count) : 0u;
  const std::size_t guarded_bytes = packed_bed_checked_add(
    owner_bytes, source_row_bytes, "packed_bed_guard_total");
  if (std::isfinite(memory_limit_bytes) &&
      static_cast<long double>(guarded_bytes) >
        static_cast<long double>(memory_limit_bytes)) {
    throw std::runtime_error(
      "Packed BED allocation guard failed before construction: "
      "packed_bed_owner=" + std::to_string(owner_bytes) +
      ", packed_bed_source_row_buffer=" +
      std::to_string(source_row_bytes) +
      ", native_memory_limit_bytes=" +
      std::to_string(memory_limit_bytes) + ".");
  }
}

Phase4aPreparedBed prepare_phase4a_bed(
    Rcpp::CharacterVector bed_files,
    int source_sample_count,
    Rcpp::List selected_columns,
    Rcpp::Nullable<Rcpp::IntegerVector> selected_rows,
    Rcpp::NumericVector allele_frequency,
    Rcpp::NumericMatrix phenotype) {
  if (source_sample_count <= 1 || bed_files.size() == 0 ||
      selected_columns.size() != bed_files.size()) {
    throw std::invalid_argument(
      "Cheng MT BED files, columns, or source sample count are invalid.");
  }
  std::vector<std::string> paths;
  paths.reserve(static_cast<std::size_t>(bed_files.size()));
  for (R_xlen_t file = 0; file < bed_files.size(); ++file) {
    if (bed_files[file] == NA_STRING) {
      throw std::invalid_argument("Cheng MT BED paths cannot be missing.");
    }
    const std::string path = Rcpp::as<std::string>(bed_files[file]);
    if (path.empty()) {
      throw std::invalid_argument("Cheng MT BED paths cannot be empty.");
    }
    paths.push_back(path);
  }

  std::vector<std::vector<int>> columns(
    static_cast<std::size_t>(selected_columns.size()));
  std::size_t marker_count = 0u;
  for (R_xlen_t file = 0; file < selected_columns.size(); ++file) {
    Rcpp::IntegerVector current = selected_columns[file];
    if (current.size() == 0) {
      throw std::invalid_argument(
        "Every Cheng MT BED file must contribute a selected marker.");
    }
    std::set<int> seen;
    for (int column : current) {
      if (column == NA_INTEGER || column <= 0 || !seen.insert(column).second) {
        throw std::invalid_argument(
          "Cheng MT BED columns must be unique positive one-based indices.");
      }
      columns[static_cast<std::size_t>(file)].push_back(column);
      ++marker_count;
    }
  }
  if (allele_frequency.size() != static_cast<R_xlen_t>(marker_count)) {
    throw std::invalid_argument(
      "Cheng MT allele frequencies must match selected markers.");
  }
  std::vector<double> frequencies(marker_count);
  for (std::size_t marker = 0u; marker < marker_count; ++marker) {
    const double value = allele_frequency[static_cast<R_xlen_t>(marker)];
    if (!std::isfinite(value) || value <= 0.0 || value >= 1.0) {
      throw std::invalid_argument(
        "Cheng MT allele frequencies must lie strictly inside (0, 1).");
    }
    frequencies[marker] = value;
  }

  std::vector<int> rows0;
  const int* row_pointer = nullptr;
  int selected_sample_count = source_sample_count;
  if (selected_rows.isNotNull()) {
    Rcpp::IntegerVector rows(selected_rows);
    if (rows.size() <= 1) {
      throw std::invalid_argument(
        "Cheng MT selected rows must contain at least two samples.");
    }
    std::set<int> seen;
    rows0.reserve(static_cast<std::size_t>(rows.size()));
    for (int row : rows) {
      if (row == NA_INTEGER || row <= 0 || row > source_sample_count ||
          !seen.insert(row).second) {
        throw std::invalid_argument(
          "Cheng MT selected rows are invalid or duplicated.");
      }
      rows0.push_back(row - 1);
    }
    row_pointer = rows0.data();
    selected_sample_count = static_cast<int>(rows0.size());
  }
  if (phenotype.nrow() != selected_sample_count || phenotype.ncol() < 2) {
    throw std::invalid_argument(
      "Cheng MT phenotype must have selected samples by T >= 2 traits.");
  }
  for (double value : phenotype) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("Cheng MT phenotype must be complete and finite.");
    }
  }

  Phase4aPreparedBed out;
  out.owner.reset(new PackedBedMatrix(read_bedfiles_to_packed_matrix(
    paths, source_sample_count, row_pointer, static_cast<int>(rows0.size()),
    columns)));
  const auto genotype = out.view();
  out.marker_maps = sblr::mt::build_mt_bed_marker_maps(genotype, frequencies);
  out.phenotype = Rcpp::as<arma::mat>(phenotype);
  return out;
}

Rcpp::NumericVector phase4a_effect_draws(
    const std::vector<arma::mat>& draws, int markers, int traits) {
  const int ndraw = static_cast<int>(draws.size());
  Rcpp::NumericVector out(
    static_cast<R_xlen_t>(ndraw) * markers * traits);
  out.attr("dim") = Rcpp::IntegerVector::create(ndraw, markers, traits);
  for (int trait = 0; trait < traits; ++trait) {
    for (int marker = 0; marker < markers; ++marker) {
      for (int draw = 0; draw < ndraw; ++draw) {
        const double value = draws[static_cast<std::size_t>(draw)](marker, trait);
        out[draw + ndraw * (marker + markers * trait)] =
          std::isnan(value) ? NA_REAL : value;
      }
    }
  }
  return out;
}

Rcpp::IntegerMatrix phase4a_state_draws(
    const std::vector<std::vector<int>>& draws, int markers) {
  const int ndraw = static_cast<int>(draws.size());
  Rcpp::IntegerMatrix out(ndraw, markers);
  for (int marker = 0; marker < markers; ++marker) {
    for (int draw = 0; draw < ndraw; ++draw) {
      out(draw, marker) = draws[static_cast<std::size_t>(draw)][
        static_cast<std::size_t>(marker)];
    }
  }
  return out;
}

Rcpp::NumericVector phase4a_matrix_draws(
    const std::vector<arma::mat>& draws, int rows, int columns) {
  const int ndraw = static_cast<int>(draws.size());
  Rcpp::NumericVector out(static_cast<R_xlen_t>(ndraw) * rows * columns);
  out.attr("dim") = Rcpp::IntegerVector::create(ndraw, rows, columns);
  for (int column = 0; column < columns; ++column) {
    for (int row = 0; row < rows; ++row) {
      for (int draw = 0; draw < ndraw; ++draw) {
        out[draw + ndraw * (row + rows * column)] =
          draws[static_cast<std::size_t>(draw)](row, column);
      }
    }
  }
  return out;
}

Rcpp::NumericMatrix phase4a_probability_draws(
    const std::vector<sblr::phase4a::ProbabilityVector>& draws,
    int pattern_count) {
  Rcpp::NumericMatrix out(static_cast<int>(draws.size()), pattern_count);
  for (std::size_t draw = 0u; draw < draws.size(); ++draw) {
    for (int state = 0; state < pattern_count; ++state) {
      out(static_cast<int>(draw), state) =
        draws[draw][static_cast<std::size_t>(state)];
    }
  }
  return out;
}

Rcpp::NumericMatrix phase4a_effect_matrix(const arma::mat& value) {
  Rcpp::NumericMatrix out = Rcpp::wrap(value);
  for (R_xlen_t index = 0; index < out.size(); ++index) {
    if (std::isnan(out[index])) out[index] = NA_REAL;
  }
  return out;
}

Rcpp::List phase4a_chain_to_list(const sblr::phase4a::ChainResult& chain) {
  const int markers = static_cast<int>(chain.final_realised.n_rows);
  const int traits = static_cast<int>(chain.final_realised.n_cols);
  const int patterns = static_cast<int>(chain.final_probability.size());
  Rcpp::NumericVector occupancy(patterns);
  for (int state = 0; state < patterns; ++state) {
    occupancy[state] = static_cast<double>(
      chain.pattern_occupancy_count[static_cast<std::size_t>(state)]);
  }
  return Rcpp::List::create(
    Rcpp::_ ["realised_effects"] = phase4a_effect_draws(
      chain.realised_draws, markers, traits),
    Rcpp::_ ["latent_effects"] = phase4a_effect_draws(
      chain.latent_draws, markers, traits),
    Rcpp::_ ["joint_states"] = phase4a_state_draws(
      chain.state_draws, markers),
    Rcpp::_ ["marker_covariance"] = phase4a_matrix_draws(
      chain.covariance_draws, traits, traits),
    Rcpp::_ ["residual_covariance"] = phase4a_matrix_draws(
      chain.residual_covariance_draws, traits, traits),
    Rcpp::_ ["activity_pattern_parameters"] = phase4a_probability_draws(
      chain.probability_draws, patterns),
    Rcpp::_ ["predictions"] = phase4a_matrix_draws(
      chain.prediction_draws, static_cast<int>(chain.final_prediction.n_rows),
      traits),
    Rcpp::_ ["convergence_marker_covariance"] = phase4a_matrix_draws(
      chain.convergence_covariance, traits, traits),
    Rcpp::_ ["convergence_residual_covariance"] = phase4a_matrix_draws(
      chain.convergence_residual_covariance, traits, traits),
    Rcpp::_ ["convergence_activity_pattern_parameters"] =
      phase4a_probability_draws(chain.convergence_probability, patterns),
    Rcpp::_ ["convergence_active_marker_count"] =
      Rcpp::wrap(chain.convergence_active_count),
    Rcpp::_ ["final_realised_effects"] = phase4a_effect_matrix(
      chain.final_realised),
    Rcpp::_ ["final_latent_effects"] = phase4a_effect_matrix(
      chain.final_latent),
    Rcpp::_ ["final_joint_states"] = Rcpp::wrap(chain.final_state),
    Rcpp::_ ["final_marker_covariance"] = Rcpp::wrap(chain.final_covariance),
    Rcpp::_ ["final_residual_covariance"] = Rcpp::wrap(
      chain.final_residual_covariance),
    Rcpp::_ ["final_activity_pattern_parameters"] =
      Rcpp::wrap(chain.final_probability),
    Rcpp::_ ["final_predictions"] = Rcpp::wrap(chain.final_prediction),
    Rcpp::_ ["last_covariance_statistic"] = Rcpp::wrap(
      chain.last_covariance_statistic),
    Rcpp::_ ["last_covariance_scale"] = Rcpp::wrap(
      chain.last_covariance_scale),
    Rcpp::_ ["last_covariance_degrees_of_freedom"] = chain.last_covariance_df,
    Rcpp::_ ["last_active_marker_count"] = chain.last_active_count,
    Rcpp::_ ["last_residual_covariance_statistic"] = Rcpp::wrap(
      chain.last_residual_statistic),
    Rcpp::_ ["last_residual_covariance_scale"] = Rcpp::wrap(
      chain.last_residual_scale),
    Rcpp::_ ["last_residual_covariance_degrees_of_freedom"] =
      chain.last_residual_df,
    Rcpp::_ ["residual_covariance_update_count"] =
      chain.residual_covariance_update_count,
    Rcpp::_ ["pattern_occupancy_counts"] = occupancy,
    Rcpp::_ ["pattern_change_count"] = static_cast<double>(
      chain.pattern_change_count));
}

}  // namespace

// Exact native packed-BED allocation contract used for R/native parity tests.
// It performs no allocation and does not read a BED file.
// [[Rcpp::export]]
Rcpp::List mtblr_phase5a_packed_bed_allocation_internal(
    double selected_sample_count,
    double marker_count,
    double source_sample_count,
    bool selected_rows_used) {
  const std::size_t selected = phase5a_exact_size(
    selected_sample_count, "selected_sample_count");
  const std::size_t markers = phase5a_exact_size(
    marker_count, "marker_count");
  const std::size_t source = phase5a_exact_size(
    source_sample_count, "source_sample_count");
  if (source < selected || (!selected_rows_used && source != selected)) {
    throw std::invalid_argument(
      "Packed BED source/selected sample counts disagree with selection policy.");
  }
  const std::size_t row_bytes = packed_bed_row_bytes(selected);
  const std::size_t stride = packed_bed_aligned_stride(selected);
  const std::size_t owner = packed_bed_owner_bytes(selected, markers);
  const std::size_t source_row = selected_rows_used ?
    packed_bed_row_bytes(source) : 0u;
  return Rcpp::List::create(
    Rcpp::_ ["selected_row_bytes"] = static_cast<double>(row_bytes),
    Rcpp::_ ["aligned_stride_bytes"] = static_cast<double>(stride),
    Rcpp::_ ["owner_bytes"] = static_cast<double>(owner),
    Rcpp::_ ["source_row_buffer_bytes"] = static_cast<double>(source_row),
    Rcpp::_ ["selected_rows_used"] = selected_rows_used);
}

// Deterministic marker-conditional oracle. It performs no draw and does not
// access a genotype resource. The optional pattern matrix permits general-T
// qualification while preserving the original two-trait call shape.
// [[Rcpp::export]]
Rcpp::List mtblr_phase4a_pattern_contract_internal(
    Rcpp::NumericVector score,
    double marker_sum_squares,
    arma::mat marker_covariance,
    arma::mat residual_covariance,
    Rcpp::NumericVector activity_pattern_probability,
    Rcpp::Nullable<Rcpp::IntegerMatrix> activity_patterns = R_NilValue) {
  const std::size_t trait_count = static_cast<std::size_t>(score.size());
  const auto patterns = activity_patterns.isNull() ?
    phase4a_canonical_patterns(trait_count) :
    phase4a_activity_patterns(Rcpp::IntegerMatrix(activity_patterns), trait_count);
  if (activity_pattern_probability.size() !=
      static_cast<R_xlen_t>(patterns.size())) {
    throw std::invalid_argument(
      "Marker inspection requires one probability per activity pattern.");
  }
  const arma::mat Vb = sblr::phase4a::require_spd(
    marker_covariance, trait_count, "marker covariance");
  const arma::mat Ve = sblr::phase4a::require_spd(
    residual_covariance, trait_count, "residual covariance");
  sblr::phase4a::ProbabilityVector probability(patterns.size(), 0.0);
  double total = 0.0;
  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    const double value = activity_pattern_probability[
      static_cast<R_xlen_t>(state)];
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument(
        "Activity-pattern probabilities must be finite and positive.");
    }
    probability[state] = value;
    total += value;
  }
  for (double& value : probability) value /= total;
  const auto kernel = sblr::phase4a::pattern_kernel(
    Rcpp::as<arma::vec>(score), marker_sum_squares, Vb,
    sblr::phase4a::inverse_spd(Ve), probability, patterns);
  Rcpp::List means(static_cast<R_xlen_t>(patterns.size()));
  Rcpp::List covariances(static_cast<R_xlen_t>(patterns.size()));
  Rcpp::List completion_coefficients(static_cast<R_xlen_t>(patterns.size()));
  Rcpp::List completion_covariances(static_cast<R_xlen_t>(patterns.size()));
  for (std::size_t state = 1u; state < patterns.size(); ++state) {
    means[static_cast<R_xlen_t>(state)] = Rcpp::wrap(kernel.active_mean[state]);
    covariances[static_cast<R_xlen_t>(state)] =
      Rcpp::wrap(kernel.active_covariance[state]);
    const arma::uvec active = sblr::phase4a::active_indices(patterns[state]);
    const arma::uvec inactive = sblr::phase4a::inactive_indices(patterns[state]);
    if (inactive.n_elem) {
      const arma::mat aa = Vb.submat(active, active);
      const arma::mat ia = Vb.submat(inactive, active);
      const arma::mat ai = Vb.submat(active, inactive);
      completion_coefficients[static_cast<R_xlen_t>(state)] =
        Rcpp::wrap(ia * sblr::phase4a::inverse_spd(aa));
      completion_covariances[static_cast<R_xlen_t>(state)] = Rcpp::wrap(
        0.5 * ((Vb.submat(inactive, inactive) -
          ia * sblr::phase4a::solve_spd(aa, ai)) +
          (Vb.submat(inactive, inactive) -
          ia * sblr::phase4a::solve_spd(aa, ai)).t()));
    }
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_ ["pattern_order"] = phase4a_pattern_ids(patterns),
    Rcpp::_ ["probability"] = Rcpp::wrap(kernel.probability),
    Rcpp::_ ["log_weight"] = Rcpp::wrap(kernel.log_weight),
    Rcpp::_ ["active_mean"] = means,
    Rcpp::_ ["active_covariance"] = covariances,
    Rcpp::_ ["completion_mean_coefficient"] = completion_coefficients,
    Rcpp::_ ["completion_covariance"] = completion_covariances);
  if (trait_count == 2u) {
    out["completion_coefficient"] = Rcpp::NumericVector::create(
      Vb(1, 0) / Vb(0, 0), Vb(0, 1) / Vb(1, 1));
    out["completion_variance"] = Rcpp::NumericVector::create(
      Vb(1, 1) - Vb(1, 0) * Vb(0, 1) / Vb(0, 0),
      Vb(0, 0) - Vb(0, 1) * Vb(1, 0) / Vb(1, 1));
  }
  return out;
}

// Qualification-only residual inverse-Wishart conditional.
// [[Rcpp::export]]
Rcpp::List mtblr_phase4b_residual_covariance_contract_internal(
    Rcpp::NumericMatrix residual,
    double prior_df,
    arma::mat prior_scale,
    int draws,
    double seed) {
  const std::size_t trait_count = static_cast<std::size_t>(residual.ncol());
  if (residual.nrow() < 1 || trait_count < 2u ||
      trait_count > sblr::phase4a::maximum_trait_count || draws < 1 ||
      !std::isfinite(prior_df) || prior_df <= trait_count - 1u ||
      !std::isfinite(seed) || seed < 0.0 || seed > 4294967295.0 ||
      seed != std::floor(seed)) {
    throw std::invalid_argument(
      "Residual inverse-Wishart conditional inputs are invalid.");
  }
  for (double value : residual) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("Residuals must be finite.");
    }
  }
  const arma::mat scale = sblr::phase4a::require_spd(
    prior_scale, trait_count, "residual-covariance prior scale");
  const arma::mat residual_matrix = Rcpp::as<arma::mat>(residual);
  arma::mat statistic = residual_matrix.t() * residual_matrix;
  statistic = 0.5 * (statistic + statistic.t());
  const double posterior_df = prior_df + residual.nrow();
  const arma::mat posterior_scale = scale + statistic;
  std::mt19937 rng(static_cast<std::uint32_t>(seed));
  std::vector<arma::mat> covariance_draws(static_cast<std::size_t>(draws));
  for (int draw = 0; draw < draws; ++draw) {
    covariance_draws[static_cast<std::size_t>(draw)] =
      sblr::phase4a::require_spd(
        sblr::mt::draw_inverse_wishart(posterior_df, posterior_scale, rng),
        trait_count, "sampled residual covariance");
  }
  return Rcpp::List::create(
    Rcpp::_ ["degrees_of_freedom"] = posterior_df,
    Rcpp::_ ["scale"] = posterior_scale,
    Rcpp::_ ["statistic"] = statistic,
    Rcpp::_ ["draws"] = phase4a_matrix_draws(
      covariance_draws, static_cast<int>(trait_count),
      static_cast<int>(trait_count)));
}

// Qualification-only general-T common-sample Cheng MT-BayesC-Pi route.
// It is deliberately not namespace-exported and does not call legacy MT code.
// [[Rcpp::export]]
Rcpp::List mtblr_phase4a_cheng_bed_internal(
    Rcpp::CharacterVector bed_files,
    int source_sample_count,
    Rcpp::List selected_columns,
    Rcpp::Nullable<Rcpp::IntegerVector> selected_rows,
    Rcpp::NumericVector allele_frequency,
    Rcpp::NumericMatrix phenotype,
    Rcpp::IntegerMatrix activity_patterns,
    arma::mat initial_residual_covariance,
    arma::mat initial_marker_covariance,
    Rcpp::NumericVector initial_activity_pattern_probability,
    Rcpp::NumericVector activity_pattern_dirichlet_prior,
    double marker_covariance_prior_df,
    arma::mat marker_covariance_prior_scale,
    bool update_marker_covariance,
    bool update_activity_pattern_probability,
    bool update_residual_covariance,
    double residual_covariance_prior_df,
    arma::mat residual_covariance_prior_scale,
    int burn_in_iterations,
    int sampling_iterations,
    int chains,
    int cores,
    Rcpp::Nullable<Rcpp::List> execution_contract,
    double native_memory_limit_bytes) {
  const std::size_t trait_count = static_cast<std::size_t>(phenotype.ncol());
  const auto patterns = phase4a_activity_patterns(activity_patterns, trait_count);
  const std::size_t pattern_count = patterns.size();
  if (chains <= 0 || cores <= 0 ||
      initial_activity_pattern_probability.size() !=
        static_cast<R_xlen_t>(pattern_count) ||
      activity_pattern_dirichlet_prior.size() !=
        static_cast<R_xlen_t>(pattern_count)) {
    throw std::invalid_argument(
      "Cheng MT requires complete patterns and positive chain controls.");
  }
  sblr::phase4a::ProbabilityVector initial_probability(pattern_count, 0.0);
  sblr::phase4a::ProbabilityVector dirichlet_prior(pattern_count, 0.0);
  double probability_total = 0.0;
  for (std::size_t state = 0u; state < pattern_count; ++state) {
    const double value = initial_activity_pattern_probability[
      static_cast<R_xlen_t>(state)];
    const double prior = activity_pattern_dirichlet_prior[
      static_cast<R_xlen_t>(state)];
    if (!std::isfinite(value) || value <= 0.0 ||
        !std::isfinite(prior) || prior <= 0.0) {
      throw std::invalid_argument(
        "Probabilities and Dirichlet shapes must be positive.");
    }
    initial_probability[state] = value;
    probability_total += value;
    dirichlet_prior[state] = prior;
  }
  for (double& value : initial_probability) value /= probability_total;

  const BlrPhase3ExecutionContract phase3 =
    parse_blr_phase3_execution_contract(
      execution_contract, chains, sampling_iterations);
  if (!phase3.active()) {
    throw std::invalid_argument(
      "The Cheng qualification route requires the active Phase 3 execution contract.");
  }
  for (int chain = 0; chain < chains; ++chain) {
    const std::string expected = "chain:" + std::to_string(chain);
    if (phase3.task_ids[static_cast<std::size_t>(chain)] != expected) {
      throw std::invalid_argument(
        "Task IDs must follow canonical joint-multitrait chain order.");
    }
  }
  sblr::phase4a::validate_native_allocation_dimensions(
    trait_count, pattern_count,
    static_cast<std::size_t>(allele_frequency.size()),
    static_cast<std::size_t>(phenotype.nrow()),
    static_cast<std::size_t>(chains),
    phase3.retained_transition_indices.size(),
    static_cast<std::size_t>(sampling_iterations));

  if (source_sample_count <= 0 || phenotype.nrow() <= 0) {
    throw std::invalid_argument(
      "Packed BED source and selected sample counts must be positive.");
  }
  phase5a_guard_packed_bed_allocation(
    static_cast<std::size_t>(phenotype.nrow()),
    static_cast<std::size_t>(source_sample_count),
    static_cast<std::size_t>(allele_frequency.size()),
    selected_rows.isNotNull(), native_memory_limit_bytes);

  Phase4aPreparedBed prepared = prepare_phase4a_bed(
    bed_files, source_sample_count, selected_columns, selected_rows,
    allele_frequency, phenotype);
  const auto genotype = prepared.view();

  const int configured_workers = blr_phase3_configured_workers(cores, chains);
  std::vector<sblr::phase4a::ChainResult> results(
    static_cast<std::size_t>(chains));
  std::vector<int> worker_ids(static_cast<std::size_t>(chains), 0);
  std::vector<int> team_sizes(static_cast<std::size_t>(chains), 1);
  std::vector<std::string> errors(static_cast<std::size_t>(chains));
  std::vector<int> failed(static_cast<std::size_t>(chains), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(configured_workers)
#endif
  for (int chain = 0; chain < chains; ++chain) {
#ifdef _OPENMP
    worker_ids[static_cast<std::size_t>(chain)] = omp_get_thread_num();
    team_sizes[static_cast<std::size_t>(chain)] = omp_get_num_threads();
#endif
    try {
      results[static_cast<std::size_t>(chain)] = sblr::phase4a::run_chain(
        genotype, prepared.marker_maps, prepared.phenotype, patterns,
        initial_residual_covariance, initial_marker_covariance,
        initial_probability, dirichlet_prior, marker_covariance_prior_df,
        marker_covariance_prior_scale, update_marker_covariance,
        update_activity_pattern_probability, update_residual_covariance,
        residual_covariance_prior_df, residual_covariance_prior_scale,
        burn_in_iterations, sampling_iterations,
        phase3.retained_transition_indices,
        phase3.task_seeds[static_cast<std::size_t>(chain)]);
    } catch (const std::exception& error) {
      failed[static_cast<std::size_t>(chain)] = 1;
      errors[static_cast<std::size_t>(chain)] = error.what();
    } catch (...) {
      failed[static_cast<std::size_t>(chain)] = 1;
      errors[static_cast<std::size_t>(chain)] = "unknown native error";
    }
  }
  for (int chain = 0; chain < chains; ++chain) {
    if (failed[static_cast<std::size_t>(chain)]) {
      throw std::runtime_error(
        "Cheng MT chain " + std::to_string(chain + 1) + " failed: " +
        errors[static_cast<std::size_t>(chain)]);
    }
  }

  Rcpp::List chain_results(chains);
  Rcpp::NumericVector native_task_seeds(chains);
  for (int chain = 0; chain < chains; ++chain) {
    chain_results[chain] = phase4a_chain_to_list(
      results[static_cast<std::size_t>(chain)]);
    native_task_seeds[chain] = static_cast<double>(
      phase3.task_seeds[static_cast<std::size_t>(chain)]);
  }
  return Rcpp::List::create(
    Rcpp::_ ["chains"] = chain_results,
    Rcpp::_ ["workers"] = blr_phase3_worker_diagnostics(
      phase3, cores, configured_workers, worker_ids, team_sizes),
    Rcpp::_ ["task_ids"] = phase3.task_ids,
    Rcpp::_ ["task_seeds"] = native_task_seeds,
    Rcpp::_ ["retained_transition_indices"] =
      phase3.retained_transition_indices,
    Rcpp::_ ["convergence_transition_indices"] =
      Rcpp::seq(1, sampling_iterations),
    Rcpp::_ ["genotype_contract"] =
      "PackedBedMatrix/BedPackedGenotypeView",
    Rcpp::_ ["activity_patterns"] = activity_patterns,
    Rcpp::_ ["qualification_only"] = true);
}
