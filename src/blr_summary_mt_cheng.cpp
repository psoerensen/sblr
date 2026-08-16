// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "blr_phase3_execution.h"
#include "blr_summary_mt_cheng.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using sblr::summary_mt::OperatorKind;
using sblr::summary_mt::OperatorResource;
using sblr::summary_mt::Provider;

std::size_t exact_size(double value, const char* field) {
  constexpr double exact_limit = 9007199254740992.0;
  if (!std::isfinite(value) || value < 0.0 || value != std::floor(value) ||
      value > exact_limit || static_cast<long double>(value) >
        static_cast<long double>(std::numeric_limits<std::size_t>::max())) {
    throw std::invalid_argument(std::string(field) +
      " must be an exactly representable nonnegative native size.");
  }
  return static_cast<std::size_t>(value);
}

std::size_t checked_add(std::size_t left, std::size_t right,
                        const char* component) {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    throw std::overflow_error(std::string(
      "Summary Cheng native allocation overflow in ") + component + ".");
  }
  return left + right;
}

std::size_t checked_multiply(std::size_t left, std::size_t right,
                             const char* component) {
  if (left != 0u && right > std::numeric_limits<std::size_t>::max() / left) {
    throw std::overflow_error(std::string(
      "Summary Cheng native allocation overflow in ") + component + ".");
  }
  return left * right;
}

sblr::phase4a::ActivityPatterns parse_patterns(
    Rcpp::IntegerMatrix value, std::size_t trait_count) {
  sblr::phase4a::ActivityPatterns out(
    static_cast<std::size_t>(value.nrow()),
    sblr::phase4a::ActivityPattern(trait_count, 0));
  if (value.ncol() != static_cast<int>(trait_count)) {
    throw std::invalid_argument(
      "Summary Cheng activity-pattern columns must match T.");
  }
  for (int state = 0; state < value.nrow(); ++state) {
    for (int trait = 0; trait < value.ncol(); ++trait) {
      const int bit = value(state, trait);
      if (bit == NA_INTEGER) {
        throw std::invalid_argument(
          "Summary Cheng activity patterns cannot be missing.");
      }
      out[static_cast<std::size_t>(state)][static_cast<std::size_t>(trait)] = bit;
    }
  }
  sblr::phase4a::validate_activity_patterns(out, trait_count);
  return out;
}

void guard_native_dimensions(
    Rcpp::List resources, Rcpp::List providers,
    std::size_t marker_count, std::size_t trait_count,
    std::size_t pattern_count, std::size_t chains,
    std::size_t retained, std::size_t convergence,
    double native_memory_limit_bytes) {
  if (std::isnan(native_memory_limit_bytes) || native_memory_limit_bytes < 0.0 ||
      (std::isinf(native_memory_limit_bytes) && native_memory_limit_bytes < 0.0)) {
    throw std::invalid_argument(
      "Summary Cheng native memory limit must be nonnegative or positive Inf.");
  }
  std::size_t bytes = 0u;
  const auto add_count = [&](std::size_t count, std::size_t width,
                             const char* component) {
    bytes = checked_add(bytes, checked_multiply(count, width, component),
                        component);
  };
  add_count(checked_multiply(marker_count, trait_count,
    "marker-trait state"), checked_multiply(chains, 16u, "chain effects"),
    "chain effect state");
  add_count(checked_multiply(marker_count, pattern_count,
    "marker-pattern output"), 16u, "two marker-pattern matrices");
  add_count(checked_multiply(checked_multiply(chains, retained,
    "chain retained"), marker_count, "retained markers"),
    checked_add(checked_multiply(trait_count, 16u, "retained effects"), 4u,
      "retained state width"), "retained output");
  add_count(checked_multiply(chains, convergence, "chain convergence"),
    checked_add(checked_multiply(trait_count, trait_count * 8u,
      "convergence covariance"), pattern_count * 8u,
      "convergence probability"), "convergence output");
  for (R_xlen_t index = 0; index < resources.size(); ++index) {
    Rcpp::List resource = resources[index];
    const std::string type = Rcpp::as<std::string>(resource["operator_type"]);
    if (type == "csr") {
      Rcpp::NumericVector row_ptr = resource["row_ptr"];
      Rcpp::IntegerVector column = resource["column_index"];
      Rcpp::NumericVector values = resource["values"];
      add_count(static_cast<std::size_t>(row_ptr.size()), 8u,
                "CSR row pointers");
      add_count(static_cast<std::size_t>(column.size()), 8u,
                "CSR columns and float values");
      if (column.size() != values.size()) {
        throw std::invalid_argument("CSR column/value lengths disagree.");
      }
    } else if (type == "full_rank_block_eigen" ||
               type == "retained_rank_block_eigen") {
      Rcpp::List blocks = resource["blocks"];
      for (R_xlen_t block_index = 0; block_index < blocks.size(); ++block_index) {
        Rcpp::List block = blocks[block_index];
        Rcpp::NumericMatrix vectors = block["eigenvectors"];
        Rcpp::NumericVector values = block["eigenvalues"];
        add_count(static_cast<std::size_t>(vectors.size()), 8u,
                  "block-eigen factors");
        add_count(static_cast<std::size_t>(values.size()), 8u,
                  "block-eigen values");
      }
    } else {
      throw std::invalid_argument(
        "Phase 6A supports only CSR and block-eigen resources.");
    }
  }
  for (R_xlen_t index = 0; index < providers.size(); ++index) {
    Rcpp::List provider = providers[index];
    Rcpp::NumericVector score = provider["score"];
    add_count(static_cast<std::size_t>(score.size()),
      checked_add(12u, checked_multiply(chains, 8u,
        "provider chain residual width"), "provider marker width"),
      "provider score, map, and residual state");
  }
  if (std::isfinite(native_memory_limit_bytes) &&
      static_cast<long double>(bytes) >
        static_cast<long double>(native_memory_limit_bytes)) {
    throw std::runtime_error(
      "Summary Cheng native allocation guard failed before operator construction: "
      "estimated_native_bytes=" + std::to_string(bytes) +
      ", native_memory_limit_bytes=" +
      std::to_string(native_memory_limit_bytes) + ".");
  }
}

OperatorResource parse_resource(Rcpp::List value) {
  OperatorResource out;
  const std::string type = Rcpp::as<std::string>(value["operator_type"]);
  Rcpp::NumericVector diagonal = value["diagonal"];
  out.marker_count = static_cast<std::size_t>(diagonal.size());
  out.diagonal = Rcpp::as<arma::rowvec>(diagonal);
  if (type == "csr") {
    out.kind = OperatorKind::csr;
    Rcpp::NumericVector row_ptr = value["row_ptr"];
    Rcpp::IntegerVector column = value["column_index"];
    Rcpp::NumericVector entries = value["values"];
    if (row_ptr.size() != static_cast<R_xlen_t>(out.marker_count + 1u) ||
        column.size() != entries.size()) {
      throw std::invalid_argument("CSR summary resource dimensions disagree.");
    }
    out.csr.ptr.assign(out.marker_count + 1u, 0u);
    for (std::size_t row = 0u; row < out.marker_count; ++row) {
      const std::size_t first = exact_size(row_ptr[static_cast<R_xlen_t>(row)],
                                           "CSR row pointer");
      const std::size_t last = exact_size(row_ptr[static_cast<R_xlen_t>(row + 1u)],
                                          "CSR row pointer");
      if (first > last || last > static_cast<std::size_t>(column.size())) {
        throw std::invalid_argument("CSR row pointers are invalid.");
      }
      bool found_diagonal = false;
      for (std::size_t position = first; position < last; ++position) {
        const int one_based = column[static_cast<R_xlen_t>(position)];
        const double entry = entries[static_cast<R_xlen_t>(position)];
        if (one_based == NA_INTEGER || one_based <= 0 ||
            static_cast<std::size_t>(one_based) > out.marker_count ||
            !std::isfinite(entry)) {
          throw std::invalid_argument("CSR summary entries are invalid.");
        }
        const int local = one_based - 1;
        if (static_cast<std::size_t>(local) == row) {
          if (found_diagonal || std::abs(entry - out.diagonal(row)) > 1e-10) {
            throw std::invalid_argument(
              "CSR explicit diagonal disagrees with its declared diagonal.");
          }
          found_diagonal = true;
        } else {
          out.csr.idx.push_back(local);
          out.csr.xij.push_back(static_cast<float>(entry));
        }
      }
      if (!found_diagonal) {
        throw std::invalid_argument("CSR summary rows require explicit diagonals.");
      }
      out.csr.ptr[row + 1u] = static_cast<std::uint64_t>(out.csr.idx.size());
    }
  } else if (type == "full_rank_block_eigen" ||
             type == "retained_rank_block_eigen") {
    out.kind = OperatorKind::block_eigen;
    Rcpp::List blocks = value["blocks"];
    out.block_of.assign(out.marker_count, -1);
    out.local_in_block.assign(out.marker_count, -1);
    out.blocks.reserve(static_cast<std::size_t>(blocks.size()));
    for (R_xlen_t group = 0; group < blocks.size(); ++group) {
      Rcpp::List current = blocks[group];
      Rcpp::IntegerVector marker = current["marker_indices"];
      Rcpp::NumericMatrix eigenvectors = current["eigenvectors"];
      Rcpp::NumericVector eigenvalues = current["eigenvalues"];
      if (eigenvectors.nrow() != marker.size() ||
          eigenvectors.ncol() != eigenvalues.size() || marker.size() == 0) {
        throw std::invalid_argument(
          "Block-eigen summary resource dimensions disagree.");
      }
      sblr::summary_mt::BlockFactor block;
      block.marker.resize(static_cast<std::size_t>(marker.size()));
      block.factor = Rcpp::as<arma::mat>(eigenvectors);
      for (R_xlen_t component = 0; component < eigenvalues.size(); ++component) {
        const double lambda = eigenvalues[component];
        if (!std::isfinite(lambda) || lambda < 0.0) {
          throw std::invalid_argument(
            "Block-eigen summary eigenvalues must be finite and nonnegative.");
        }
        block.factor.col(static_cast<arma::uword>(component)) *=
          std::sqrt(lambda);
      }
      for (R_xlen_t local = 0; local < marker.size(); ++local) {
        const int one_based = marker[local];
        if (one_based == NA_INTEGER || one_based <= 0 ||
            static_cast<std::size_t>(one_based) > out.marker_count) {
          throw std::invalid_argument(
            "Block-eigen summary marker indices are invalid.");
        }
        const std::size_t local_marker = static_cast<std::size_t>(one_based - 1);
        block.marker[static_cast<std::size_t>(local)] = one_based - 1;
        if (out.block_of[local_marker] != -1) {
          throw std::invalid_argument(
            "Block-eigen summary blocks duplicate a marker.");
        }
        out.block_of[local_marker] = static_cast<int>(group);
        out.local_in_block[local_marker] = static_cast<int>(local);
        const double reconstructed = arma::dot(
          block.factor.row(static_cast<arma::uword>(local)),
          block.factor.row(static_cast<arma::uword>(local)));
        if (std::abs(reconstructed - out.diagonal(local_marker)) >
            1e-9 * std::max(1.0, std::abs(out.diagonal(local_marker)))) {
          throw std::invalid_argument(
            "Block-eigen reconstructed diagonal disagrees with metadata.");
        }
      }
      out.blocks.push_back(std::move(block));
    }
  } else {
    throw std::invalid_argument(
      "Phase 6A supports only CSR and block-eigen resources.");
  }
  sblr::summary_mt::validate_operator_resource(out);
  return out;
}

Provider parse_provider(Rcpp::List value) {
  Provider out;
  out.id = Rcpp::as<std::string>(value["provider_id"]);
  out.trait = exact_size(Rcpp::as<double>(value["trait_index"]),
                         "provider trait index");
  out.resource = exact_size(Rcpp::as<double>(value["resource_index"]),
                            "provider resource index");
  Rcpp::IntegerVector map = value["local_to_global"];
  Rcpp::NumericVector score = value["score"];
  out.local_to_global.resize(static_cast<std::size_t>(map.size()));
  out.score.resize(static_cast<std::size_t>(score.size()));
  for (R_xlen_t index = 0; index < map.size(); ++index) {
    const int current = map[index];
    if (current == NA_INTEGER || current < 0) {
      throw std::invalid_argument("Provider marker maps must be zero based.");
    }
    out.local_to_global[static_cast<std::size_t>(index)] =
      static_cast<std::size_t>(current);
  }
  for (R_xlen_t index = 0; index < score.size(); ++index) {
    out.score[static_cast<std::size_t>(index)] = score[index];
  }
  out.residual_scale = Rcpp::as<double>(value["residual_scale"]);
  out.sample_size = Rcpp::as<double>(value["sample_size"]);
  return out;
}

Rcpp::NumericVector effect_draws(
    const std::vector<arma::mat>& draws, int markers, int traits) {
  const int ndraw = static_cast<int>(draws.size());
  Rcpp::NumericVector out(static_cast<R_xlen_t>(ndraw) * markers * traits);
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

Rcpp::IntegerMatrix state_draws(
    const std::vector<std::vector<int>>& draws, int markers) {
  Rcpp::IntegerMatrix out(static_cast<int>(draws.size()), markers);
  for (int marker = 0; marker < markers; ++marker) {
    for (std::size_t draw = 0u; draw < draws.size(); ++draw) {
      out(static_cast<int>(draw), marker) = draws[draw][static_cast<std::size_t>(marker)];
    }
  }
  return out;
}

Rcpp::NumericVector matrix_draws(
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

Rcpp::NumericMatrix probability_draws(
    const std::vector<sblr::phase4a::ProbabilityVector>& draws,
    int patterns) {
  Rcpp::NumericMatrix out(static_cast<int>(draws.size()), patterns);
  for (std::size_t draw = 0u; draw < draws.size(); ++draw) {
    for (int state = 0; state < patterns; ++state) {
      out(static_cast<int>(draw), state) =
        draws[draw][static_cast<std::size_t>(state)];
    }
  }
  return out;
}

Rcpp::NumericMatrix effect_matrix(const arma::mat& value) {
  Rcpp::NumericMatrix out = Rcpp::wrap(value);
  for (R_xlen_t index = 0; index < out.size(); ++index) {
    if (std::isnan(out[index])) out[index] = NA_REAL;
  }
  return out;
}

Rcpp::List chain_to_list(const sblr::summary_mt::ChainResult& chain) {
  const int markers = static_cast<int>(chain.final_realised.n_rows);
  const int traits = static_cast<int>(chain.final_realised.n_cols);
  const int patterns = static_cast<int>(chain.final_probability.size());
  Rcpp::NumericVector occupancy(patterns);
  for (int state = 0; state < patterns; ++state) {
    occupancy[state] = static_cast<double>(
      chain.pattern_occupancy_count[static_cast<std::size_t>(state)]);
  }
  Rcpp::List provider_residual(
    static_cast<R_xlen_t>(chain.final_provider_residual_score.size()));
  for (std::size_t provider = 0u;
       provider < chain.final_provider_residual_score.size(); ++provider) {
    provider_residual[static_cast<R_xlen_t>(provider)] =
      Rcpp::wrap(chain.final_provider_residual_score[provider]);
  }
  return Rcpp::List::create(
    Rcpp::_ ["realised_effects"] = effect_draws(
      chain.realised_draws, markers, traits),
    Rcpp::_ ["latent_effects"] = effect_draws(
      chain.latent_draws, markers, traits),
    Rcpp::_ ["joint_states"] = state_draws(chain.state_draws, markers),
    Rcpp::_ ["marker_covariance"] = matrix_draws(
      chain.covariance_draws, traits, traits),
    Rcpp::_ ["activity_pattern_parameters"] = probability_draws(
      chain.probability_draws, patterns),
    Rcpp::_ ["convergence_marker_covariance"] = matrix_draws(
      chain.convergence_covariance, traits, traits),
    Rcpp::_ ["convergence_activity_pattern_parameters"] = probability_draws(
      chain.convergence_probability, patterns),
    Rcpp::_ ["convergence_active_marker_count"] =
      Rcpp::wrap(chain.convergence_active_count),
    Rcpp::_ ["final_realised_effects"] = effect_matrix(chain.final_realised),
    Rcpp::_ ["final_latent_effects"] = effect_matrix(chain.final_latent),
    Rcpp::_ ["final_joint_states"] = Rcpp::wrap(chain.final_state),
    Rcpp::_ ["final_marker_covariance"] = Rcpp::wrap(chain.final_covariance),
    Rcpp::_ ["final_activity_pattern_parameters"] =
      Rcpp::wrap(chain.final_probability),
    Rcpp::_ ["final_provider_residual_score"] = provider_residual,
    Rcpp::_ ["last_covariance_statistic"] =
      Rcpp::wrap(chain.last_covariance_statistic),
    Rcpp::_ ["last_covariance_scale"] = Rcpp::wrap(chain.last_covariance_scale),
    Rcpp::_ ["last_covariance_degrees_of_freedom"] = chain.last_covariance_df,
    Rcpp::_ ["last_active_marker_count"] = chain.last_active_count,
    Rcpp::_ ["pattern_occupancy_counts"] = occupancy,
    Rcpp::_ ["pattern_change_count"] =
      static_cast<double>(chain.pattern_change_count));
}

}  // namespace

// Deterministic independent-summary marker conditional used by qualification
// references. It performs no RNG draw.
// [[Rcpp::export]]
Rcpp::List mtblr_phase6a_summary_pattern_contract_internal(
    Rcpp::NumericVector aggregated_score,
    Rcpp::NumericVector aggregated_diagonal,
    arma::mat marker_covariance,
    Rcpp::NumericVector activity_pattern_probability,
    Rcpp::IntegerMatrix activity_patterns) {
  const std::size_t trait_count = static_cast<std::size_t>(aggregated_score.size());
  const auto patterns = parse_patterns(activity_patterns, trait_count);
  if (aggregated_diagonal.size() != aggregated_score.size() ||
      activity_pattern_probability.size() !=
        static_cast<R_xlen_t>(patterns.size())) {
    throw std::invalid_argument(
      "Summary marker-contract vectors have inconsistent dimensions.");
  }
  sblr::phase4a::ProbabilityVector probability(patterns.size(), 0.0);
  double total = 0.0;
  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    const double value = activity_pattern_probability[
      static_cast<R_xlen_t>(state)];
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument(
        "Summary activity-pattern probabilities must be positive.");
    }
    probability[state] = value;
    total += value;
  }
  for (double& value : probability) value /= total;
  const auto kernel = sblr::summary_mt::pattern_kernel(
    Rcpp::as<arma::vec>(aggregated_score),
    Rcpp::as<arma::vec>(aggregated_diagonal),
    sblr::phase4a::require_spd(
      marker_covariance, trait_count, "marker covariance"),
    probability, patterns);
  Rcpp::List means(static_cast<R_xlen_t>(patterns.size()));
  Rcpp::List covariances(static_cast<R_xlen_t>(patterns.size()));
  for (std::size_t state = 1u; state < patterns.size(); ++state) {
    means[static_cast<R_xlen_t>(state)] = Rcpp::wrap(kernel.active_mean[state]);
    covariances[static_cast<R_xlen_t>(state)] =
      Rcpp::wrap(kernel.active_covariance[state]);
  }
  return Rcpp::List::create(
    Rcpp::_ ["probability"] = Rcpp::wrap(kernel.probability),
    Rcpp::_ ["log_weight"] = Rcpp::wrap(kernel.log_weight),
    Rcpp::_ ["active_mean"] = means,
    Rcpp::_ ["active_covariance"] = covariances);
}

// Qualification-only general-T independent-summary Cheng MT-BayesC-Pi route.
// Operator resources are parsed once and shared immutably by all chain tasks.
// [[Rcpp::export]]
Rcpp::List mtblr_phase6a_summary_cheng_internal(
    Rcpp::List operator_resources,
    Rcpp::List likelihood_providers,
    int global_marker_count,
    int trait_count,
    Rcpp::IntegerMatrix activity_patterns,
    arma::mat initial_marker_covariance,
    Rcpp::NumericVector initial_activity_pattern_probability,
    Rcpp::NumericVector activity_pattern_dirichlet_prior,
    double marker_covariance_prior_df,
    arma::mat marker_covariance_prior_scale,
    bool update_marker_covariance,
    bool update_activity_pattern_probability,
    int burn_in_iterations,
    int sampling_iterations,
    int chains,
    int cores,
    Rcpp::Nullable<Rcpp::List> execution_contract,
    double native_memory_limit_bytes) {
  if (global_marker_count <= 0 || trait_count < 2 || chains <= 0 || cores <= 0) {
    throw std::invalid_argument(
      "Summary Cheng dimensions and execution controls must be positive.");
  }
  const std::size_t markers = static_cast<std::size_t>(global_marker_count);
  const std::size_t traits = static_cast<std::size_t>(trait_count);
  const auto patterns = parse_patterns(activity_patterns, traits);
  const BlrPhase3ExecutionContract phase3 =
    parse_blr_phase3_execution_contract(
      execution_contract, chains, sampling_iterations);
  if (!phase3.active()) {
    throw std::invalid_argument(
      "The summary Cheng route requires the active Phase 3 contract.");
  }
  for (int chain = 0; chain < chains; ++chain) {
    if (phase3.task_ids[static_cast<std::size_t>(chain)] !=
        "chain:" + std::to_string(chain)) {
      throw std::invalid_argument(
        "Summary Cheng task IDs must follow canonical chain order.");
    }
  }
  sblr::phase4a::validate_native_allocation_dimensions(
    traits, patterns.size(), markers, 0u, static_cast<std::size_t>(chains),
    phase3.retained_transition_indices.size(),
    static_cast<std::size_t>(sampling_iterations));
  guard_native_dimensions(
    operator_resources, likelihood_providers, markers, traits,
    patterns.size(), static_cast<std::size_t>(chains),
    phase3.retained_transition_indices.size(),
    static_cast<std::size_t>(sampling_iterations), native_memory_limit_bytes);

  std::vector<OperatorResource> resources;
  resources.reserve(static_cast<std::size_t>(operator_resources.size()));
  for (R_xlen_t index = 0; index < operator_resources.size(); ++index) {
    resources.push_back(parse_resource(operator_resources[index]));
  }
  std::vector<Provider> providers;
  providers.reserve(static_cast<std::size_t>(likelihood_providers.size()));
  std::vector<std::vector<sblr::summary_mt::ProviderMarkerReference>>
    marker_providers(markers);
  std::vector<int> trait_seen(traits, 0);
  for (R_xlen_t index = 0; index < likelihood_providers.size(); ++index) {
    providers.push_back(parse_provider(likelihood_providers[index]));
    const std::size_t provider_index = static_cast<std::size_t>(index);
    sblr::summary_mt::validate_provider(
      providers.back(), resources, traits, markers);
    ++trait_seen[providers.back().trait];
    for (std::size_t local = 0u;
         local < providers.back().local_to_global.size(); ++local) {
      marker_providers[providers.back().local_to_global[local]].push_back(
        sblr::summary_mt::ProviderMarkerReference{provider_index, local});
    }
  }
  for (int count : trait_seen) {
    if (count == 0) {
      throw std::invalid_argument(
        "Every summary Cheng trait requires at least one provider.");
    }
  }
  if (initial_activity_pattern_probability.size() !=
        static_cast<R_xlen_t>(patterns.size()) ||
      activity_pattern_dirichlet_prior.size() !=
        static_cast<R_xlen_t>(patterns.size())) {
    throw std::invalid_argument(
      "Summary Cheng probability dimensions disagree with 2^T.");
  }
  sblr::phase4a::ProbabilityVector initial_probability(patterns.size(), 0.0);
  sblr::phase4a::ProbabilityVector dirichlet_prior(patterns.size(), 0.0);
  double probability_total = 0.0;
  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    const double value = initial_activity_pattern_probability[
      static_cast<R_xlen_t>(state)];
    const double prior = activity_pattern_dirichlet_prior[
      static_cast<R_xlen_t>(state)];
    if (!std::isfinite(value) || value <= 0.0 ||
        !std::isfinite(prior) || prior <= 0.0) {
      throw std::invalid_argument(
        "Summary Cheng probabilities and Dirichlet shapes must be positive.");
    }
    initial_probability[state] = value;
    probability_total += value;
    dirichlet_prior[state] = prior;
  }
  for (double& value : initial_probability) value /= probability_total;

  const int configured_workers = blr_phase3_configured_workers(cores, chains);
  std::vector<sblr::summary_mt::ChainResult> results(
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
      results[static_cast<std::size_t>(chain)] = sblr::summary_mt::run_chain(
        resources, providers, marker_providers, markers, traits, patterns,
        initial_marker_covariance, initial_probability, dirichlet_prior,
        marker_covariance_prior_df, marker_covariance_prior_scale,
        update_marker_covariance, update_activity_pattern_probability,
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
        "Summary Cheng chain " + std::to_string(chain + 1) +
        " failed: " + errors[static_cast<std::size_t>(chain)]);
    }
  }
  Rcpp::List chain_results(chains);
  Rcpp::NumericVector native_task_seeds(chains);
  for (int chain = 0; chain < chains; ++chain) {
    chain_results[chain] = chain_to_list(results[static_cast<std::size_t>(chain)]);
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
    Rcpp::_ ["operator_contract"] =
      "shared_summary_provider_operator_v1",
    Rcpp::_ ["qualification_only"] = true);
}
