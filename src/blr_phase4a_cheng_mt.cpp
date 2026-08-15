// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "blr_phase3_execution.h"
#include "blr_phase4a_cheng_mt.h"
#include "packed_bed.h"

#include <array>
#include <cmath>
#include <cstdint>
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
      "Phase 4a BED files, columns, or source sample count are invalid.");
  }
  std::vector<std::string> paths;
  paths.reserve(static_cast<std::size_t>(bed_files.size()));
  for (R_xlen_t file = 0; file < bed_files.size(); ++file) {
    if (bed_files[file] == NA_STRING) {
      throw std::invalid_argument("Phase 4a BED paths cannot be missing.");
    }
    const std::string path = Rcpp::as<std::string>(bed_files[file]);
    if (path.empty()) {
      throw std::invalid_argument("Phase 4a BED paths cannot be empty.");
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
        "Every Phase 4a BED file must contribute a selected marker.");
    }
    std::set<int> seen;
    for (int column : current) {
      if (column == NA_INTEGER || column <= 0 || !seen.insert(column).second) {
        throw std::invalid_argument(
          "Phase 4a BED columns must be unique positive one-based indices.");
      }
      columns[static_cast<std::size_t>(file)].push_back(column);
      ++marker_count;
    }
  }
  if (allele_frequency.size() != static_cast<R_xlen_t>(marker_count)) {
    throw std::invalid_argument(
      "Phase 4a allele frequencies must match selected markers.");
  }
  std::vector<double> frequencies(marker_count);
  for (std::size_t marker = 0u; marker < marker_count; ++marker) {
    const double value = allele_frequency[static_cast<R_xlen_t>(marker)];
    if (!std::isfinite(value) || value <= 0.0 || value >= 1.0) {
      throw std::invalid_argument(
        "Phase 4a allele frequencies must lie strictly inside (0, 1).");
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
        "Phase 4a selected rows must contain at least two samples.");
    }
    std::set<int> seen;
    rows0.reserve(static_cast<std::size_t>(rows.size()));
    for (int row : rows) {
      if (row == NA_INTEGER || row <= 0 || row > source_sample_count ||
          !seen.insert(row).second) {
        throw std::invalid_argument(
          "Phase 4a selected rows are invalid or duplicated.");
      }
      rows0.push_back(row - 1);
    }
    row_pointer = rows0.data();
    selected_sample_count = static_cast<int>(rows0.size());
  }
  if (phenotype.nrow() != selected_sample_count || phenotype.ncol() != 2) {
    throw std::invalid_argument(
      "Phase 4a phenotype must have selected samples by exactly two traits.");
  }
  for (double value : phenotype) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("Phase 4a phenotype must be complete and finite.");
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
    const std::vector<arma::mat>& draws) {
  const int ndraw = static_cast<int>(draws.size());
  const int markers = ndraw == 0 ? 0 : static_cast<int>(draws[0].n_rows);
  Rcpp::NumericVector out(static_cast<R_xlen_t>(ndraw) * markers * 2);
  out.attr("dim") = Rcpp::IntegerVector::create(ndraw, markers, 2);
  for (int trait = 0; trait < 2; ++trait) {
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
    const std::vector<std::vector<int>>& draws) {
  const int ndraw = static_cast<int>(draws.size());
  const int markers = ndraw == 0 ? 0 : static_cast<int>(draws[0].size());
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
    const std::vector<std::array<double, 4>>& draws) {
  Rcpp::NumericMatrix out(static_cast<int>(draws.size()), 4);
  for (std::size_t draw = 0u; draw < draws.size(); ++draw) {
    for (std::size_t state = 0u; state < 4u; ++state) {
      out(static_cast<int>(draw), static_cast<int>(state)) = draws[draw][state];
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
  Rcpp::NumericMatrix transitions(4, 4);
  for (std::size_t from = 0u; from < 4u; ++from) {
    for (std::size_t to = 0u; to < 4u; ++to) {
      transitions(static_cast<int>(from), static_cast<int>(to)) =
        static_cast<double>(chain.transition_count[from][to]);
    }
  }
  return Rcpp::List::create(
    Rcpp::_ ["realised_effects"] = phase4a_effect_draws(chain.realised_draws),
    Rcpp::_ ["latent_effects"] = phase4a_effect_draws(chain.latent_draws),
    Rcpp::_ ["joint_states"] = phase4a_state_draws(chain.state_draws),
    Rcpp::_ ["marker_covariance"] = phase4a_matrix_draws(
      chain.covariance_draws, 2, 2),
    Rcpp::_ ["activity_pattern_parameters"] = phase4a_probability_draws(
      chain.probability_draws),
    Rcpp::_ ["predictions"] = phase4a_matrix_draws(
      chain.prediction_draws, static_cast<int>(chain.final_prediction.n_rows), 2),
    Rcpp::_ ["convergence_marker_covariance"] = phase4a_matrix_draws(
      chain.convergence_covariance, 2, 2),
    Rcpp::_ ["convergence_activity_pattern_parameters"] =
      phase4a_probability_draws(chain.convergence_probability),
    Rcpp::_ ["convergence_active_marker_count"] =
      Rcpp::wrap(chain.convergence_active_count),
    Rcpp::_ ["final_realised_effects"] = phase4a_effect_matrix(
      chain.final_realised),
    Rcpp::_ ["final_latent_effects"] = phase4a_effect_matrix(
      chain.final_latent),
    Rcpp::_ ["final_joint_states"] = Rcpp::wrap(chain.final_state),
    Rcpp::_ ["final_marker_covariance"] = Rcpp::wrap(chain.final_covariance),
    Rcpp::_ ["final_activity_pattern_parameters"] =
      Rcpp::wrap(chain.final_probability),
    Rcpp::_ ["final_predictions"] = Rcpp::wrap(chain.final_prediction),
    Rcpp::_ ["last_covariance_statistic"] = Rcpp::wrap(
      chain.last_covariance_statistic),
    Rcpp::_ ["last_covariance_scale"] = Rcpp::wrap(
      chain.last_covariance_scale),
    Rcpp::_ ["last_covariance_degrees_of_freedom"] =
      chain.last_covariance_df,
    Rcpp::_ ["last_active_marker_count"] = chain.last_active_count,
    Rcpp::_ ["transition_counts"] = transitions);
}

}  // namespace

// Deterministic Phase 4a marker conditional oracle. It performs no draw and
// does not access a genotype resource.
// [[Rcpp::export]]
Rcpp::List mtblr_phase4a_pattern_contract_internal(
    Rcpp::NumericVector score,
    double marker_sum_squares,
    arma::mat marker_covariance,
    arma::mat residual_covariance,
    Rcpp::NumericVector activity_pattern_probability) {
  if (score.size() != 2 || activity_pattern_probability.size() != 4) {
    throw std::invalid_argument(
      "Phase 4a marker inspection requires two scores and four probabilities.");
  }
  const arma::mat Vb = sblr::phase4a::require_spd(
    marker_covariance, "marker covariance");
  const arma::mat Ve = sblr::phase4a::require_spd(
    residual_covariance, "residual covariance");
  std::array<double, 4> probability;
  double total = 0.0;
  for (std::size_t state = 0u; state < 4u; ++state) {
    const double value = activity_pattern_probability[
      static_cast<R_xlen_t>(state)];
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument(
        "Phase 4a activity-pattern probabilities must be finite and positive.");
    }
    probability[state] = value;
    total += value;
  }
  for (double& value : probability) value /= total;
  const auto kernel = sblr::phase4a::pattern_kernel(
    Rcpp::as<arma::vec>(score), marker_sum_squares, Vb,
    sblr::phase4a::inverse_spd(Ve), probability);
  Rcpp::List means(4), covariances(4);
  for (std::size_t state = 1u; state < 4u; ++state) {
    means[static_cast<R_xlen_t>(state)] = Rcpp::wrap(kernel.active_mean[state]);
    covariances[static_cast<R_xlen_t>(state)] =
      Rcpp::wrap(kernel.active_covariance[state]);
  }
  return Rcpp::List::create(
    Rcpp::_ ["pattern_order"] = Rcpp::CharacterVector::create(
      "0_0", "1_0", "0_1", "1_1"),
    Rcpp::_ ["probability"] = Rcpp::wrap(kernel.probability),
    Rcpp::_ ["log_weight"] = Rcpp::wrap(kernel.log_weight),
    Rcpp::_ ["active_mean"] = means,
    Rcpp::_ ["active_covariance"] = covariances,
    Rcpp::_ ["completion_coefficient"] = Rcpp::NumericVector::create(
      Vb(1, 0) / Vb(0, 0), Vb(0, 1) / Vb(1, 1)),
    Rcpp::_ ["completion_variance"] = Rcpp::NumericVector::create(
      Vb(1, 1) - Vb(1, 0) * Vb(0, 1) / Vb(0, 0),
      Vb(0, 0) - Vb(0, 1) * Vb(1, 0) / Vb(1, 1)));
}

// Qualification-only two-trait common-sample Cheng MT-BayesC-Pi route.
// It is deliberately not namespace-exported and does not call the legacy MT
// covariance implementation.
// [[Rcpp::export]]
Rcpp::List mtblr_phase4a_cheng_bed_internal(
    Rcpp::CharacterVector bed_files,
    int source_sample_count,
    Rcpp::List selected_columns,
    Rcpp::Nullable<Rcpp::IntegerVector> selected_rows,
    Rcpp::NumericVector allele_frequency,
    Rcpp::NumericMatrix phenotype,
    arma::mat fixed_residual_covariance,
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
    Rcpp::Nullable<Rcpp::List> execution_contract) {
  if (phenotype.ncol() != 2 || chains <= 0 || cores <= 0 ||
      initial_activity_pattern_probability.size() != 4 ||
      activity_pattern_dirichlet_prior.size() != 4) {
    throw std::invalid_argument(
      "Phase 4a requires two traits, four patterns, and positive chain controls.");
  }
  std::array<double, 4> initial_probability;
  std::array<double, 4> dirichlet_prior;
  double probability_total = 0.0;
  for (std::size_t state = 0u; state < 4u; ++state) {
    const double value = initial_activity_pattern_probability[
      static_cast<R_xlen_t>(state)];
    const double prior = activity_pattern_dirichlet_prior[
      static_cast<R_xlen_t>(state)];
    if (!std::isfinite(value) || value <= 0.0 ||
        !std::isfinite(prior) || prior <= 0.0) {
      throw std::invalid_argument(
        "Phase 4a probabilities and Dirichlet shapes must be positive.");
    }
    initial_probability[state] = value;
    probability_total += value;
    dirichlet_prior[state] = prior;
  }
  for (double& value : initial_probability) value /= probability_total;

  Phase4aPreparedBed prepared = prepare_phase4a_bed(
    bed_files, source_sample_count, selected_columns, selected_rows,
    allele_frequency, phenotype);
  const auto genotype = prepared.view();
  const BlrPhase3ExecutionContract phase3 =
    parse_blr_phase3_execution_contract(
      execution_contract, chains, sampling_iterations);
  if (!phase3.active()) {
    throw std::invalid_argument(
      "Phase 4a requires the active Phase 3 execution contract.");
  }
  for (int chain = 0; chain < chains; ++chain) {
    const std::string expected = "chain:" + std::to_string(chain);
    if (phase3.task_ids[static_cast<std::size_t>(chain)] != expected) {
      throw std::invalid_argument(
        "Phase 4a task IDs must follow canonical joint-multitrait chain order.");
    }
  }

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
        genotype, prepared.marker_maps, prepared.phenotype,
        fixed_residual_covariance, initial_marker_covariance,
        initial_probability, dirichlet_prior, marker_covariance_prior_df,
        marker_covariance_prior_scale, update_marker_covariance,
        update_activity_pattern_probability, burn_in_iterations,
        sampling_iterations, phase3.retained_transition_indices,
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
        "Phase 4a chain " + std::to_string(chain + 1) + " failed: " +
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
    Rcpp::_ ["qualification_only"] = true);
}
