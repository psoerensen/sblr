#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "blr_result.h"
#include "blr_spec.h"

namespace {

using namespace sblr::core;

Rcpp::List list_field(const Rcpp::List& value, const char* field) {
  if (!value.containsElementNamed(field) || Rf_isNull(value[field])) {
    throw std::invalid_argument(std::string(field) + " must be a list");
  }
  return Rcpp::as<Rcpp::List>(value[field]);
}

std::string string_field(const Rcpp::List& value, const char* field,
                         const std::string& path) {
  if (!value.containsElementNamed(field) || Rf_isNull(value[field])) {
    throw std::invalid_argument(path + " must be present");
  }
  Rcpp::CharacterVector input(value[field]);
  if (input.size() != 1 || input[0] == NA_STRING) {
    throw std::invalid_argument(path + " must be a character scalar");
  }
  const std::string output = Rcpp::as<std::string>(input[0]);
  if (output.empty()) {
    throw std::invalid_argument(path + " must not be empty");
  }
  return output;
}

int integer_field(const Rcpp::List& value, const char* field,
                  const std::string& path) {
  if (!value.containsElementNamed(field) || Rf_isNull(value[field])) {
    throw std::invalid_argument(path + " must be present");
  }
  Rcpp::NumericVector input(value[field]);
  if (input.size() != 1 || Rcpp::NumericVector::is_na(input[0]) ||
      !std::isfinite(input[0]) || input[0] != std::floor(input[0]) ||
      input[0] < static_cast<double>(std::numeric_limits<int>::min()) ||
      input[0] > static_cast<double>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(path + " must be an integer scalar");
  }
  return static_cast<int>(input[0]);
}

std::size_t size_field(const Rcpp::List& value, const char* field,
                       const std::string& path) {
  const int input = integer_field(value, field, path);
  if (input < 0) throw std::invalid_argument(path + " must be non-negative");
  return static_cast<std::size_t>(input);
}

bool bool_field(const Rcpp::List& value, const char* field,
                const std::string& path) {
  if (!value.containsElementNamed(field) || Rf_isNull(value[field])) {
    throw std::invalid_argument(path + " must be present");
  }
  Rcpp::LogicalVector input(value[field]);
  if (input.size() != 1 || input[0] == NA_LOGICAL) {
    throw std::invalid_argument(path + " must be true or false");
  }
  return input[0] == TRUE;
}

void require_tag(const std::string& actual, const char* expected,
                 const std::string& path) {
  if (actual != expected) {
    throw std::invalid_argument(path + " has an unsupported Phase 1 tag");
  }
}

ResolvedSpec resolved_spec_from_r(const Rcpp::List& input) {
  const Rcpp::List schema = list_field(input, "schema");
  const Rcpp::List data = list_field(input, "data");
  const Rcpp::List model = list_field(input, "model");
  const Rcpp::List mcmc = list_field(input, "mcmc");
  const Rcpp::List output = list_field(input, "output");
  const Rcpp::List execution = list_field(input, "execution");
  const Rcpp::List csr = list_field(data, "csr");

  ResolvedSpec spec;
  spec.schema_name = string_field(schema, "name", "schema$name");
  spec.schema_version = integer_field(schema, "version", "schema$version");

  require_tag(string_field(data, "representation", "data$representation"),
              "csr", "data$representation");
  require_tag(string_field(data, "design", "data$design"),
              "independent_traits", "data$design");
  require_tag(string_field(data, "scaling", "data$scaling"),
              "standardized_genotype", "data$scaling");
  spec.data.marker_count = size_field(data, "n_markers", "data$n_markers");
  spec.data.trait_count = size_field(data, "n_traits", "data$n_traits");
  spec.data.marker_ids = Rcpp::as<std::vector<std::string>>(data["marker_ids"]);
  spec.data.trait_ids = Rcpp::as<std::vector<std::string>>(data["trait_ids"]);
  spec.data.sample_size = Rcpp::as<std::vector<int>>(data["sample_size"]);
  spec.data.csr.resource_id =
    string_field(csr, "resource_id", "data$csr$resource_id");
  spec.data.csr.marker_count =
    size_field(csr, "marker_count", "data$csr$marker_count");
  spec.data.csr.shared_read_only =
    bool_field(csr, "shared_read_only", "data$csr$shared_read_only");
  spec.data.csr.per_chain_data =
    bool_field(csr, "per_chain_data", "data$csr$per_chain_data");
  spec.data.csr.lifetime_exceeds_chains = bool_field(
    csr, "lifetime_exceeds_chains", "data$csr$lifetime_exceeds_chains"
  );

  require_tag(string_field(model, "kernel", "model$kernel"),
              "scalar", "model$kernel");
  require_tag(string_field(model, "family", "model$family"),
              "bayesc", "model$family");
  require_tag(string_field(model, "state", "model$state"),
              "binary", "model$state");
  require_tag(string_field(model, "probability", "model$probability"),
              "global_binary", "model$probability");
  require_tag(string_field(model, "scale", "model$scale"),
              "unit", "model$scale");
  require_tag(string_field(model, "trait_covariance",
                           "model$trait_covariance"),
              "scalar_independent", "model$trait_covariance");
  require_tag(string_field(model, "residual_covariance",
                           "model$residual_covariance"),
              "scalar_independent", "model$residual_covariance");

  spec.mcmc.nit = integer_field(mcmc, "nit", "mcmc$nit");
  spec.mcmc.nburn = integer_field(mcmc, "nburn", "mcmc$nburn");
  spec.mcmc.nthin = integer_field(mcmc, "nthin", "mcmc$nthin");
  spec.mcmc.nchains = integer_field(mcmc, "nchains", "mcmc$nchains");
  spec.mcmc.ncores = integer_field(mcmc, "ncores", "mcmc$ncores");
  spec.mcmc.seed = integer_field(mcmc, "seed", "mcmc$seed");
  spec.mcmc.has_explicit_chain_seeds =
    mcmc.containsElementNamed("chain_seeds") &&
    !Rf_isNull(mcmc["chain_seeds"]);
  if (spec.mcmc.has_explicit_chain_seeds) {
    spec.mcmc.chain_seeds =
      Rcpp::as<std::vector<int>>(mcmc["chain_seeds"]);
  }

  spec.output.marker_mean =
    bool_field(output, "marker_mean", "output$marker_mean");
  spec.output.marker_pip =
    bool_field(output, "marker_pip", "output$marker_pip");
  spec.output.parameter_traces = bool_field(
    output, "parameter_traces", "output$parameter_traces"
  );
  spec.output.final_state =
    bool_field(output, "final_state", "output$final_state");
  spec.output.keep_chain_summaries = bool_field(
    output, "keep_chain_summaries", "output$keep_chain_summaries"
  );

  require_tag(string_field(execution, "operator", "execution$operator"),
              "csr", "execution$operator");
  require_tag(string_field(execution, "backend_reference",
                           "execution$backend_reference"),
              "stblr_cpg_omp_csr", "execution$backend_reference");
  spec.execution.scheduled =
    bool_field(execution, "scheduled", "execution$scheduled");

  validate_resolved_spec(spec);
  return spec;
}

Rcpp::List resolved_spec_to_r(const ResolvedSpec& spec) {
  Rcpp::RObject chain_seeds = R_NilValue;
  if (spec.mcmc.has_explicit_chain_seeds) {
    chain_seeds = Rcpp::wrap(spec.mcmc.chain_seeds);
  }
  return Rcpp::List::create(
    Rcpp::_["schema"] = Rcpp::List::create(
      Rcpp::_["name"] = spec.schema_name,
      Rcpp::_["version"] = spec.schema_version
    ),
    Rcpp::_["data"] = Rcpp::List::create(
      Rcpp::_["representation"] = "csr",
      Rcpp::_["design"] = "independent_traits",
      Rcpp::_["n_markers"] = static_cast<int>(spec.data.marker_count),
      Rcpp::_["n_traits"] = static_cast<int>(spec.data.trait_count),
      Rcpp::_["marker_ids"] = spec.data.marker_ids,
      Rcpp::_["trait_ids"] = spec.data.trait_ids,
      Rcpp::_["sample_size"] = spec.data.sample_size,
      Rcpp::_["scaling"] = "standardized_genotype",
      Rcpp::_["csr"] = Rcpp::List::create(
        Rcpp::_["resource_id"] = spec.data.csr.resource_id,
        Rcpp::_["marker_count"] =
          static_cast<int>(spec.data.csr.marker_count),
        Rcpp::_["shared_read_only"] = spec.data.csr.shared_read_only,
        Rcpp::_["per_chain_data"] = spec.data.csr.per_chain_data,
        Rcpp::_["lifetime_exceeds_chains"] =
          spec.data.csr.lifetime_exceeds_chains
      )
    ),
    Rcpp::_["model"] = Rcpp::List::create(
      Rcpp::_["kernel"] = "scalar",
      Rcpp::_["family"] = "bayesc",
      Rcpp::_["state"] = "binary",
      Rcpp::_["probability"] = "global_binary",
      Rcpp::_["scale"] = "unit",
      Rcpp::_["trait_covariance"] = "scalar_independent",
      Rcpp::_["residual_covariance"] = "scalar_independent"
    ),
    Rcpp::_["mcmc"] = Rcpp::List::create(
      Rcpp::_["nit"] = spec.mcmc.nit,
      Rcpp::_["nburn"] = spec.mcmc.nburn,
      Rcpp::_["nthin"] = spec.mcmc.nthin,
      Rcpp::_["nchains"] = spec.mcmc.nchains,
      Rcpp::_["ncores"] = spec.mcmc.ncores,
      Rcpp::_["seed"] = spec.mcmc.seed,
      Rcpp::_["chain_seeds"] = chain_seeds
    ),
    Rcpp::_["output"] = Rcpp::List::create(
      Rcpp::_["marker_mean"] = spec.output.marker_mean,
      Rcpp::_["marker_pip"] = spec.output.marker_pip,
      Rcpp::_["parameter_traces"] = spec.output.parameter_traces,
      Rcpp::_["final_state"] = spec.output.final_state,
      Rcpp::_["keep_chain_summaries"] = spec.output.keep_chain_summaries
    ),
    Rcpp::_["execution"] = Rcpp::List::create(
      Rcpp::_["operator"] = "csr",
      Rcpp::_["backend_reference"] = "stblr_cpg_omp_csr",
      Rcpp::_["scheduled"] = spec.execution.scheduled
    )
  );
}

}  // namespace

// Specification-only conversion boundary. No execution backend is called.
// [[Rcpp::export]]
Rcpp::List blr_phase1_validate_spec_cpp(Rcpp::List spec) {
  return resolved_spec_to_r(resolved_spec_from_r(spec));
}

// Constructs a typed result vocabulary solely to validate canonical shapes.
// [[Rcpp::export]]
Rcpp::List blr_phase1_validate_result_dimensions_cpp(Rcpp::List dimensions) {
  using namespace sblr::core;
  const std::size_t markers =
    size_field(dimensions, "markers", "result$markers");
  const std::size_t traits =
    size_field(dimensions, "traits", "result$traits");
  const std::size_t retained =
    size_field(dimensions, "retained_samples", "result$retained_samples");
  const std::size_t parameters =
    size_field(dimensions, "parameter_dimension",
               "result$parameter_dimension");
  const std::size_t marker_effect_length = size_field(
    dimensions, "marker_effect_length", "result$marker_effect_length"
  );
  const std::size_t marker_pip_length = size_field(
    dimensions, "marker_pip_length", "result$marker_pip_length"
  );
  const std::size_t final_effect_length = size_field(
    dimensions, "final_effect_length", "result$final_effect_length"
  );
  const std::size_t final_state_length = size_field(
    dimensions, "final_state_length", "result$final_state_length"
  );
  const std::size_t trace_length =
    size_field(dimensions, "trace_length", "result$trace_length");
  const std::size_t trait_covariance_length = size_field(
    dimensions, "trait_covariance_length",
    "result$trait_covariance_length"
  );
  const std::size_t residual_covariance_length = size_field(
    dimensions, "residual_covariance_length",
    "result$residual_covariance_length"
  );
  const bool has_component = bool_field(
    dimensions, "has_component_probability",
    "result$has_component_probability"
  );
  const std::size_t components = size_field(
    dimensions, "components", "result$components"
  );
  const bool has_pattern = bool_field(
    dimensions, "has_pattern_probability",
    "result$has_pattern_probability"
  );
  const std::size_t patterns =
    size_field(dimensions, "patterns", "result$patterns");

  BlrResult result;
  result.marker.marker_count = markers;
  result.marker.trait_count = traits;
  result.marker.effect_mean.resize(marker_effect_length);
  result.marker.pip.resize(marker_pip_length);
  result.state.final_effect.resize(final_effect_length);
  result.state.final_state.resize(final_state_length);
  result.trace.retained_samples = retained;
  result.trace.parameter_count = parameters;
  result.trace.values.resize(trace_length);
  result.variance.trait_count = traits;
  result.variance.trait_covariance.resize(trait_covariance_length);
  result.variance.residual_covariance.resize(residual_covariance_length);
  result.diagnostics.retained_samples = retained;
  result.optional.has_component_probability = has_component;
  result.optional.component_count = components;
  result.optional.has_pattern_probability = has_pattern;
  result.optional.pattern_count = patterns;
  validate_blr_result(result);

  Rcpp::RObject component_dimensions = R_NilValue;
  if (has_component) {
    component_dimensions = Rcpp::IntegerVector::create(
      static_cast<int>(markers), static_cast<int>(components),
      static_cast<int>(traits)
    );
  }
  Rcpp::RObject pattern_dimensions = R_NilValue;
  if (has_pattern) {
    pattern_dimensions = Rcpp::IntegerVector::create(
      static_cast<int>(markers), static_cast<int>(patterns)
    );
  }
  return Rcpp::List::create(
    Rcpp::_["marker_effect"] = Rcpp::IntegerVector::create(
      static_cast<int>(markers), static_cast<int>(traits)
    ),
    Rcpp::_["marker_pip"] = Rcpp::IntegerVector::create(
      static_cast<int>(markers), static_cast<int>(traits)
    ),
    Rcpp::_["trace"] = Rcpp::IntegerVector::create(
      static_cast<int>(retained), static_cast<int>(parameters)
    ),
    Rcpp::_["trait_covariance"] = Rcpp::IntegerVector::create(
      static_cast<int>(traits), static_cast<int>(traits)
    ),
    Rcpp::_["residual_covariance"] = Rcpp::IntegerVector::create(
      static_cast<int>(traits), static_cast<int>(traits)
    ),
    Rcpp::_["component_probability"] = component_dimensions,
    Rcpp::_["pattern_probability"] = pattern_dimensions
  );
}
