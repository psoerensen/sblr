#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "blr_result.h"
#include "blr_scalar_execution.h"
#include "blr_csr_bayesr_types.h"
#include "blr_csr_sbayesrc_types.h"
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

// Test-only binding for the binding-neutral Phase 4 task contract.
// [[Rcpp::export]]
Rcpp::DataFrame blr_phase4_scalar_tasks_cpp(int ntraits, int nchains) {
  using namespace sblr::core;
  if (ntraits < 0 || nchains < 0) {
    throw std::invalid_argument("task dimensions must be non-negative");
  }
  const std::vector<ScalarChainTask> tasks = make_scalar_chain_tasks(
    static_cast<std::size_t>(ntraits), static_cast<std::size_t>(nchains)
  );
  Rcpp::IntegerVector trait(tasks.size());
  Rcpp::IntegerVector chain(tasks.size());
  Rcpp::IntegerVector task(tasks.size());
  for (std::size_t index = 0; index < tasks.size(); ++index) {
    trait[index] = static_cast<int>(tasks[index].trait_index);
    chain[index] = static_cast<int>(tasks[index].chain_index);
    task[index] = static_cast<int>(tasks[index].task_index);
  }
  return Rcpp::DataFrame::create(
    Rcpp::_["trait"] = trait,
    Rcpp::_["chain"] = chain,
    Rcpp::_["task"] = task
  );
}

// Test-only binding for exact Phase 4 seed resolution.
// [[Rcpp::export]]
Rcpp::IntegerVector blr_phase4_scalar_seeds_cpp(
  int seed,
  int ntraits,
  int nchains,
  Rcpp::IntegerVector chain_seeds
) {
  using namespace sblr::core;
  if (ntraits < 0 || nchains < 0) {
    throw std::invalid_argument("seed dimensions must be non-negative");
  }
  const std::vector<int> explicit_seeds =
    Rcpp::as<std::vector<int>>(chain_seeds);
  const std::vector<ScalarChainTask> tasks = make_scalar_chain_tasks(
    static_cast<std::size_t>(ntraits), static_cast<std::size_t>(nchains)
  );
  validate_scalar_chain_seeds(static_cast<std::size_t>(nchains), explicit_seeds);
  Rcpp::IntegerVector resolved(tasks.size());
  for (std::size_t index = 0; index < tasks.size(); ++index) {
    resolved[index] = static_cast<int>(resolve_scalar_chain_seed(
      seed, static_cast<std::size_t>(nchains), explicit_seeds, tasks[index]
    ));
  }
  return resolved;
}

// Test-only binding for the common retained-iteration predicate.
// [[Rcpp::export]]
Rcpp::LogicalVector blr_phase4_retained_iterations_cpp(
  int trace_length,
  int burnin,
  int thinning
) {
  using namespace sblr::core;
  if (trace_length < 0 || burnin < 0 || thinning <= 0) {
    throw std::invalid_argument("retained-iteration controls are invalid");
  }
  Rcpp::LogicalVector retained(trace_length);
  for (int iteration = 0; iteration < trace_length; ++iteration) {
    retained[iteration] = scalar_iteration_is_retained(
      iteration, burnin, thinning
    );
  }
  return retained;
}

// Internal Phase 5A contract bridge. It validates and round-trips metadata;
// it deliberately has no sampler or CSR-reader call.
// [[Rcpp::export]]
Rcpp::List blr_phase5a_validate_bayesr_contract_cpp(Rcpp::List spec) {
  using namespace sblr::core;
  Rcpp::List data = spec["data"], component = spec["component"];
  Rcpp::List controls = spec["controls"], output = spec["output"];
  const std::size_t m = Rcpp::as<std::size_t>(data["marker_count"]);
  const std::size_t nt = Rcpp::as<std::size_t>(data["trait_count"]);
  std::vector<std::uint64_t> row_ptr(m + 1, 0);
  std::vector<std::uint32_t> col_index;
  std::vector<float> values;
  std::vector<double> diagonal(m, 1.0);
  std::vector<int> sample_size = Rcpp::as<std::vector<int>>(data["sample_size"]);
  CsrBayesRExecutionInput x;
  x.data.marker_count=m; x.data.trait_count=nt; x.data.row_ptr=row_ptr.data();
  x.data.row_ptr_count=row_ptr.size(); x.data.column_index=col_index.data();
  x.data.value=values.data(); x.data.nonzero_count=0; x.data.diagonal=diagonal.data();
  x.data.diagonal_count=diagonal.size(); x.data.sample_size=sample_size.data();
  x.data.sample_size_count=sample_size.size();
  x.data.shared_read_only=Rcpp::as<bool>(data["shared_read_only"]);
  x.data.per_chain_payload=Rcpp::as<bool>(data["per_chain_payload"]);
  x.data.storage_outlives_execution=Rcpp::as<bool>(data["storage_outlives_execution"]);
  x.marker_order=Rcpp::as<std::vector<std::string>>(data["marker_order"]);
  x.trait_order=Rcpp::as<std::vector<std::string>>(data["trait_order"]);
  x.components.scales=Rcpp::as<std::vector<double>>(component["scales"]);
  x.components.initial_probability=Rcpp::as<std::vector<double>>(component["initial_probability"]);
  x.components.dirichlet_prior=Rcpp::as<std::vector<double>>(component["dirichlet_prior"]);
  x.components.null_component=Rcpp::as<std::size_t>(component["null_component"]);
  x.components.update_probability=Rcpp::as<bool>(component["update_probability"]);
  x.components.scale_interpretation=Rcpp::as<std::string>(component["scale_interpretation"]);
  x.controls.iterations=Rcpp::as<int>(controls["iterations"]); x.controls.burnin=Rcpp::as<int>(controls["burnin"]);
  x.controls.thinning=Rcpp::as<int>(controls["thinning"]); x.controls.chains=Rcpp::as<int>(controls["chains"]);
  x.controls.cores=Rcpp::as<int>(controls["cores"]); x.controls.seed=Rcpp::as<int>(controls["seed"]);
  x.controls.chain_seeds=Rcpp::as<std::vector<int>>(controls["chain_seeds"]);
  x.controls.keep_chains=Rcpp::as<bool>(controls["keep_chains"]);
  x.controls.update_marker_variance=Rcpp::as<bool>(controls["update_marker_variance"]);
  x.controls.update_residual_variance=Rcpp::as<bool>(controls["update_residual_variance"]);
  x.controls.update_ld_swap=Rcpp::as<bool>(controls["update_ld_swap"]);
  x.controls.ld_swap_probability=Rcpp::as<double>(controls["ld_swap_probability"]);
  x.controls.ld_swap_r2=Rcpp::as<double>(controls["ld_swap_r2"]);
  x.controls.ld_swap_max_friends=Rcpp::as<int>(controls["ld_swap_max_friends"]);
  x.controls.ld_swap_moves=Rcpp::as<int>(controls["ld_swap_moves"]);
  x.output.keep_chains=Rcpp::as<bool>(output["keep_chains"]);
  validate_csr_bayesr_execution_input(x);
  return Rcpp::List::create(Rcpp::Named("data")=data,Rcpp::Named("component")=component,
    Rcpp::Named("controls")=controls,Rcpp::Named("output")=output,
    Rcpp::Named("validated")=true,Rcpp::Named("invokes_sampler")=false);
}

// Internal Phase 7A validation-only bridge. It owns small test buffers while
// validating borrowed views and never calls the SBayesRC sampler.
// [[Rcpp::export]]
Rcpp::List blr_phase7a_validate_sbayesrc_contract_cpp(Rcpp::List spec) {
  using namespace sblr::core;
  Rcpp::List data=spec["data"], annotation=spec["annotation"], component=spec["component"];
  Rcpp::List alpha=spec["alpha"], probability=spec["probability"], controls=spec["controls"], output=spec["output"];
  const std::size_t m=Rcpp::as<std::size_t>(data["marker_count"]), nt=Rcpp::as<std::size_t>(data["trait_count"]);
  std::vector<std::uint64_t> row_ptr(m+1,0); std::vector<double> diagonal(m,1.0);
  std::vector<int> sample_size=Rcpp::as<std::vector<int>>(data["sample_size"]);
  Rcpp::NumericMatrix A=annotation["values"];
  CsrSBayesRCExecutionInput x;
  x.data.marker_count=m; x.data.trait_count=nt; x.data.row_ptr=row_ptr.data(); x.data.row_ptr_count=row_ptr.size();
  x.data.diagonal=diagonal.data(); x.data.diagonal_count=diagonal.size(); x.data.sample_size=sample_size.data(); x.data.sample_size_count=sample_size.size();
  x.data.shared_read_only=Rcpp::as<bool>(data["shared_read_only"]); x.data.per_chain_payload=Rcpp::as<bool>(data["per_chain_payload"]); x.data.storage_outlives_execution=Rcpp::as<bool>(data["storage_outlives_execution"]);
  x.marker_order=Rcpp::as<std::vector<std::string>>(data["marker_order"]); x.trait_order=Rcpp::as<std::vector<std::string>>(data["trait_order"]);
  x.annotation.marker_count=Rcpp::as<std::size_t>(annotation["marker_count"]); x.annotation.annotation_count=Rcpp::as<std::size_t>(annotation["annotation_count"]);
  x.annotation.values=A.begin(); x.annotation.value_count=static_cast<std::size_t>(A.size()); x.annotation.annotation_order=Rcpp::as<std::vector<std::string>>(annotation["annotation_order"]);
  x.annotation.layout=Rcpp::as<std::string>(annotation["layout"]); x.annotation.includes_intercept=Rcpp::as<bool>(annotation["includes_intercept"]);
  x.annotation.standardized=Rcpp::as<bool>(annotation["standardized"]); x.annotation.centered_binary=Rcpp::as<bool>(annotation["centered_binary"]);
  x.annotation.shared_read_only=Rcpp::as<bool>(annotation["shared_read_only"]); x.annotation.per_chain_payload=Rcpp::as<bool>(annotation["per_chain_payload"]); x.annotation.storage_outlives_execution=Rcpp::as<bool>(annotation["storage_outlives_execution"]);
  x.components.scales=Rcpp::as<std::vector<double>>(component["scales"]); x.components.null_component=Rcpp::as<std::size_t>(component["null_component"]); x.components.scale_interpretation=Rcpp::as<std::string>(component["scale_interpretation"]);
  x.alpha.annotation_count=Rcpp::as<std::size_t>(alpha["annotation_count"]); x.alpha.step_count=Rcpp::as<std::size_t>(alpha["step_count"]);
  x.alpha.initial_values=Rcpp::as<std::vector<double>>(alpha["initial_values"]); x.alpha.initial_variance=Rcpp::as<std::vector<double>>(alpha["initial_variance"]);
  x.alpha.intercept_flat=Rcpp::as<bool>(alpha["intercept_flat"]); x.alpha.update=Rcpp::as<bool>(alpha["update"]); x.alpha.variance_prior_a=Rcpp::as<double>(alpha["variance_prior_a"]); x.alpha.variance_prior_b=Rcpp::as<double>(alpha["variance_prior_b"]); x.alpha.update_every=Rcpp::as<int>(alpha["update_every"]);
  x.probability.transformation=Rcpp::as<std::string>(probability["transformation"]); x.probability.reference_component=Rcpp::as<std::size_t>(probability["reference_component"]); x.probability.stick_order=Rcpp::as<std::vector<std::size_t>>(probability["stick_order"]); x.probability.probability_floor=Rcpp::as<double>(probability["probability_floor"]); x.probability.floor_then_normalize=Rcpp::as<bool>(probability["floor_then_normalize"]);
  x.controls.iterations=Rcpp::as<int>(controls["iterations"]); x.controls.burnin=Rcpp::as<int>(controls["burnin"]); x.controls.thinning=Rcpp::as<int>(controls["thinning"]); x.controls.chains=Rcpp::as<int>(controls["chains"]); x.controls.cores=Rcpp::as<int>(controls["cores"]); x.controls.seed=Rcpp::as<int>(controls["seed"]); x.controls.chain_seeds=Rcpp::as<std::vector<int>>(controls["chain_seeds"]); x.controls.keep_chains=Rcpp::as<bool>(controls["keep_chains"]); x.controls.update_ld_swap=Rcpp::as<bool>(controls["update_ld_swap"]); x.controls.ld_swap_probability=Rcpp::as<double>(controls["ld_swap_probability"]); x.controls.ld_swap_r2=Rcpp::as<double>(controls["ld_swap_r2"]); x.controls.ld_swap_max_friends=Rcpp::as<int>(controls["ld_swap_max_friends"]); x.controls.ld_swap_moves=Rcpp::as<int>(controls["ld_swap_moves"]);
  x.output.keep_chains=Rcpp::as<bool>(output["keep_chains"]); x.output.diagnostics=Rcpp::as<bool>(output["diagnostics"]);
  validate_csr_sbayesrc_execution_input(x);
  return Rcpp::List::create(Rcpp::Named("data")=data,Rcpp::Named("annotation")=annotation,Rcpp::Named("component")=component,Rcpp::Named("alpha")=alpha,Rcpp::Named("probability")=probability,Rcpp::Named("controls")=controls,Rcpp::Named("output")=output,Rcpp::Named("validated")=true,Rcpp::Named("invokes_sampler")=false);
}
