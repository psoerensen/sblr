#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "blr_result.h"
#include "blr_scalar_execution.h"
#include "blr_csr_bayesr_types.h"
#include "blr_csr_sbayesrc_types.h"
#include "blr_csr_annotation_bayesc_types.h"
#include "blr_scheduled_execution_types.h"
#include "blr_spec.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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

namespace {
struct Phase9ACommonBuffers {
 std::vector<std::uint64_t> row_ptr; std::vector<double> diagonal; std::vector<int> sample_size;
 sblr::core::AnnotationBayesCDataView data; sblr::core::AnnotationBayesCControls controls;
};
Phase9ACommonBuffers phase9a_common(Rcpp::List data, Rcpp::List controls) {
 Phase9ACommonBuffers b; const std::size_t m=Rcpp::as<std::size_t>(data["marker_count"]), nt=Rcpp::as<std::size_t>(data["trait_count"]);
 b.row_ptr.assign(m+1,0); b.diagonal.assign(m,1.0); b.sample_size=Rcpp::as<std::vector<int>>(data["sample_size"]);
 b.data.marker_count=m; b.data.trait_count=nt; b.data.row_ptr=b.row_ptr.data(); b.data.row_ptr_count=b.row_ptr.size(); b.data.diagonal=b.diagonal.data(); b.data.diagonal_count=b.diagonal.size(); b.data.sample_size=b.sample_size.data(); b.data.sample_size_count=b.sample_size.size();
 b.data.shared_read_only=Rcpp::as<bool>(data["shared_read_only"]); b.data.per_chain_payload=Rcpp::as<bool>(data["per_chain_payload"]); b.data.storage_outlives_execution=Rcpp::as<bool>(data["storage_outlives_execution"]);
 b.controls.iterations=Rcpp::as<int>(controls["iterations"]); b.controls.burnin=Rcpp::as<int>(controls["burnin"]); b.controls.thinning=Rcpp::as<int>(controls["thinning"]); b.controls.chains=Rcpp::as<int>(controls["chains"]); b.controls.cores=Rcpp::as<int>(controls["cores"]); b.controls.seed=Rcpp::as<int>(controls["seed"]); b.controls.chain_seeds=Rcpp::as<std::vector<int>>(controls["chain_seeds"]); b.controls.keep_chains=Rcpp::as<bool>(controls["keep_chains"]); b.controls.update_marker_variance=Rcpp::as<bool>(controls["update_marker_variance"]); b.controls.update_residual_variance=Rcpp::as<bool>(controls["update_residual_variance"]); b.controls.update_global_probability=Rcpp::as<bool>(controls["update_global_probability"]); b.controls.update_ld_swap=Rcpp::as<bool>(controls["update_ld_swap"]); b.controls.ld_swap_probability=Rcpp::as<double>(controls["ld_swap_probability"]); b.controls.ld_swap_r2=Rcpp::as<double>(controls["ld_swap_r2"]); b.controls.ld_swap_max_friends=Rcpp::as<int>(controls["ld_swap_max_friends"]); b.controls.ld_swap_moves=Rcpp::as<int>(controls["ld_swap_moves"]);
 return b;
}
Rcpp::List phase9a_result(Rcpp::List spec, const char* policy) { return Rcpp::List::create(Rcpp::Named("schema")="blr_annotation_bayesc_contract_v1",Rcpp::Named("policy")=policy,Rcpp::Named("spec")=spec,Rcpp::Named("validated")=true,Rcpp::Named("invokes_sampler")=false); }
}

// [[Rcpp::export]]
Rcpp::List blr_phase9a_validate_fixed_prior_bayesc_cpp(Rcpp::List spec) {
 using namespace sblr::core; Rcpp::List data=spec["data"], controls=spec["controls"], p=spec["policy"];
 Phase9ACommonBuffers b=phase9a_common(data,controls); std::vector<double> pi=Rcpp::as<std::vector<double>>(p["marker_probability"]), mult=Rcpp::as<std::vector<double>>(p["marker_multiplier"]);
 FixedPriorBayesCPolicyView v; v.marker_probability=pi.data(); v.probability_count=pi.size(); v.marker_multiplier=mult.data(); v.multiplier_count=mult.size(); v.use_marker_probability=Rcpp::as<bool>(p["use_marker_probability"]); v.use_marker_multiplier=Rcpp::as<bool>(p["use_marker_multiplier"]); v.shared_read_only=Rcpp::as<bool>(p["shared_read_only"]); v.per_chain_payload=Rcpp::as<bool>(p["per_chain_payload"]); v.storage_outlives_execution=Rcpp::as<bool>(p["storage_outlives_execution"]);
 validate_annotation_bayesc_common(b.data,b.controls); validate_fixed_prior_policy(v,b.data.marker_count,b.data.trait_count); return phase9a_result(spec,"fixed_marker_prior");
}

// [[Rcpp::export]]
Rcpp::List blr_phase9a_validate_group_bayesc_cpp(Rcpp::List spec) {
 using namespace sblr::core; Rcpp::List data=spec["data"], controls=spec["controls"], p=spec["policy"];
 Phase9ACommonBuffers b=phase9a_common(data,controls); std::vector<int> map=Rcpp::as<std::vector<int>>(p["marker_group"]); std::vector<double> pi=Rcpp::as<std::vector<double>>(p["initial_probability"]), mult=Rcpp::as<std::vector<double>>(p["initial_multiplier"]), pa=Rcpp::as<std::vector<double>>(p["prior_a"]), pb=Rcpp::as<std::vector<double>>(p["prior_b"]);
 GroupBayesCPolicyView v; v.marker_group=map.data(); v.marker_group_count=map.size(); v.group_count=Rcpp::as<std::size_t>(p["group_count"]); v.group_order=Rcpp::as<std::vector<std::string>>(p["group_order"]); v.initial_probability=pi.data(); v.probability_count=pi.size(); v.initial_multiplier=mult.data(); v.multiplier_count=mult.size(); v.probability_prior_a=pa.data(); v.probability_prior_b=pb.data(); v.probability_prior_count=pa.size(); v.update_probability=Rcpp::as<bool>(p["update_probability"]); v.update_multiplier=Rcpp::as<bool>(p["update_multiplier"]); v.normalize_multiplier=Rcpp::as<bool>(p["normalize_multiplier"]); v.zero_based_index=Rcpp::as<bool>(p["zero_based_index"]); v.shared_read_only=Rcpp::as<bool>(p["shared_read_only"]); v.per_chain_payload=Rcpp::as<bool>(p["per_chain_payload"]); v.storage_outlives_execution=Rcpp::as<bool>(p["storage_outlives_execution"]);
 validate_annotation_bayesc_common(b.data,b.controls); validate_group_policy(v,b.data.marker_count,b.data.trait_count); return phase9a_result(spec,"group");
}

// [[Rcpp::export]]
Rcpp::List blr_phase9a_validate_learned_annotation_bayesc_cpp(Rcpp::List spec) {
 using namespace sblr::core; Rcpp::List data=spec["data"], controls=spec["controls"], p=spec["policy"]; Phase9ACommonBuffers b=phase9a_common(data,controls);
 Rcpp::NumericMatrix A=p["annotation"]; std::vector<double> ep=Rcpp::as<std::vector<double>>(p["eta_probability_init"]), ev=Rcpp::as<std::vector<double>>(p["eta_multiplier_init"]);
 LearnedAnnotationBayesCPolicyView v; v.marker_count=Rcpp::as<std::size_t>(p["marker_count"]); v.annotation_count=Rcpp::as<std::size_t>(p["annotation_count"]); v.trait_count=Rcpp::as<std::size_t>(p["trait_count"]); v.annotation=A.begin(); v.annotation_value_count=A.size(); v.annotation_order=Rcpp::as<std::vector<std::string>>(p["annotation_order"]); v.layout=Rcpp::as<std::string>(p["layout"]); v.includes_intercept=Rcpp::as<bool>(p["includes_intercept"]); v.eta_probability_init=ep.data(); v.eta_probability_count=ep.size(); v.eta_multiplier_init=ev.data(); v.eta_multiplier_count=ev.size(); v.learn_probability=Rcpp::as<bool>(p["learn_probability"]); v.learn_multiplier=Rcpp::as<bool>(p["learn_multiplier"]); v.probability_prior_sd=Rcpp::as<double>(p["probability_prior_sd"]); v.multiplier_prior_sd=Rcpp::as<double>(p["multiplier_prior_sd"]); v.probability_proposal_sd=Rcpp::as<double>(p["probability_proposal_sd"]); v.multiplier_proposal_sd=Rcpp::as<double>(p["multiplier_proposal_sd"]); v.update_every=Rcpp::as<int>(p["update_every"]); v.probability_min=Rcpp::as<double>(p["probability_min"]); v.probability_max=Rcpp::as<double>(p["probability_max"]); v.multiplier_min=Rcpp::as<double>(p["multiplier_min"]); v.multiplier_max=Rcpp::as<double>(p["multiplier_max"]); v.probability_link=Rcpp::as<std::string>(p["probability_link"]); v.multiplier_link=Rcpp::as<std::string>(p["multiplier_link"]); v.shared_read_only=Rcpp::as<bool>(p["shared_read_only"]); v.per_chain_payload=Rcpp::as<bool>(p["per_chain_payload"]); v.storage_outlives_execution=Rcpp::as<bool>(p["storage_outlives_execution"]);
 validate_annotation_bayesc_common(b.data,b.controls); validate_learned_annotation_policy(v,b.data.marker_count,b.data.trait_count); return phase9a_result(spec,"learned_annotation");
}

// Internal Phase 10A validation-only scheduled contract bridge. It constructs
// metadata/state vocabulary only and invokes neither sampler nor RNG.
// [[Rcpp::export]]
Rcpp::List blr_phase10a_validate_scheduled_execution_cpp(Rcpp::List spec) {
 using namespace sblr::core;
 Rcpp::List execution=spec["execution"], sweep=spec["sweep"], skip=spec["skip"];
 Rcpp::List candidate=spec["candidate"], neighbor=spec["neighbor"], state=spec["state"];
 ScheduledExecutionControl x;
 x.marker_count=Rcpp::as<std::size_t>(execution["marker_count"]);
 x.trait_count=Rcpp::as<std::size_t>(execution["trait_count"]);
 x.iterations=Rcpp::as<int>(execution["iterations"]); x.burnin=Rcpp::as<int>(execution["burnin"]);
 x.thinning=Rcpp::as<int>(execution["thinning"]); x.chains=Rcpp::as<int>(execution["chains"]);
 x.cores=Rcpp::as<int>(execution["cores"]); x.seed=Rcpp::as<int>(execution["seed"]);
 x.chain_seeds=Rcpp::as<std::vector<int>>(execution["chain_seeds"]);
 x.keep_chains=Rcpp::as<bool>(execution["keep_chains"]);
 Rcpp::List rng_ownership=spec["rng_ownership"];
 x.rng_ownership.engine_owner=Rcpp::as<std::string>(rng_ownership["engine_owner"]);
 x.rng_ownership.distribution_owner=Rcpp::as<std::string>(rng_ownership["distribution_owner"]);
 x.rng_ownership.lifetime=Rcpp::as<std::string>(rng_ownership["lifetime"]);
 x.rng_ownership.worker_thread_owner=Rcpp::as<std::string>(rng_ownership["worker_thread_owner"]);
 x.rng_ownership.fit_persistent_distribution_state=
  Rcpp::as<bool>(rng_ownership["fit_persistent_distribution_state"]);
 x.sweep.full_sweep_every=Rcpp::as<int>(sweep["full_sweep_every"]);
 x.sweep.iteration_zero_is_full=Rcpp::as<bool>(sweep["iteration_zero_is_full"]);
 x.skip.base_interval=Rcpp::as<int>(skip["null_skip_base"]);
 x.skip.maximum_interval=Rcpp::as<int>(skip["null_skip_max"]);
 x.skip.burnin_only=Rcpp::as<bool>(skip["burnin_only"]);
 x.skip.growth_rule=Rcpp::as<std::string>(skip["growth_rule"]);
 x.candidate.probability_threshold=Rcpp::as<double>(candidate["threshold"]);
 x.candidate.lifetime=Rcpp::as<int>(candidate["lifetime"]);
 x.neighbor.enabled=Rcpp::as<bool>(neighbor["enabled"]);
 x.neighbor.effect_difference_threshold=Rcpp::as<double>(neighbor["difference_threshold"]);
 x.neighbor.maximum_neighbors=Rcpp::as<int>(neighbor["maximum_neighbors"]);
 x.neighbor.friend_marker_count=Rcpp::as<std::size_t>(neighbor["friend_marker_count"]);
 x.neighbor.shared_read_only=Rcpp::as<bool>(neighbor["shared_read_only"]);
 x.neighbor.storage_outlives_execution=Rcpp::as<bool>(neighbor["storage_outlives_execution"]);
 // The non-null sentinel represents borrowed friend storage during validation.
 const int friend_sentinel=0; x.neighbor.friend_data=&friend_sentinel;
 validate_scheduled_execution_control(x);
 ScheduledMarkerState s;
 s.scheduled_at=Rcpp::as<std::vector<int>>(state["scheduled_at"]);
 s.last_updated=Rcpp::as<std::vector<int>>(state["last_updated"]);
 s.candidate=Rcpp::as<std::vector<unsigned char>>(state["candidate"]);
 s.in_candidate_list=Rcpp::as<std::vector<unsigned char>>(state["in_candidate_list"]);
 s.in_active_list=Rcpp::as<std::vector<unsigned char>>(state["in_active_list"]);
 s.last_interesting=Rcpp::as<std::vector<int>>(state["last_interesting"]);
 validate_scheduled_marker_state(s,x.marker_count);
 return Rcpp::List::create(Rcpp::Named("schema")="blr_scheduled_execution_contract_v1",
  Rcpp::Named("spec")=spec,Rcpp::Named("validated")=true,
  Rcpp::Named("invokes_sampler")=false,Rcpp::Named("consumes_rng")=false);
}

namespace {
double phase10a_persistent_normal(std::mt19937& engine, bool reset) {
 static thread_local std::normal_distribution<double> distribution(0.0,1.0);
 if (reset) distribution.reset();
 return distribution(engine);
}
}

// Internal Phase 10A diagnostic for the exact persistent-distribution pattern.
// [[Rcpp::export]]
Rcpp::List blr_phase10a_distribution_cache_diagnostic_cpp(int seed, int threads=2) {
 if (threads<=0) throw std::invalid_argument("threads must be positive");
 Rcpp::NumericVector cached(threads), fresh(threads), first(threads);
 std::vector<double> first_native(static_cast<std::size_t>(threads));
 std::vector<double> cached_native(static_cast<std::size_t>(threads));
 std::vector<double> fresh_native(static_cast<std::size_t>(threads));
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
 for (int thread=0; thread<threads; ++thread) {
  const unsigned int thread_seed=static_cast<unsigned int>(seed+1009*thread);
  std::mt19937 engine(thread_seed);
  first_native[static_cast<std::size_t>(thread)]=phase10a_persistent_normal(engine,true);
  std::mt19937 reseeded(thread_seed);
  cached_native[static_cast<std::size_t>(thread)]=phase10a_persistent_normal(reseeded,false);
  std::mt19937 reset_engine(thread_seed);
  fresh_native[static_cast<std::size_t>(thread)]=phase10a_persistent_normal(reset_engine,true);
 }
 for (int thread=0;thread<threads;++thread) {
  first[thread]=first_native[static_cast<std::size_t>(thread)];
  cached[thread]=cached_native[static_cast<std::size_t>(thread)];
  fresh[thread]=fresh_native[static_cast<std::size_t>(thread)];
 }
 bool cached_state_survives_engine_reseed=false;
 for (int thread=0;thread<threads;++thread) {
  if (cached_native[static_cast<std::size_t>(thread)] !=
      fresh_native[static_cast<std::size_t>(thread)]) {
   cached_state_survives_engine_reseed=true;
   break;
  }
 }
 return Rcpp::List::create(Rcpp::Named("first")=first,
  Rcpp::Named("after_engine_reseed_without_distribution_reset")=cached,
  Rcpp::Named("after_distribution_reset")=fresh,
  Rcpp::Named("cached_state_survives_engine_reseed")=
   cached_state_survives_engine_reseed,
  Rcpp::Named("threads")=threads);
}

// Internal Phase 10B diagnostic for fit-bounded chain RNG ownership.
// [[Rcpp::export]]
Rcpp::List blr_phase10b_chain_rng_diagnostic_cpp(int seed, int draws=7) {
 if (draws<=0) throw std::invalid_argument("draws must be positive");
 using sblr::core::ScheduledChainRng;
 ScheduledChainRng first(static_cast<std::uint64_t>(seed));
 ScheduledChainRng second(static_cast<std::uint64_t>(seed));
 Rcpp::NumericVector a(draws),b(draws),after_other(draws);
 for(int i=0;i<draws;++i) {
  a[i]=first.normal(first.engine);
  b[i]=second.normal(second.engine);
 }
 ScheduledChainRng other(static_cast<std::uint64_t>(seed+1));
 for(int i=0;i<2*draws+1;++i) (void)other.normal(other.engine);
 ScheduledChainRng reconstructed(static_cast<std::uint64_t>(seed));
 for(int i=0;i<draws;++i) after_other[i]=reconstructed.normal(reconstructed.engine);
 return Rcpp::List::create(
  Rcpp::Named("first")=a,Rcpp::Named("identical_seed")=b,
  Rcpp::Named("after_odd_other_chain")=after_other,
  Rcpp::Named("owner")="chain",Rcpp::Named("lifetime")="one_chain_execution",
  Rcpp::Named("worker_thread_owner")="none");
}
