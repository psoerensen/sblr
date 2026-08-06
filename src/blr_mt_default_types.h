#ifndef SBLR_BLR_MT_DEFAULT_TYPES_H
#define SBLR_BLR_MT_DEFAULT_TYPES_H

#include <armadillo>
#include "blr_mt_bayesr_types.h"
#include "blr_mt_bayesrc_types.h"
#include <utility>
#include <vector>

namespace sblr {
namespace mt {

struct MtDefaultDataView {
 const std::vector<std::vector<double>>& wy;
 const std::vector<std::vector<double>>& ww;
 const std::vector<double>& yy;
 const std::vector<std::vector<std::vector<double>>>& XXvalues;
 const std::vector<std::vector<int>>& XXindices;
 const std::vector<int>& n;
};

struct MtExtendedTraceSpec {
 bool covariance=false;
 bool probability=false;
 bool annotations=false;
 bool full_probability_states=false;
 bool selected_b=false;
 bool selected_d=false;
 bool selected_component=false;
 bool aggregate_component_states=false;
 std::vector<int> selected_markers;
};

struct MtExtendedTraceResult {
 std::vector<std::vector<double>> cov_b;
 std::vector<std::vector<double>> cov_g;
 std::vector<std::vector<double>> cov_e;
 std::vector<std::vector<double>> component_pi;
 std::vector<std::vector<double>> pattern_pi;
 std::vector<std::vector<double>> joint_pi;
 std::vector<std::vector<double>> annotation_alpha;
 std::vector<std::vector<double>> annotation_sigma;
 std::vector<std::vector<double>> selected_b;
 std::vector<std::vector<int>> selected_d;
 std::vector<std::vector<int>> selected_component;
 std::vector<std::vector<int>> component_count;
 std::vector<std::vector<int>> realized_active_count;
 std::vector<std::vector<int>> stick_eligible_count;
 std::vector<std::vector<int>> stick_continue_count;
 std::vector<std::vector<int>> stick_stop_count;
};

inline std::vector<int> mt_component_counts(
 const std::vector<int>& component, int component_count
) {
 if (component_count<2)
  throw std::invalid_argument("aggregate component tracing requires at least two components");
 std::vector<int> counts(static_cast<std::size_t>(component_count),0);
 for (int value:component) {
  if (value<0 || value>=component_count)
   throw std::invalid_argument("component state is outside the aggregate trace domain");
  ++counts[static_cast<std::size_t>(value)];
 }
 return counts;
}

inline void mt_allocate_aggregate_component_traces(
 MtExtendedTraceResult& trace, int component_count, int iterations
) {
 if (component_count<2 || iterations<1)
  throw std::invalid_argument("invalid aggregate component trace dimensions");
 trace.component_count.assign(static_cast<std::size_t>(component_count),
  std::vector<int>(static_cast<std::size_t>(iterations)));
 trace.realized_active_count.assign(1,
  std::vector<int>(static_cast<std::size_t>(iterations)));
 const std::size_t sticks=static_cast<std::size_t>(component_count-1);
 trace.stick_eligible_count.assign(sticks,
  std::vector<int>(static_cast<std::size_t>(iterations)));
 trace.stick_continue_count.assign(sticks,
  std::vector<int>(static_cast<std::size_t>(iterations)));
 trace.stick_stop_count.assign(sticks,
  std::vector<int>(static_cast<std::size_t>(iterations)));
}

inline void mt_capture_aggregate_component_state(
 const std::vector<int>& counts, int marker_count, int iteration,
 MtExtendedTraceResult& trace
) {
 if (counts.size()<2 || iteration<0 || marker_count<0 ||
     trace.component_count.size()!=counts.size())
  throw std::invalid_argument("invalid aggregate component capture state");
 int total=0;
 for (std::size_t component=0;component<counts.size();++component) {
  const int count=counts[component];
  if (count<0) throw std::invalid_argument("component count cannot be negative");
  total+=count;
  trace.component_count[component][static_cast<std::size_t>(iteration)]=count;
 }
 if (total!=marker_count)
  throw std::invalid_argument("aggregate component counts do not sum to marker count");
 const int active=marker_count-counts[0];
 trace.realized_active_count[0][static_cast<std::size_t>(iteration)]=active;
 int eligible=marker_count;
 for (std::size_t stick=0;stick+1<counts.size();++stick) {
  const int stop=counts[stick];
  const int continued=eligible-stop;
  trace.stick_eligible_count[stick][static_cast<std::size_t>(iteration)]=eligible;
  trace.stick_continue_count[stick][static_cast<std::size_t>(iteration)]=continued;
  trace.stick_stop_count[stick][static_cast<std::size_t>(iteration)]=stop;
  eligible=continued;
 }
}

struct MtDefaultModelSpec {
 const std::vector<std::vector<int>>& models;
 const std::vector<std::vector<int>>& sets;
 int method;
 const MtJointStateSpec* joint=nullptr;
 const std::vector<double>* marker_scale=nullptr;
 const std::vector<double>* pi_prior=nullptr;
 const MtBayesRCSpec* bayesrc=nullptr;
 const MtExtendedTraceSpec* convergence=nullptr;
};

struct MtDefaultCovariancePriorView {
 const std::vector<std::vector<double>>& ssb_prior;
 const std::vector<std::vector<double>>& sse_prior;
 double nub;
 double nue;
};

struct MtDefaultExecutionSpec {
 bool updateB;
 bool updateE;
 bool updatePi;
 int nit;
 int nburn;
 int nthin;
 int seed;
};

struct MtDefaultInitialState {
 std::vector<std::vector<double>> b;
 arma::mat B;
 arma::mat E;
 std::vector<double> pi;
 std::vector<int> component;
 std::vector<std::vector<double>> beta;
 std::vector<std::vector<int>> state;
 arma::mat annotation_alpha;
 arma::vec annotation_sigma;
 std::vector<double> pattern_pi;
};

struct MtDefaultCoreResult {
 int nt=0;
 int m=0;
 int nmodels=0;

 double marker_retained_count=0.0;
 double covb_retained_count=0.0;
 double covg_retained_count=0.0;
 double cove_retained_count=0.0;
 double pi_retained_count=0.0;

 std::vector<std::vector<double>> bm;
 std::vector<std::vector<double>> dm;
 std::vector<std::vector<double>> r;
 std::vector<std::vector<double>> b;
 std::vector<std::vector<int>> d;
 std::vector<int> component;
 std::vector<std::vector<double>> component_counts;
 std::vector<int> order;

 std::vector<std::vector<double>> vbs;
 std::vector<std::vector<double>> vgs;
 std::vector<std::vector<double>> ves;
 std::vector<std::vector<double>> vle;
 std::vector<std::vector<double>> vld;

 std::vector<std::vector<double>> cvbm;
 std::vector<std::vector<double>> cvgm;
 std::vector<std::vector<double>> cvem;

 arma::mat B;
 arma::mat G;
 arma::mat E;
 std::vector<double> pi;
 std::vector<double> pis;
 std::vector<std::vector<double>> pi_trace;

 std::vector<std::vector<double>> pistrait;
 std::vector<double> pismarker;

 arma::mat annotation_alpha_final;
 arma::mat annotation_alpha_mean;
 arma::vec annotation_sigma_final;
 arma::vec annotation_sigma_mean;
 std::vector<double> pattern_pi_final;
 std::vector<double> pattern_pi_mean;
 std::vector<std::vector<double>> pattern_pi_trace;
 arma::mat prior_component_probabilities;
 int annotation_updates_attempted=0;
 int annotation_updates_completed=0;
 MtExtendedTraceResult convergence;
};

struct MtDefaultFinalResult {
 int nt=0;
 int m=0;
 int nmodels=0;

 double marker_retained_count=0.0;
 double covb_retained_count=0.0;
 double covg_retained_count=0.0;
 double cove_retained_count=0.0;
 double pi_retained_count=0.0;

 std::vector<std::vector<double>> bm;
 std::vector<std::vector<double>> dm;
 std::vector<std::vector<double>> r;
 std::vector<std::vector<double>> b;
 std::vector<std::vector<int>> d;
 std::vector<int> component;
 std::vector<std::vector<double>> component_probabilities;
 std::vector<int> marker_order;

 std::vector<std::vector<double>> vbs;
 std::vector<std::vector<double>> vgs;
 std::vector<std::vector<double>> ves;
 std::vector<std::vector<double>> vle;
 std::vector<std::vector<double>> vld;

 std::vector<std::vector<double>> covb;
 std::vector<std::vector<double>> covg;
 std::vector<std::vector<double>> cove;

 arma::mat vb;
 arma::mat vg;
 arma::mat ve;
 std::vector<double> pi_final;
 std::vector<double> pi_mean;
 std::vector<std::vector<double>> pi_trace;

 std::vector<std::vector<double>> pitrait;
 std::vector<double> pimarker;

 arma::mat annotation_alpha_final;
 arma::mat annotation_alpha_mean;
 arma::vec annotation_sigma_final;
 arma::vec annotation_sigma_mean;
 std::vector<double> pattern_pi_final;
 std::vector<double> pattern_pi_mean;
 std::vector<std::vector<double>> pattern_pi_trace;
 arma::mat prior_component_probabilities;
 int annotation_updates_attempted=0;
 int annotation_updates_completed=0;
 MtExtendedTraceResult convergence;
};

}  // namespace mt
}  // namespace sblr

#endif
