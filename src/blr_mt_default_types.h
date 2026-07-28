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
};

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
