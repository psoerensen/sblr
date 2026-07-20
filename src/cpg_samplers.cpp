#include "cpg_samplers.h"

#include <algorithm>
#include <numeric>

// Shared by the retained internal mtblr_cpg_omp_csr() research route.
void samplePi_cpg(std::vector<double>& cmodel,
                  std::vector<double>& pi,
                  std::mt19937& gen) {
 for (size_t k = 0; k < cmodel.size(); k++) {
  std::gamma_distribution<double> rgamma(cmodel[k], 1.0);
  pi[k] = rgamma(gen);
 }
 double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
 for (size_t k = 0; k < cmodel.size(); k++) {
  pi[k] = pi[k] / psum;
 }
 std::fill(cmodel.begin(), cmodel.end(), 1.0);
}
