#pragma once

#include <random>
#include <vector>

// Shared by the retained internal mtblr_cpg_omp_csr() research route.
void samplePi_cpg(
  std::vector<double>& cmodel,
  std::vector<double>& pi,
  std::mt19937& gen);
