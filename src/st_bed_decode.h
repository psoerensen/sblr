#ifndef ST_BED_DECODE_H
#define ST_BED_DECODE_H

// Shared standardized BED block-decode helpers.

#include <cmath>
#include <cstdint>
#include <vector>

#include "packed_bed.h"

#ifdef _OPENMP
#include <omp.h>
#endif

inline int effective_nthreads(int nthreads) {
#ifdef _OPENMP
  if (nthreads > 0) return nthreads;
  return omp_get_max_threads();
#else
  (void)nthreads;
  return 1;
#endif
}

inline void standardized_lut(double p, float lut[4]) {
  const double denom = std::sqrt(2.0 * p * (1.0 - p));
  if (denom <= 0.0 || !std::isfinite(denom)) {
    lut[0] = lut[1] = lut[2] = lut[3] = 0.0f;
  } else {
    lut[0] = static_cast<float>((2.0 - 2.0 * p) / denom);  // BED 00 -> dosage 2
    lut[1] = 0.0f;                                          // BED 01 -> missing
    lut[2] = static_cast<float>((1.0 - 2.0 * p) / denom);  // BED 10 -> dosage 1
    lut[3] = static_cast<float>((0.0 - 2.0 * p) / denom);  // BED 11 -> dosage 0
  }
}

// Decode markers [marker_start, marker_start+marker_len) into standardized
// float columns Z (row-major per marker: Z + ii*n), length marker_len * n.
inline void decode_packed_block_float(
    const PackedBedMatrix& G, int marker_start, int marker_len,
    const double* af, float* Z, int nthreads) {
  const int n = G.n;
  const std::size_t nbytes = G.nbytes;
  const int nth = effective_nthreads(nthreads);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
  for (int ii = 0; ii < marker_len; ++ii) {
    const int global_i = marker_start + ii;
    const uint8_t* packed = G.row(global_i);
    float* z = Z + static_cast<std::size_t>(ii) * n;
    float lut[4];
    standardized_lut(af[global_i], lut);
    for (std::size_t kb = 0; kb < nbytes; ++kb) {
      const uint8_t x = packed[kb];
      const int jbase = static_cast<int>(kb << 2);
      if (jbase + 0 < n) z[jbase + 0] = lut[(x >> 0) & 3u];
      if (jbase + 1 < n) z[jbase + 1] = lut[(x >> 2) & 3u];
      if (jbase + 2 < n) z[jbase + 2] = lut[(x >> 4) & 3u];
      if (jbase + 3 < n) z[jbase + 3] = lut[(x >> 6) & 3u];
    }
  }
}

// x_i'x_i for every marker under the same standardization.
inline std::vector<double> compute_xx_from_packed_standardized(
    const PackedBedMatrix& G, const double* af, int nthreads) {
  const int n = G.n;
  const int m = G.m;
  const std::size_t nbytes = G.nbytes;
  const int nth = effective_nthreads(nthreads);
  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
  for (int marker = 0; marker < m; ++marker) {
    const uint8_t* packed = G.row(marker);
    float lut[4];
    standardized_lut(af[marker], lut);
    double s = 0.0;
    for (std::size_t kb = 0; kb < nbytes; ++kb) {
      const uint8_t x = packed[kb];
      const int jbase = static_cast<int>(kb << 2);
      if (jbase + 0 < n) { const double z = lut[(x >> 0) & 3u]; s += z * z; }
      if (jbase + 1 < n) { const double z = lut[(x >> 2) & 3u]; s += z * z; }
      if (jbase + 2 < n) { const double z = lut[(x >> 4) & 3u]; s += z * z; }
      if (jbase + 3 < n) { const double z = lut[(x >> 6) & 3u]; s += z * z; }
    }
    xx[static_cast<std::size_t>(marker)] = s;
  }
  return xx;
}

#endif  // ST_BED_DECODE_H
