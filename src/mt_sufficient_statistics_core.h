// =============================================================================
// mt_sufficient_statistics_core.h
// =============================================================================

#pragma once

#include <string>
#include <vector>

struct BedXtXytyResult {
 int n = 0;        // number of individuals used for y / sufficient statistics
 int n_bed = 0;    // number of individuals in the BED file
 int m = 0;
 int nt = 0;
 bool scale = true;
 bool rows_used = false;

 // Trait-major layout: wy[t][i], ww[t][i]
 std::vector<std::vector<double>> wy;
 std::vector<std::vector<double>> ww;

 // yy[t] = y_t' y_t over selected rows
 std::vector<double> yy;

 // xx[i] = x_i' x_i over selected rows. Same values are repeated in ww[t][i].
 std::vector<double> xx;

 // Allele frequencies used for standardization.
 // If af_computed = true, these are computed over the selected rows.
 std::vector<double> af;
 bool af_computed = false;

 // 0-based selected BED rows used by the core. Empty means all rows in order.
 std::vector<int> rows0;

 double mean_xx = 0.0;
 double min_xx = 0.0;
 double max_xx = 0.0;
};

void bed_xtx_diag_xty_core(
  const char* file,
  int n_bed,
  const int* cls,
  const double* af,
  bool use_af,
  double* af_used,
  const int* rows0,
  int n_used,
  int m,
  int nt,
  bool scale,
  const double* y_flat,   // n_used x nt, row-major: y_flat[k * nt + t]
  double* xty_flat,       // m x nt, row-major: xty_flat[i * nt + t]
  double* xx,             // length m
  int nthreads,
  int MG = 64,
  int JB = 1024,
  int TB = 32
);

BedXtXytyResult bed_xtx_xty_cpp(
  const std::string& bed_file,
  int n_bed,
  const std::vector<int>& cls,
  const std::vector<double>& af,
  bool use_af,
  const std::vector<int>& rows0,
  bool use_rows,
  const std::vector<double>& y_flat,
  int nt,
  bool scale,
  int nthreads,
  int MG = 64,
  int JB = 1024,
  int TB = 32
);

