#include <Rcpp.h>

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "packed_bed.h"

static std::vector<std::string> rcpp_to_string_vector(
  Rcpp::CharacterVector x
) {
 std::vector<std::string> out(static_cast<std::size_t>(x.size()));

 for (int i = 0; i < x.size(); ++i) {
  out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(x[i]);
 }

 return out;
}

static std::vector<int> rcpp_rows0_or_empty(
  Rcpp::Nullable<Rcpp::IntegerVector> rows
) {
 std::vector<int> out;

 if (rows.isNotNull()) {
  Rcpp::IntegerVector r(rows);
  out.resize(static_cast<std::size_t>(r.size()));

  for (int i = 0; i < r.size(); ++i) {
   if (r[i] == NA_INTEGER) {
    throw std::runtime_error("rows contains NA.");
   }

   if (r[i] <= 0) {
    throw std::runtime_error("rows must contain positive 1-based indices.");
   }

   out[static_cast<std::size_t>(i)] = r[i] - 1;
  }
 }

 return out;
}

static std::vector<std::vector<int>> rcpp_to_int_list(
  Rcpp::List cls
) {
 std::vector<std::vector<int>> out(static_cast<std::size_t>(cls.size()));

 for (int f = 0; f < cls.size(); ++f) {
  Rcpp::IntegerVector x = cls[f];
  out[static_cast<std::size_t>(f)].resize(static_cast<std::size_t>(x.size()));

  for (int k = 0; k < x.size(); ++k) {
   if (x[k] == NA_INTEGER) {
    throw std::runtime_error("cls contains NA.");
   }

   if (x[k] <= 0) {
    throw std::runtime_error("cls must contain positive 1-based marker indices.");
   }

   out[static_cast<std::size_t>(f)][static_cast<std::size_t>(k)] = x[k];
  }
 }

 return out;
}

static std::vector<double> rcpp_flatten_numeric_list_or_empty(
  Rcpp::Nullable<Rcpp::List> xlist
) {
 std::vector<double> out;

 if (xlist.isNotNull()) {
  Rcpp::List xl = Rcpp::as<Rcpp::List>(xlist.get());

  for (int f = 0; f < xl.size(); ++f) {
   Rcpp::NumericVector x = xl[f];

   for (int i = 0; i < x.size(); ++i) {
    out.push_back(x[i]);
   }
  }
 }

 return out;
}

static inline double decode_standardized_or_raw_bed_code(
  unsigned int code,
  double p,
  bool scale
) {
 if (scale) {
  const double denom = std::sqrt(2.0 * p * (1.0 - p));

  if (denom <= 0.0 || !std::isfinite(denom)) {
   return 0.0;
  }

  switch (code) {
  case 0u: return (2.0 - 2.0 * p) / denom;  // BED 00 -> dosage 2
  case 1u: return 0.0;                      // BED 01 -> missing, mean-imputed after standardization
  case 2u: return (1.0 - 2.0 * p) / denom;  // BED 10 -> dosage 1
  default: return (0.0 - 2.0 * p) / denom;  // BED 11 -> dosage 0
  }
 }

 switch (code) {
 case 0u: return 2.0;
 case 1u: return 2.0 * p;  // missing -> mean dosage
 case 2u: return 1.0;
 default: return 0.0;
 }
}

// [[Rcpp::export]]
Rcpp::List test_read_bedfiles_to_packed_matrix(
  Rcpp::CharacterVector bedfiles,
  int n,
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::List cls
) {
 std::vector<std::string> bedfiles_cpp =
  rcpp_to_string_vector(bedfiles);

 std::vector<int> rows0 =
  rcpp_rows0_or_empty(rows);

 std::vector<std::vector<int>> cls_cpp =
  rcpp_to_int_list(cls);

 PackedBedMatrix G = read_bedfiles_to_packed_matrix(
  bedfiles_cpp,
  n,
  rows0.empty() ? nullptr : rows0.data(),
  rows0.empty() ? 0 : static_cast<int>(rows0.size()),
  cls_cpp
 );

 return Rcpp::List::create(
  Rcpp::Named("n_bed") = n,
  Rcpp::Named("n_used") = G.n,
  Rcpp::Named("m") = G.m,
  Rcpp::Named("nbytes_per_variant") = static_cast<double>(G.nbytes),
  Rcpp::Named("stride") = static_cast<double>(G.stride),
  Rcpp::Named("packed_bytes") =
   static_cast<double>(static_cast<std::size_t>(G.m) * G.nbytes),
   Rcpp::Named("packed_stride_bytes") =
    static_cast<double>(static_cast<std::size_t>(G.m) * G.stride)
 );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix test_read_bedfiles_to_dense_matrix(
  Rcpp::CharacterVector bedfiles,
  int n,
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::List cls
) {
 std::vector<std::string> bedfiles_cpp =
  rcpp_to_string_vector(bedfiles);

 std::vector<int> rows0 =
  rcpp_rows0_or_empty(rows);

 std::vector<std::vector<int>> cls_cpp =
  rcpp_to_int_list(cls);

 PackedBedMatrix G = read_bedfiles_to_packed_matrix(
  bedfiles_cpp,
  n,
  rows0.empty() ? nullptr : rows0.data(),
  rows0.empty() ? 0 : static_cast<int>(rows0.size()),
  cls_cpp
 );

 Rcpp::NumericMatrix X(G.n, G.m);

 for (int j = 0; j < G.m; ++j) {
  const uint8_t* row = G.row(j);

  for (int i = 0; i < G.n; ++i) {
   const double value = code_to_dosage(get_bed_code(row, i));

   if (std::isnan(value)) {
    X(i, j) = NA_REAL;
   } else {
    X(i, j) = value;
   }
  }
 }

 return X;
}

// [[Rcpp::export]]
Rcpp::DataFrame bed_pairwise_xtx_check(
  Rcpp::CharacterVector bed_files,
  int n,
  Rcpp::List cls,
  Rcpp::IntegerMatrix pairs,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::Nullable<Rcpp::List> af = R_NilValue,
  bool scale = true
) {
 if (n <= 0) {
  throw std::runtime_error("bed_pairwise_xtx_check: n must be positive.");
 }

 if (pairs.ncol() != 2) {
  throw std::runtime_error("bed_pairwise_xtx_check: pairs must have exactly two columns.");
 }

 std::vector<std::string> bed_files_cpp =
  rcpp_to_string_vector(bed_files);

 std::vector<std::vector<int>> cls_by_file =
  rcpp_to_int_list(cls);

 std::vector<int> rows0 =
  rcpp_rows0_or_empty(rows);

 const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
 const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());

 PackedBedMatrix G = read_bedfiles_to_packed_matrix(
  bed_files_cpp,
  n,
  rows0_ptr,
  n_rows,
  cls_by_file
 );

 const int n_used = G.n;
 const int m = G.m;

 std::vector<double> af_cpp =
  rcpp_flatten_numeric_list_or_empty(af);

 const bool af_computed = af_cpp.empty();

 if (af_computed) {
  af_cpp = compute_af_from_packed(G);
 }

 if (static_cast<int>(af_cpp.size()) != m) {
  throw std::runtime_error(
    "bed_pairwise_xtx_check: af must have one value per selected marker."
  );
 }

 for (int i = 0; i < m; ++i) {
  const double p = af_cpp[static_cast<std::size_t>(i)];

  if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
   throw std::runtime_error(
     "bed_pairwise_xtx_check: af contains invalid values; expected finite values in [0, 1]."
   );
  }
 }

 const int npairs = pairs.nrow();

 Rcpp::IntegerVector out_i(npairs);
 Rcpp::IntegerVector out_j(npairs);
 Rcpp::NumericVector out_xx_i(npairs);
 Rcpp::NumericVector out_xx_j(npairs);
 Rcpp::NumericVector out_xij(npairs);
 Rcpp::NumericVector out_r(npairs);

 for (int k = 0; k < npairs; ++k) {
  if (pairs(k, 0) == NA_INTEGER || pairs(k, 1) == NA_INTEGER) {
   throw std::runtime_error("bed_pairwise_xtx_check: pairs contains NA.");
  }

  const int i = pairs(k, 0) - 1;
  const int j = pairs(k, 1) - 1;

  if (i < 0 || i >= m || j < 0 || j >= m) {
   throw std::runtime_error(
     "bed_pairwise_xtx_check: pairs contains marker index out of range."
   );
  }

  const uint8_t* row_i = G.row(i);
  const uint8_t* row_j = G.row(j);

  const double p_i = af_cpp[static_cast<std::size_t>(i)];
  const double p_j = af_cpp[static_cast<std::size_t>(j)];

  double xx_i = 0.0;
  double xx_j = 0.0;
  double xij = 0.0;

  for (std::size_t kb = 0; kb < G.nbytes; ++kb) {
   const uint8_t byte_i = row_i[kb];
   const uint8_t byte_j = row_j[kb];
   const int jbase = static_cast<int>(kb << 2);

   for (int pos = 0; pos < 4; ++pos) {
    const int sample = jbase + pos;

    if (sample >= n_used) {
     break;
    }

    const unsigned int code_i = (byte_i >> (2 * pos)) & 3u;
    const unsigned int code_j = (byte_j >> (2 * pos)) & 3u;

    const double xi = decode_standardized_or_raw_bed_code(
     code_i,
     p_i,
     scale
    );

    const double xj = decode_standardized_or_raw_bed_code(
     code_j,
     p_j,
     scale
    );

    xx_i += xi * xi;
    xx_j += xj * xj;
    xij += xi * xj;
   }
  }

  const double denom = std::sqrt(xx_i * xx_j);
  const double rij = denom > 0.0 ? xij / denom : NA_REAL;

  out_i[k] = i + 1;
  out_j[k] = j + 1;
  out_xx_i[k] = xx_i;
  out_xx_j[k] = xx_j;
  out_xij[k] = xij;
  out_r[k] = rij;
 }

 return Rcpp::DataFrame::create(
  Rcpp::Named("i") = out_i,
  Rcpp::Named("j") = out_j,
  Rcpp::Named("xx_i_exact") = out_xx_i,
  Rcpp::Named("xx_j_exact") = out_xx_j,
  Rcpp::Named("xij_exact") = out_xij,
  Rcpp::Named("r_exact") = out_r,
  Rcpp::Named("n_used") = n_used,
  Rcpp::Named("af_computed") = af_computed
 );
}


// [[Rcpp::export]]
Rcpp::DataFrame bed_score_from_b_check(
  Rcpp::CharacterVector bed_files,
  int n,
  Rcpp::List cls,
  Rcpp::NumericVector y,
  Rcpp::NumericVector b,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::Nullable<Rcpp::List> af = R_NilValue,
  bool scale = true
) {
 std::vector<std::string> bed_files_cpp =
  rcpp_to_string_vector(bed_files);

 std::vector<std::vector<int>> cls_by_file =
  rcpp_to_int_list(cls);

 std::vector<int> rows0 =
  rcpp_rows0_or_empty(rows);

 const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
 const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());

 PackedBedMatrix G = read_bedfiles_to_packed_matrix(
  bed_files_cpp,
  n,
  rows0_ptr,
  n_rows,
  cls_by_file
 );

 const int n_used = G.n;
 const int m = G.m;

 if (y.size() != n_used) {
  throw std::runtime_error(
    "bed_score_from_b_check: length(y) must equal number of samples used."
  );
 }

 if (b.size() != m) {
  throw std::runtime_error(
    "bed_score_from_b_check: length(b) must equal number of selected markers."
  );
 }

 std::vector<double> af_cpp =
  rcpp_flatten_numeric_list_or_empty(af);

 const bool af_computed = af_cpp.empty();

 if (af_computed) {
  af_cpp = compute_af_from_packed(G);
 }

 if (static_cast<int>(af_cpp.size()) != m) {
  throw std::runtime_error(
    "bed_score_from_b_check: af must have one value per selected marker."
  );
 }

 auto decode_value = [&](unsigned int code, double p) -> double {
  if (scale) {
   const double denom = std::sqrt(2.0 * p * (1.0 - p));

   if (denom <= 0.0 || !std::isfinite(denom)) {
    return 0.0;
   }

   switch (code) {
   case 0u: return (2.0 - 2.0 * p) / denom;
   case 1u: return 0.0;
   case 2u: return (1.0 - 2.0 * p) / denom;
   default: return (0.0 - 2.0 * p) / denom;
   }
  }

  switch (code) {
  case 0u: return 2.0;
  case 1u: return 2.0 * p;
  case 2u: return 1.0;
  default: return 0.0;
  }
 };

 std::vector<double> xb(static_cast<std::size_t>(n_used), 0.0);

 for (int j = 0; j < m; ++j) {
  const double bj = b[j];

  if (bj == 0.0) {
   continue;
  }

  const uint8_t* row = G.row(j);
  const double p = af_cpp[static_cast<std::size_t>(j)];

  for (std::size_t kb = 0; kb < G.nbytes; ++kb) {
   const uint8_t byte = row[kb];
   const int jbase = static_cast<int>(kb << 2);

   for (int pos = 0; pos < 4; ++pos) {
    const int sample = jbase + pos;

    if (sample >= n_used) {
     break;
    }

    const unsigned int code = (byte >> (2 * pos)) & 3u;
    const double x = decode_value(code, p);

    xb[static_cast<std::size_t>(sample)] += x * bj;
   }
  }
 }

 std::vector<double> e(static_cast<std::size_t>(n_used), 0.0);

 for (int i = 0; i < n_used; ++i) {
  e[static_cast<std::size_t>(i)] =
   y[i] - xb[static_cast<std::size_t>(i)];
 }

 Rcpp::NumericVector r_exact(m);
 Rcpp::NumericVector wy_exact(m);
 Rcpp::NumericVector ww_exact(m);

 for (int j = 0; j < m; ++j) {
  const uint8_t* row = G.row(j);
  const double p = af_cpp[static_cast<std::size_t>(j)];

  double xte = 0.0;
  double xty = 0.0;
  double xtx = 0.0;

  for (std::size_t kb = 0; kb < G.nbytes; ++kb) {
   const uint8_t byte = row[kb];
   const int jbase = static_cast<int>(kb << 2);

   for (int pos = 0; pos < 4; ++pos) {
    const int sample = jbase + pos;

    if (sample >= n_used) {
     break;
    }

    const unsigned int code = (byte >> (2 * pos)) & 3u;
    const double x = decode_value(code, p);

    xte += x * e[static_cast<std::size_t>(sample)];
    xty += x * y[sample];
    xtx += x * x;
   }
  }

  r_exact[j] = xte;
  wy_exact[j] = xty;
  ww_exact[j] = xtx;
 }

 return Rcpp::DataFrame::create(
  Rcpp::Named("marker") = Rcpp::seq(1, m),
  Rcpp::Named("r_exact") = r_exact,
  Rcpp::Named("wy_exact") = wy_exact,
  Rcpp::Named("ww_exact") = ww_exact,
  Rcpp::Named("n_used") = n_used,
  Rcpp::Named("af_computed") = af_computed
 );
}


// #include <Rcpp.h>
// #include <vector>
// #include <string>
//
// #include "packed_bed.h"
//
// static std::vector<std::string> rcpp_to_string_vector(
//   Rcpp::CharacterVector x
// ) {
//  std::vector<std::string> out(x.size());
//  for (int i = 0; i < x.size(); ++i) {
//   out[i] = Rcpp::as<std::string>(x[i]);
//  }
//  return out;
// }
//
// static std::vector<int> rcpp_rows0_or_empty(
//   Rcpp::Nullable<Rcpp::IntegerVector> rows
// ) {
//  std::vector<int> out;
//
//  if (rows.isNotNull()) {
//   Rcpp::IntegerVector r(rows);
//   out.resize(r.size());
//
//   for (int i = 0; i < r.size(); ++i) {
//    out[i] = r[i] - 1;
//   }
//  }
//
//  return out;
// }
//
// static std::vector<std::vector<int>> rcpp_to_int_list(
//   Rcpp::List cls
// ) {
//  std::vector<std::vector<int>> out(cls.size());
//
//  for (int f = 0; f < cls.size(); ++f) {
//   Rcpp::IntegerVector x = cls[f];
//   out[f].resize(x.size());
//
//   for (int k = 0; k < x.size(); ++k) {
//    out[f][k] = x[k];
//   }
//  }
//
//  return out;
// }
//
// // [[Rcpp::export]]
// Rcpp::List test_read_bedfiles_to_packed_matrix(
//   Rcpp::CharacterVector bedfiles,
//   int n,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   Rcpp::List cls
// ) {
//  std::vector<std::string> bedfiles_cpp =
//   rcpp_to_string_vector(bedfiles);
//
//  std::vector<int> rows0 =
//   rcpp_rows0_or_empty(rows);
//
//  std::vector<std::vector<int>> cls_cpp =
//   rcpp_to_int_list(cls);
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bedfiles_cpp,
//   n,
//   rows0.empty() ? nullptr : rows0.data(),
//   rows0.empty() ? 0 : static_cast<int>(rows0.size()),
//   cls_cpp
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("n_bed") = n,
//   Rcpp::Named("n_used") = G.n,
//   Rcpp::Named("m") = G.m,
//   Rcpp::Named("nbytes_per_variant") = static_cast<double>(G.nbytes),
//   Rcpp::Named("stride") = static_cast<double>(G.stride),
//   Rcpp::Named("packed_bytes") =
//    static_cast<double>(static_cast<std::size_t>(G.m) * G.nbytes),
//    Rcpp::Named("packed_stride_bytes") =
//     static_cast<double>(static_cast<std::size_t>(G.m) * G.stride)
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::NumericMatrix test_read_bedfiles_to_dense_matrix(
//   Rcpp::CharacterVector bedfiles,
//   int n,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   Rcpp::List cls
// ) {
//  std::vector<std::string> bedfiles_cpp =
//   rcpp_to_string_vector(bedfiles);
//
//  std::vector<int> rows0 =
//   rcpp_rows0_or_empty(rows);
//
//  std::vector<std::vector<int>> cls_cpp =
//   rcpp_to_int_list(cls);
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bedfiles_cpp,
//   n,
//   rows0.empty() ? nullptr : rows0.data(),
//   rows0.empty() ? 0 : static_cast<int>(rows0.size()),
//   cls_cpp
//  );
//
//  Rcpp::NumericMatrix X(G.n, G.m);
//
//  for (int j = 0; j < G.m; ++j) {
//   const uint8_t* row = G.row(j);
//
//   for (int i = 0; i < G.n; ++i) {
//    const double value = code_to_dosage(get_bed_code(row, i));
//
//    if (std::isnan(value)) {
//     X(i, j) = NA_REAL;
//    } else {
//     X(i, j) = value;
//    }
//   }
//  }
//
//  return X;
// }
//
//
// // [[Rcpp::export]]
// Rcpp::DataFrame bed_pairwise_xtx_check(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::IntegerMatrix pairs,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   bool scale = true
// ) {
//  std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
//  std::vector<int> rows0 = copy_rows0_or_empty(rows, n);
//
//  const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
//  const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bed_files_cpp,
//   n,
//   rows0_ptr,
//   n_rows,
//   cls_by_file
//  );
//
//  const int n_used = G.n;
//  const int m = G.m;
//
//  std::vector<double> af_cpp = flatten_af_list(af);
//  const bool af_computed = af_cpp.empty();
//
//  if (af_computed) {
//   af_cpp = compute_af_from_packed(G);
//  }
//
//  if (static_cast<int>(af_cpp.size()) != m) {
//   throw std::runtime_error("bed_pairwise_xtx_check: af must have one value per selected marker.");
//  }
//
//  const int npairs = pairs.nrow();
//
//  Rcpp::IntegerVector out_i(npairs);
//  Rcpp::IntegerVector out_j(npairs);
//  Rcpp::NumericVector out_xx_i(npairs);
//  Rcpp::NumericVector out_xx_j(npairs);
//  Rcpp::NumericVector out_xij(npairs);
//  Rcpp::NumericVector out_r(npairs);
//
//  auto decode_value = [&](unsigned int code, double p) -> double {
//   if (scale) {
//    const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//    if (denom <= 0.0 || !std::isfinite(denom)) {
//     return 0.0;
//    }
//
//    switch (code) {
//    case 0u: return (2.0 - 2.0 * p) / denom;
//    case 1u: return 0.0;
//    case 2u: return (1.0 - 2.0 * p) / denom;
//    default: return (0.0 - 2.0 * p) / denom;
//    }
//   }
//
//   switch (code) {
//   case 0u: return 2.0;
//   case 1u: return 2.0 * p;
//   case 2u: return 1.0;
//   default: return 0.0;
//   }
//  };
//
//  for (int k = 0; k < npairs; ++k) {
//   const int i = pairs(k, 0) - 1;
//   const int j = pairs(k, 1) - 1;
//
//   if (i < 0 || i >= m || j < 0 || j >= m) {
//    throw std::runtime_error("bed_pairwise_xtx_check: pairs contains marker index out of range.");
//   }
//
//   const uint8_t* row_i = G.row(i);
//   const uint8_t* row_j = G.row(j);
//
//   const double p_i = af_cpp[static_cast<std::size_t>(i)];
//   const double p_j = af_cpp[static_cast<std::size_t>(j)];
//
//   double xx_i = 0.0;
//   double xx_j = 0.0;
//   double xij = 0.0;
//
//   for (std::size_t kb = 0; kb < G.nbytes; ++kb) {
//    const uint8_t byte_i = row_i[kb];
//    const uint8_t byte_j = row_j[kb];
//    const int jbase = static_cast<int>(kb << 2);
//
//    for (int pos = 0; pos < 4; ++pos) {
//     const int sample = jbase + pos;
//
//     if (sample >= n_used) {
//      break;
//     }
//
//     const unsigned int code_i = (byte_i >> (2 * pos)) & 3u;
//     const unsigned int code_j = (byte_j >> (2 * pos)) & 3u;
//
//     const double xi = decode_value(code_i, p_i);
//     const double xj = decode_value(code_j, p_j);
//
//     xx_i += xi * xi;
//     xx_j += xj * xj;
//     xij += xi * xj;
//    }
//   }
//
//   const double denom = std::sqrt(xx_i * xx_j);
//   const double rij = denom > 0.0 ? xij / denom : NA_REAL;
//
//   out_i[k] = i + 1;
//   out_j[k] = j + 1;
//   out_xx_i[k] = xx_i;
//   out_xx_j[k] = xx_j;
//   out_xij[k] = xij;
//   out_r[k] = rij;
//  }
//
//  return Rcpp::DataFrame::create(
//   Rcpp::Named("i") = out_i,
//   Rcpp::Named("j") = out_j,
//   Rcpp::Named("xx_i_exact") = out_xx_i,
//   Rcpp::Named("xx_j_exact") = out_xx_j,
//   Rcpp::Named("xij_exact") = out_xij,
//   Rcpp::Named("r_exact") = out_r,
//   Rcpp::Named("n_used") = n_used,
//   Rcpp::Named("af_computed") = af_computed
//  );
// }
