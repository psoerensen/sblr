#ifndef SBLR_BLR_CSR_BAYESC_OPERATOR_ADAPTER_H
#define SBLR_BLR_CSR_BAYESC_OPERATOR_ADAPTER_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include <Rcpp.h>

#include "st_csr_common.h"

struct LDLDFriends {
 std::vector<uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

template <class OpT>
struct BayescOperatorContext {
 OpT op;
 LDLDFriends ld_swap_friends;
 Rcpp::List diagnostics;

 BayescOperatorContext(
   const OpT& op_, const LDLDFriends& friends_, const Rcpp::List& diagnostics_
 ) : op(op_), ld_swap_friends(friends_), diagnostics(diagnostics_) {}

 BayescOperatorContext(
   OpT&& op_, LDLDFriends&& friends_, const Rcpp::List& diagnostics_
 ) : op(std::move(op_)), ld_swap_friends(std::move(friends_)),
     diagnostics(diagnostics_) {}
};

inline LDLDFriends build_ld_swap_friends_st_csr(
  int m, const STLDCSR& ld, const std::vector<double>& xx,
  double min_r2, int max_friends
) {
 std::vector<std::vector<std::pair<int, double>>> rows(
  static_cast<std::size_t>(m));
 for (int i = 0; i < m; ++i) {
  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end = ld.ptr[static_cast<std::size_t>(i + 1)];
  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   if (j <= i) continue;
   const double denom = xx[static_cast<std::size_t>(i)] *
    xx[static_cast<std::size_t>(j)];
   if (!std::isfinite(denom) || denom <= 0.0) continue;
   const double xij = static_cast<double>(ld.xij[static_cast<std::size_t>(p)]);
   const double r2 = (xij * xij) / denom;
   if (!std::isfinite(r2) || r2 < min_r2) continue;
   rows[static_cast<std::size_t>(i)].push_back(std::make_pair(j, r2));
   rows[static_cast<std::size_t>(j)].push_back(std::make_pair(i, r2));
  }
 }
 LDLDFriends friends;
 friends.ptr.resize(static_cast<std::size_t>(m) + 1);
 friends.ptr[0] = 0;
 for (int i = 0; i < m; ++i) {
  std::vector<std::pair<int, double>>& row = rows[static_cast<std::size_t>(i)];
  std::sort(row.begin(), row.end(),
   [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
    if (a.second == b.second) return a.first < b.first;
    return a.second > b.second;
   });
  std::vector<std::pair<int, double>> unique_row;
  unique_row.reserve(row.size());
  for (const auto& value : row) {
   if (!unique_row.empty() && unique_row.back().first == value.first) continue;
   unique_row.push_back(value);
  }
  if (static_cast<int>(unique_row.size()) > max_friends) {
   unique_row.resize(static_cast<std::size_t>(max_friends));
  }
  row.swap(unique_row);
  friends.ptr[static_cast<std::size_t>(i + 1)] =
   friends.ptr[static_cast<std::size_t>(i)] + row.size();
 }
 const uint64_t nfriend = friends.ptr[static_cast<std::size_t>(m)];
 friends.idx.resize(static_cast<std::size_t>(nfriend));
 friends.r2.resize(static_cast<std::size_t>(nfriend));
 for (int i = 0; i < m; ++i) {
  const uint64_t offset = friends.ptr[static_cast<std::size_t>(i)];
  const auto& row = rows[static_cast<std::size_t>(i)];
  for (std::size_t k = 0; k < row.size(); ++k) {
   friends.idx[static_cast<std::size_t>(offset + k)] = row[k].first;
   friends.r2[static_cast<std::size_t>(offset + k)] = row[k].second;
  }
 }
 return friends;
}

#endif
