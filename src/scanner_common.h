#ifndef NPWBS_SCANNER_COMMON_H
#define NPWBS_SCANNER_COMMON_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

struct NpwbsBest {
  double value;
  int split;
  NpwbsBest() : value(-std::numeric_limits<double>::infinity()),
                split(NA_INTEGER) {}
};

inline void npwbs_validate(const Rcpp::NumericVector& x,
                           const Rcpp::IntegerVector& starts,
                           const Rcpp::IntegerVector& ends, int d) {
  if (starts.size() != ends.size()) Rcpp::stop("starts and ends must have equal length");
  if (d < 1) Rcpp::stop("d must be positive");
  for (int i = 0; i < x.size(); ++i) {
    if (!R_finite(x[i])) Rcpp::stop("y must contain only finite values");
  }
  for (int i = 0; i < starts.size(); ++i) {
    if (starts[i] == NA_INTEGER || ends[i] == NA_INTEGER ||
        starts[i] < 1 || ends[i] < starts[i] || ends[i] > x.size()) {
      Rcpp::stop("invalid scan interval");
    }
  }
}

inline std::vector<double> npwbs_ranks(const Rcpp::NumericVector& x,
                                       int start0, int end0) {
  int n = end0 - start0 + 1;
  std::vector<std::pair<double, int> > ordered;
  ordered.reserve(n);
  for (int i = 0; i < n; ++i) ordered.push_back(std::make_pair(x[start0 + i], i));
  std::sort(ordered.begin(), ordered.end(),
            [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
              if (a.first != b.first) return a.first < b.first;
              return a.second < b.second;
            });
  std::vector<double> ranks(n);
  int pos = 0;
  while (pos < n) {
    int next = pos + 1;
    while (next < n && ordered[next].first == ordered[pos].first) ++next;
    double rank = (static_cast<double>(pos + 1) + static_cast<double>(next)) / 2.0;
    for (int j = pos; j < next; ++j) ranks[ordered[j].second] = rank;
    pos = next;
  }
  return ranks;
}

inline void npwbs_update(NpwbsBest& best, double value, int split) {
  if (R_finite(value) && value > best.value) {
    best.value = value;
    best.split = split;
  }
}

inline Rcpp::List npwbs_result(const NpwbsBest& best) {
  return Rcpp::List::create(Rcpp::Named("value") = best.value,
                            Rcpp::Named("split") = best.split);
}

#endif
