#include "scanner_common.h"

double npwbs_baum_component(const std::vector<double>& ranks,
                            int other_size, int total_n) {
  int size = ranks.size();
  double total = 0.0;
  double scale = (total_n + 1.0) / (size + 1.0);
  double variance_scale = other_size * (total_n + 1.0) / (size + 2.0);
  for (int index = 0; index < size; ++index) {
    double i = index + 1.0;
    double p = i / (size + 1.0);
    double difference = ranks[index] - scale * i;
    total += difference * difference / (p * (1.0 - p) * variance_scale);
  }
  return total / size;
}

// [[Rcpp::export]]
Rcpp::List cpp_scan_baumgartner(Rcpp::NumericVector x, Rcpp::IntegerVector starts,
                                Rcpp::IntegerVector ends, int d = 2) {
  npwbs_validate(x, starts, ends, d);
  NpwbsBest overall;
  for (int interval = 0; interval < starts.size(); ++interval) {
    int n = ends[interval] - starts[interval] + 1;
    if (n < 2 * d) continue;
    std::vector<double> ranks = npwbs_ranks(x, starts[interval] - 1, ends[interval] - 1);
    std::vector<double> left, right = ranks;
    left.reserve(n);
    std::sort(right.begin(), right.end());
    double best = -std::numeric_limits<double>::infinity();
    int best_split = NA_INTEGER;
    for (int split = 1; split <= n - d; ++split) {
      double moved = ranks[split - 1];
      left.insert(std::lower_bound(left.begin(), left.end(), moved), moved);
      right.erase(std::lower_bound(right.begin(), right.end(), moved));
      if (split < d) continue;
      double value = 0.5 * (npwbs_baum_component(left, n - split, n) +
                            npwbs_baum_component(right, split, n));
      if (value > best) {
        best = value;
        best_split = starts[interval] + split - 1;
      }
    }
    npwbs_update(overall, best, best_split);
  }
  return npwbs_result(overall);
}
