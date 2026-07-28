#include "scanner_common.h"

// [[Rcpp::export]]
Rcpp::List cpp_scan_mw(Rcpp::NumericVector x, Rcpp::IntegerVector starts,
                       Rcpp::IntegerVector ends, int d = 2) {
  npwbs_validate(x, starts, ends, d);
  NpwbsBest overall;
  for (int interval = 0; interval < starts.size(); ++interval) {
    int n = ends[interval] - starts[interval] + 1;
    if (n < 2 * d) continue;
    std::vector<double> ranks = npwbs_ranks(x, starts[interval] - 1, ends[interval] - 1);
    double sum = 0.0;
    double best = -std::numeric_limits<double>::infinity();
    int best_split = NA_INTEGER;
    for (int split = 1; split <= n - d; ++split) {
      sum += ranks[split - 1];
      if (split < d) continue;
      double n1 = split, n2 = n - split;
      double u = sum - n1 * (n1 + 1.0) / 2.0;
      double z = u - n1 * n2 / 2.0;
      double value = z * z / (n1 * n2 * (n + 1.0) / 12.0);
      if (value > best) {
        best = value;
        best_split = starts[interval] + split - 1;
      }
    }
    npwbs_update(overall, best, best_split);
  }
  return npwbs_result(overall);
}
