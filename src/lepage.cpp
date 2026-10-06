#include "scanner_common.h"

// [[Rcpp::export]]
Rcpp::List cpp_scan_lepage(Rcpp::NumericVector x, Rcpp::IntegerVector starts,
                           Rcpp::IntegerVector ends, int d = 2,
                           std::string combination = "sum") {
  npwbs_validate(x, starts, ends, d);
  if (combination != "sum" && combination != "max") {
    Rcpp::stop("combination must be exactly 'sum' or 'max'");
  }
  NpwbsBest overall;
  for (int interval = 0; interval < starts.size(); ++interval) {
    int start0 = starts[interval] - 1;
    int n = ends[interval] - starts[interval] + 1;
    if (n < 2 * d) continue;
    std::vector<double> ranks = npwbs_ranks(x, start0, ends[interval] - 1);
    std::vector<double> rank_sum(n), mood_sum(n);
    double nd = static_cast<double>(n);
    double center = (nd + 1.0) / 2.0;
    for (int i = 0; i < n; ++i) {
      double mood = (ranks[i] - center) * (ranks[i] - center);
      rank_sum[i] = ranks[i] + (i ? rank_sum[i - 1] : 0.0);
      mood_sum[i] = mood + (i ? mood_sum[i - 1] : 0.0);
    }
    double best = -std::numeric_limits<double>::infinity();
    int best_split = NA_INTEGER;
    for (int split = d; split <= n - d; ++split) {
      double n1 = split, n2 = n - split;
      double u = rank_sum[split - 1] - n1 * (n1 + 1.0) / 2.0;
      double mw = u - n1 * n2 / 2.0;
      mw = mw * mw / (n1 * n2 * (nd + 1.0) / 12.0);
      double mood = mood_sum[split - 1] - n1 * (nd * nd - 1.0) / 12.0;
      mood = mood * mood / (n1 * n2 * (nd + 1.0) * (nd * nd - 4.0) / 180.0);
      double value = combination == "sum" ? mw + mood : std::max(mw, mood);
      if (value > best) {
        best = value;
        best_split = starts[interval] + split - 1;
      }
    }
    npwbs_update(overall, best, best_split);
  }
  return npwbs_result(overall);
}
