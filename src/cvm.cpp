#include "scanner_common.h"

class NpwbsFenwick {
  int n;
  std::vector<double> tree;
public:
  explicit NpwbsFenwick(int size) : n(size), tree(size + 1, 0.0) {}
  void add(int index, double value) {
    for (int i = index; i <= n; i += i & -i) tree[i] += value;
  }
  double sum(int index) const {
    double out = 0.0;
    for (int i = index; i > 0; i -= i & -i) out += tree[i];
    return out;
  }
};

// [[Rcpp::export]]
Rcpp::List cpp_scan_cvm(Rcpp::NumericVector x, Rcpp::IntegerVector starts,
                        Rcpp::IntegerVector ends, int d = 2) {
  npwbs_validate(x, starts, ends, d);
  NpwbsBest overall;
  for (int interval = 0; interval < starts.size(); ++interval) {
    int n = ends[interval] - starts[interval] + 1;
    if (n < 2 * d) continue;
    std::vector<double> ranks = npwbs_ranks(x, starts[interval] - 1, ends[interval] - 1);
    NpwbsFenwick counts(n), weights(n);
    double saa = 0.0, sja = 0.0, total_weight = 0.0;
    double best = -std::numeric_limits<double>::infinity();
    int best_split = NA_INTEGER;
    for (int split = 1; split <= n - d; ++split) {
      int rank = static_cast<int>(std::round(ranks[split - 1]));
      double suffix = n - rank + 1.0;
      double suffix_a = counts.sum(rank) * suffix + total_weight - weights.sum(rank);
      saa += 2.0 * suffix_a + suffix;
      sja += (rank + static_cast<double>(n)) * suffix / 2.0;
      counts.add(rank, 1.0);
      weights.add(rank, suffix);
      total_weight += suffix;
      if (split < d) continue;
      double nd = n, m = split, n2 = nd - m;
      double sj2 = nd * (nd + 1.0) * (2.0 * nd + 1.0) / 6.0;
      double w = (nd * nd * saa - 2.0 * nd * m * sja + m * m * sj2) /
        (nd * nd * m * n2);
      double mean = (nd + 1.0) / (6.0 * nd);
      double variance = (nd + 1.0) *
        (4.0 * m * n2 * nd - 3.0 * m * m - 3.0 * n2 * n2 - 2.0 * m * n2) /
        (180.0 * m * n2 * nd * nd);
      double value = (w - mean) / std::sqrt(variance);
      if (R_finite(value) && value > best) {
        best = value;
        best_split = starts[interval] + split - 1;
      }
    }
    npwbs_update(overall, best, best_split);
  }
  return npwbs_result(overall);
}
