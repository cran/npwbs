#include "scanner_common.h"

class NpwbsLongDoubleFenwick {
  int n;
  std::vector<long double> tree;
public:
  explicit NpwbsLongDoubleFenwick(int size) : n(size), tree(size + 1, 0.0L) {}
  void add(int index, long double value) {
    for (int i = index; i <= n; i += i & -i) tree[i] += value;
  }
  long double sum(int index) const {
    long double out = 0.0L;
    for (int i = index; i > 0; i -= i & -i) out += tree[i];
    return out;
  }
};

static std::vector<long double> ad_harmonic(1, 0.0L);
static std::vector<long double> ad_g(3, 0.0L);

static void ensure_ad_cache(int n) {
  while (static_cast<int>(ad_harmonic.size()) <= n) {
    int k = static_cast<int>(ad_harmonic.size());
    ad_harmonic.push_back(ad_harmonic.back() + 1.0L / static_cast<long double>(k));
  }
  while (static_cast<int>(ad_g.size()) <= n) {
    int current_n = static_cast<int>(ad_g.size()) - 1;
    long double nd = static_cast<long double>(current_n);
    ad_g.push_back(ad_g[current_n] + 1.0L / (nd * nd) +
      2.0L * (ad_harmonic[current_n - 1] - 1.0L) / (nd * (nd + 1.0L)));
  }
}

static long double ad_null_variance(int n, int split) {
  if (n < 4 || split < 1 || split >= n) return NA_REAL;
  ensure_ad_cache(n);
  long double nd = static_cast<long double>(n);
  long double m = static_cast<long double>(split);
  long double q = nd - m;
  long double H = 1.0L / m + 1.0L / q;
  long double h = ad_harmonic[n - 1];
  long double g = ad_g[n];
  long double a = 4.0L * g - 6.0L + (10.0L - 6.0L * g) * H;
  long double b = 12.0L * g + 8.0L * h - 22.0L +
    (2.0L * g - 14.0L * h - 4.0L) * H;
  long double c = 36.0L * h + 4.0L + (2.0L * h - 6.0L) * H;
  return (a * nd * nd * nd + b * nd * nd + c * nd + 24.0L) /
    ((nd - 1.0L) * (nd - 2.0L) * (nd - 3.0L));
}

// [[Rcpp::export]]
Rcpp::List cpp_scan_anderson_darling(Rcpp::NumericVector x,
                                     Rcpp::IntegerVector starts,
                                     Rcpp::IntegerVector ends, int d = 2) {
  npwbs_validate(x, starts, ends, d);
  NpwbsBest overall;
  for (int interval = 0; interval < starts.size(); ++interval) {
    int n = ends[interval] - starts[interval] + 1;
    if (n < 2 * d || n < 4) continue;
    std::vector<double> ranks = npwbs_ranks(x, starts[interval] - 1, ends[interval] - 1);
    ensure_ad_cache(n);
    long double nd = static_cast<long double>(n);
    long double harmonic_n_minus_1 = ad_harmonic[n - 1];
    std::vector<long double> K(n + 1, 0.0L), L1(n + 1, 0.0L);
    for (int r = 1; r <= n; ++r) {
      K[r] = (harmonic_n_minus_1 - ad_harmonic[r - 1] + ad_harmonic[n - r]) / nd;
      L1[r] = ad_harmonic[n - r];
    }
    long double C0 = nd * harmonic_n_minus_1 - (nd - 1.0L);
    NpwbsLongDoubleFenwick counts(n), weights(n);
    long double total_weight = 0.0L, T1 = 0.0L, T2 = 0.0L;
    double best = -std::numeric_limits<double>::infinity();
    int best_split = NA_INTEGER;
    for (int split = 1; split <= n - 1; ++split) {
      int rank = static_cast<int>(std::round(ranks[split - 1]));
      if (rank < 1 || rank > n ||
          std::fabs(ranks[split - 1] - static_cast<double>(rank)) > 1e-8) {
        Rcpp::stop("Anderson-Darling requires integer ranks 1:N; ties are not supported");
      }
      long double count_le = counts.sum(rank);
      long double weight_le = weights.sum(rank);
      T1 += 2.0L * (count_le * K[rank] + total_weight - weight_le) + K[rank];
      T2 += L1[rank];
      counts.add(rank, 1.0L);
      weights.add(rank, K[rank]);
      total_weight += K[rank];
      if (split < d || split > n - d) continue;
      long double m = static_cast<long double>(split), q = nd - m;
      long double raw = (nd * nd * T1 - 2.0L * nd * m * T2 + m * m * C0) / (m * q);
      long double variance = ad_null_variance(n, split);
      double value = -std::numeric_limits<double>::infinity();
      if (std::isfinite(static_cast<double>(raw)) &&
          std::isfinite(static_cast<double>(variance)) && variance > 0.0L) {
        value = static_cast<double>((raw - 1.0L) / std::sqrt(variance));
      }
      if (R_finite(value) && value > best) {
        best = value;
        best_split = starts[interval] + split - 1;
      }
    }
    npwbs_update(overall, best, best_split);
  }
  return npwbs_result(overall);
}
