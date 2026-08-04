#include "scanner_common.h"
#include <memory>

struct NpwbsZhangMoments {
  int min_n;
  int max_n;
  Rcpp::IntegerVector offsets;
  Rcpp::NumericVector mean_zplus;
  Rcpp::NumericVector sd_zplus;

  explicit NpwbsZhangMoments(const Rcpp::List& artifact) {
    Rcpp::List metadata = artifact["metadata"];
    Rcpp::List data = artifact["data"];
    min_n = Rcpp::as<int>(metadata["min_N"]);
    max_n = Rcpp::as<int>(metadata["max_N"]);
    offsets = data["offsets"];
    mean_zplus = data["mean_zplus"];
    sd_zplus = data["sd_zplus"];
  }

  void lookup(int n, int split, double& mean, double& sd) const {
    if (n < min_n || n > max_n) Rcpp::stop("Zhang moment table does not cover interval length");
    int position = n - min_n;
    int m_star = std::min(split, n - split);
    int index = offsets[position] + m_star - 4;
    mean = mean_zplus[index];
    sd = sd_zplus[index];
  }
};

class NpwbsZhangWeightCache {
private:
  std::vector<std::vector<double> > weights_;
  std::vector<unsigned char> ready_;

public:
  explicit NpwbsZhangWeightCache(int max_n)
    : weights_(max_n + 1), ready_(max_n + 1, 0) {}

  const std::vector<double>& get(int size) {
    if (!ready_[size]) {
      weights_[size].resize(size);
      for (int j = 1; j <= size; ++j) {
        weights_[size][j - 1] = std::log(
          static_cast<double>(size) /
            (static_cast<double>(j) - 0.5) - 1.0
        );
      }
      ready_[size] = 1;
    }
    return weights_[size];
  }
};

double npwbs_zhang_component(const std::vector<int>& sorted_ranks,
                             const std::vector<double>& a,
                             const std::vector<double>& b) {
  double total = 0.0;
  for (int i = 0; i < static_cast<int>(sorted_ranks.size()); ++i) {
    total += a[i] * b[sorted_ranks[i] - 1];
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::List cpp_scan_zhang(Rcpp::NumericVector x, Rcpp::IntegerVector starts,
                          Rcpp::IntegerVector ends, int d,
                          Rcpp::List bundled,
                          Rcpp::Nullable<Rcpp::List> extended = R_NilValue) {
  npwbs_validate(x, starts, ends, d);
  if (d != 4) Rcpp::stop("Zhang requires d=4");

  NpwbsZhangMoments bundled_moments(bundled);
  std::unique_ptr<NpwbsZhangMoments> extended_moments;
  if (extended.isNotNull()) {
    extended_moments.reset(new NpwbsZhangMoments(Rcpp::List(extended)));
  }

  int max_n = 0;
  for (int interval = 0; interval < starts.size(); ++interval) {
    max_n = std::max(max_n, ends[interval] - starts[interval] + 1);
  }
  NpwbsZhangWeightCache weights(max_n);

  NpwbsBest overall;
  for (int interval = 0; interval < starts.size(); ++interval) {
    int start0 = starts[interval] - 1;
    int n = ends[interval] - starts[interval] + 1;
    if (n < 2 * d) continue;

    const NpwbsZhangMoments* moments = &bundled_moments;
    if (n > bundled_moments.max_n) {
      if (!extended_moments) Rcpp::stop("Extended Zhang moments are not loaded");
      moments = extended_moments.get();
    }

    std::vector<double> ranks = npwbs_ranks(x, start0, ends[interval] - 1);
    std::vector<int> int_ranks(n);
    std::vector<bool> seen(n + 1, false);
    for (int i = 0; i < n; ++i) {
      int rank = static_cast<int>(std::round(ranks[i]));
      if (rank < 1 || rank > n ||
          std::fabs(ranks[i] - static_cast<double>(rank)) > 1e-8 ||
          seen[rank]) {
        Rcpp::stop("Zhang requires a strict rank ordering");
      }
      int_ranks[i] = rank;
      seen[rank] = true;
    }

    const std::vector<double>& b = weights.get(n);

    std::vector<int> left;
    std::vector<int> right;
    left.reserve(n);
    right.reserve(n);
    for (int i = 0; i < d; ++i) left.push_back(int_ranks[i]);
    for (int i = d; i < n; ++i) right.push_back(int_ranks[i]);
    std::sort(left.begin(), left.end());
    std::sort(right.begin(), right.end());

    double best = -std::numeric_limits<double>::infinity();
    int best_split = NA_INTEGER;
    for (int split = d; split <= n - d; ++split) {
      const std::vector<double>& a_left = weights.get(split);
      const std::vector<double>& a_right = weights.get(n - split);
      double zc = (
        npwbs_zhang_component(left, a_left, b) +
        npwbs_zhang_component(right, a_right, b)
      ) / static_cast<double>(n);
      double mean, sd;
      moments->lookup(n, split, mean, sd);
      double value = ((4.0 - zc) - mean) / sd;
      if (R_finite(value) && value > best) {
        best = value;
        best_split = starts[interval] + split - 1;
      }

      if (split < n - d) {
        int moved = int_ranks[split];
        left.insert(std::lower_bound(left.begin(), left.end(), moved), moved);
        std::vector<int>::iterator position =
          std::lower_bound(right.begin(), right.end(), moved);
        if (position == right.end() || *position != moved) {
          Rcpp::stop("Internal Zhang rank bookkeeping failure");
        }
        right.erase(position);
      }
    }
    npwbs_update(overall, best, best_split);
  }
  return npwbs_result(overall);
}
