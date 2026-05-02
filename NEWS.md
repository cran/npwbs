## npwbs 0.3.0

- Updated built-in Lepage thresholds to match the revised manuscript calibration for the default settings `M = 10000` and `d = 2`.
- Clarified that the packaged thresholds support `alpha = 0.05` and `alpha = 0.01`.
- For segment lengths above `n = 10000`, the package now uses the `n = 10000` threshold with a warning instead of silently extrapolating.
- Clarified the returned changepoint convention: a returned value `k` denotes a split between observations `k` and `k + 1`.
