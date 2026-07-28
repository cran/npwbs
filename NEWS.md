## npwbs 0.4.0

- Added Mann--Whitney, Mood, Cramér--von Mises, and modified Baumgartner detection alongside the existing Lepage method, using compact batch C++ scanners.
- Added whole-input random strict ordering for ties, controlled by `breakTies`.
- Lepage `alpha = 0.05` now uses the final registered book-chapter calibration and may produce small changepoint differences from version 0.3.0.
- Lepage `alpha = 0.01` retains the established package calibration because no corresponding book-chapter calibration was generated.
- Both bundled Lepage threshold vectors are now stored only through `n = 10000`.

## npwbs 0.3.0

- Updated built-in Lepage thresholds to match the revised manuscript calibration for the default settings `M = 10000` and `d = 2`.
- Clarified that the packaged thresholds support `alpha = 0.05` and `alpha = 0.01`.
- For segment lengths above `n = 10000`, the package now uses the `n = 10000` threshold with a warning instead of silently extrapolating.
- Clarified the returned changepoint convention: a returned value `k` denotes a split between observations `k` and `k + 1`.
