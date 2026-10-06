minimumSegmentLength <- function(test) {
  if (!test %in% c("lepage", "mw", "mood", "cvm", "baumgartner", "zhang", "anderson_darling")) {
    stop(sprintf("Error: unsupported method '%s'", test))
  }
  10L
}

thresholdCalibrationMaxN <- function() {
  10000L
}

prepareThresholds <- function(entry, yLength) {
  thresholds <- entry$thresholds
  if (entry$method == "zhang") {
    if (length(thresholds) < entry$max_n) {
      stop("Error: bundled Zhang threshold vector is shorter than 3000")
    }
    if (yLength > entry$max_n) {
      stop("The Zhang method currently supports sequences of length at most 3000.")
    }
    return(thresholds[seq_len(yLength)])
  }
  if (entry$method == "anderson_darling") {
    if (length(thresholds) < entry$max_n) {
      stop("Error: bundled Anderson-Darling threshold vector is shorter than 3000")
    }
    if (yLength > entry$max_n) {
      stop("The Anderson-Darling method currently supports sequences of length at most 3000.")
    }
    return(thresholds[seq_len(yLength)])
  }
  calibratedMaxN <- thresholdCalibrationMaxN()
  if (length(thresholds) < calibratedMaxN) {
    stop("Error: bundled threshold vector is shorter than the supported calibration range")
  }
  if (yLength > calibratedMaxN) {
    warning(
      "Thresholds are provided through n = 10000; the n = 10000 threshold will be used for longer segments."
    )
  }
  if (length(thresholds) > calibratedMaxN) {
    thresholds[(calibratedMaxN + 1L):length(thresholds)] <- thresholds[calibratedMaxN]
  }
  if (length(thresholds) < yLength) {
    thresholds <- c(thresholds, rep(thresholds[calibratedMaxN], yLength - length(thresholds)))
  }
  thresholds
}

sampleWBSInterval <- function(n, minlength) {
  if (n < minlength) stop("Error: n must be at least minlength")
  maxStart <- n - minlength + 1L
  start <- sample.int(maxStart, 1L)
  minEnd <- start + minlength - 1L
  end <- sample.int(n - minEnd + 1L, 1L) + minEnd - 1L
  c(start = start, end = end)
}

sampleWBSIntervals <- function(n, minlength, sims) {
  if (n < minlength) stop("Error: n must be at least minlength")
  if (length(sims) != 1L || !is.finite(sims) || sims < 1 || sims != floor(sims)) {
    stop("Error: sims must be a positive integer")
  }
  sims <- as.integer(sims)

  starts <- integer(sims)
  ends <- integer(sims)
  for (i in seq_len(sims)) {
    interval <- sampleWBSInterval(n, minlength)
    starts[i] <- interval[["start"]]
    ends[i] <- interval[["end"]]
  }
  list(starts = starts, ends = ends)
}

scanWBSBatch <- function(y, starts, ends, test, d = 2, combination = "sum") {
  switch(test,
    lepage = cpp_scan_lepage(y, starts, ends, d, combination),
    mw = cpp_scan_mw(y, starts, ends, d),
    mood = cpp_scan_mood(y, starts, ends, d),
    cvm = cpp_scan_cvm(y, starts, ends, d),
    baumgartner = cpp_scan_baumgartner(y, starts, ends, d),
    anderson_darling = cpp_scan_anderson_darling(y, starts, ends, d),
    zhang = {
      moments <- .zhang_moment_artifacts(max(ends - starts + 1L))
      cpp_scan_zhang(y, starts, ends, d, moments$bundled, moments$extended)
    },
    stop(sprintf("Error: unsupported method '%s'", test))
  )
}

WBS <- function(y, sims = 10000, test = "lepage", d = 2,
                combination = "sum") {
  minlength <- minimumSegmentLength(test)
  intervals <- sampleWBSIntervals(length(y), minlength, sims)
  scanWBSBatch(y, intervals$starts, intervals$ends, test, d, combination)
}

.thresholdKey <- function(method, combination, M, alpha, d) {
  alpha_key <- if (identical(alpha, 0.05)) {
    "0.05"
  } else if (identical(alpha, 0.01)) {
    "0.01"
  } else {
    stop("Error: alpha must be exactly 0.05 or 0.01")
  }
  paste(method, if (method == "lepage") combination else "scalar",
        as.integer(M), alpha_key,
        as.integer(d), sep = "|")
}

.thresholdRegistry <- function() {
  metadata <- function(method, combination, M, alpha, d, thresholds, max_n,
                       provenance, version) {
    list(
      method = method, combination = combination, M = as.integer(M),
      alpha = alpha, d = as.integer(d), minimum_interval_length = 10L,
      sampler = "legacy_asymmetric_with_replacement", forced_full_interval = FALSE,
      fresh_intervals_at_recursion = TRUE, fresh_intervals_at_pruning = TRUE,
      comparison = ">", min_n = 10L, max_n = as.integer(max_n),
      interpolation = "precomputed_integer_table",
      calibration_provenance = provenance, calibration_version = version,
      thresholds = thresholds
    )
  }
  entries <- list(
    metadata("lepage", "sum", 1000L, 0.05, 2L, LPthresholds05M1000,
             10000L, "verified WBS-Lepage revision calibration", "revision_20260913"),
    metadata("lepage", "sum", 1000L, 0.01, 2L, LPthresholds01M1000,
             10000L, "derived from verified M=1000 raw maxima; order statistic 9901", "revision_20260913_alpha01"),
    metadata("lepage", "max", 1000L, 0.05, 2L, LPMaxThresholds05M1000,
             10000L, "verified final WBS-max M=1000 calibration", "revision_20260913"),
    metadata("lepage", "max", 10000L, 0.05, 2L, LPMaxThresholds05M10000,
             10000L, "verified WBS-Lepage revision calibration", "revision_20260913"),
    metadata("lepage", "sum", 10000L, 0.05, 2L, LPthresholds05,
             10000L, "exact frozen npwbs 0.5.0 table", "npwbs_0.5.0"),
    metadata("lepage", "sum", 10000L, 0.01, 2L, LPthresholds01,
             10000L, "exact frozen npwbs 0.5.0 table", "npwbs_0.5.0"),
    metadata("mw", "scalar", 1000L, 0.05, 2L, MWthresholds05M1000,
             10000L, "verified M=1000 calibration", "revision_20260917"),
    metadata("mw", "scalar", 10000L, 0.05, 2L, MWthresholds05,
             10000L, "exact frozen npwbs 0.5.0 table", "npwbs_0.5.0"),
    metadata("mood", "scalar", 1000L, 0.05, 2L, Moodthresholds05M1000,
             10000L, "verified M=1000 calibration", "revision_20260917"),
    metadata("mood", "scalar", 10000L, 0.05, 2L, Moodthresholds05,
             10000L, "exact frozen npwbs 0.5.0 table", "npwbs_0.5.0"),
    metadata("cvm", "scalar", 1000L, 0.05, 2L, CVMthresholds05M1000,
             10000L, "verified M=1000 calibration", "revision_20260917"),
    metadata("cvm", "scalar", 10000L, 0.05, 2L, CVMthresholds05,
             10000L, "exact frozen npwbs 0.5.0 table", "npwbs_0.5.0"),
    metadata("baumgartner", "scalar", 1000L, 0.05, 2L, BaumgartnerThresholds05M1000,
             10000L, "verified M=1000 calibration", "revision_20260917"),
    metadata("baumgartner", "scalar", 10000L, 0.05, 2L, BaumgartnerThresholds05,
             10000L, "exact frozen npwbs 0.5.0 table", "npwbs_0.5.0"),
    metadata("zhang", "scalar", 1000L, 0.05, 4L, ZhangThresholds05M1000,
             3000L, "verified M=1000 calibration", "revision_20260917"),
    metadata("zhang", "scalar", 10000L, 0.05, 4L, ZhangThresholds05,
             3000L, "exact frozen npwbs 0.5.0 table", "npwbs_0.5.0"),
    metadata("anderson_darling", "scalar", 1000L, 0.05, 2L,
             AndersonDarlingThresholds05M1000, 3000L,
             "verified returned AD M=1000 calibration", "revision_20260923")
  )
  keys <- vapply(entries, function(x) {
    .thresholdKey(x$method, x$combination, x$M, x$alpha, x$d)
  }, character(1L))
  stats::setNames(entries, keys)
}

selectThresholdEntry <- function(method, combination, M, alpha, d) {
  key <- .thresholdKey(method, combination, M, alpha, d)
  entry <- .thresholdRegistry()[[key]]
  if (is.null(entry)) {
    stop(sprintf(
      paste0("Error: no built-in threshold calibration for method='%s', ",
             "combination='%s', M=%s, alpha=%s, d=%s"),
      method, combination, format(M), format(alpha), format(d)
    ))
  }
  entry
}

detectChanges <- function(y, alpha = 0.05, prune = TRUE, M = 1000, d = 2,
                          displayOutput = FALSE, method = "lepage", breakTies = TRUE,
                          combination = "sum") {
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) == 0L || any(!is.finite(y))) {
    stop("Error: y must be a non-empty finite numeric vector")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("lepage", "mw", "mood", "cvm", "baumgartner", "zhang", "anderson_darling")) {
    stop("Error: method must be exactly one of 'lepage', 'mw', 'mood', 'cvm', 'baumgartner', 'zhang', or 'anderson_darling'")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || !is.finite(alpha)) {
    stop("Error: alpha must be a single supported numeric value")
  }
  if (!alpha %in% c(0.05, 0.01)) {
    stop("Error: alpha must be exactly 0.05 or 0.01")
  }
  if (!is.character(combination) || length(combination) != 1L || is.na(combination) ||
      !combination %in% c("sum", "max")) {
    stop("Error: combination must be exactly 'sum' or 'max'")
  }
  if (method != "lepage" && combination != "sum") {
    stop("Error: combination='max' is supported only for method='lepage'")
  }
  if (!is.numeric(M) || length(M) != 1L || is.na(M) || !is.finite(M) ||
      M != floor(M) || !M %in% c(1000, 10000)) {
    stop("Error: M must be exactly 1000 or 10000 with a matching built-in calibration")
  }
  if (method == "zhang") {
    if (!missing(d) &&
        (!is.numeric(d) || length(d) != 1L || is.na(d) || !is.finite(d) || d != 4)) {
      stop("Error: method='zhang' supports only d=4")
    }
    d <- 4L
  } else if (!is.numeric(d) || length(d) != 1L || is.na(d) ||
             !is.finite(d) || d != 2) {
    stop("Error: bundled thresholds support only d=2")
  }
  if (!is.logical(prune) || length(prune) != 1L || is.na(prune)) {
    stop("Error: prune must be TRUE or FALSE")
  }
  if (!is.logical(breakTies) || length(breakTies) != 1L || is.na(breakTies)) {
    stop("Error: breakTies must be TRUE or FALSE")
  }
  if (!is.logical(displayOutput) || length(displayOutput) != 1L || is.na(displayOutput)) {
    stop("Error: displayOutput must be TRUE or FALSE")
  }

  entry <- selectThresholdEntry(method, combination, M, alpha, d)
  if (method == "zhang") {
    .zhang_moment_artifacts(length(y))
  }
  thresholds <- prepareThresholds(entry, length(y))
  if (anyDuplicated(y)) {
    if (!breakTies) stop("Error: ties are present and breakTies=FALSE")
    warning("Ties detected; applying one random strict ordering to the complete input.")
    y <- rank(y, ties.method = "random")
  }

  cps <- detectChanges_aux(y, 1, length(y), method, thresholds, M, d,
                           displayOutput, combination = combination)
  if (prune) cps <- prunecps(y, cps, method, thresholds, M = M, d = d,
                             combination = combination)
  cps
}

detectChanges_aux <- function(y, start, end, test, thresholds, M = 10000, d = 2,
                              displayOutput = FALSE,
                              combination = "sum") {
  n <- end - start + 1L
  if (n < minimumSegmentLength(test)) return(numeric())
  val <- WBS(y[start:end], M, test, d, combination)
  if (!is.finite(val$value) || val$value <= thresholds[n]) return(numeric())
  k <- start + val$split - 1L
  if (displayOutput) {
    print(sprintf("found change point at %s on [%s,%s] with statistic %s",
                  k, start, end, val$value))
  }
  c(
    detectChanges_aux(y, start, k, test, thresholds, M, d,
                      displayOutput = displayOutput, combination = combination),
    k,
    detectChanges_aux(y, k + 1L, end, test, thresholds, M, d,
                      displayOutput = displayOutput, combination = combination)
  )
}

prunecps <- function(y, cps, test, thresholds, M = 10000, d = 2,
                     combination = "sum") {
  minlength <- minimumSegmentLength(test)
  cps <- c(0, cps, length(y))
  repeat {
    deleted <- FALSE
    if (length(cps) < 3L) break
    i <- 2L
    repeat {
      if (i >= length(cps)) break
      n <- cps[i + 1L] - cps[i - 1L]
      if (n < minlength) {
        cps <- cps[-i]
        deleted <- TRUE
        next
      }
      val <- WBS(y[(cps[i - 1L] + 1L):cps[i + 1L]], M, test, d,
                 combination)
      if (!is.finite(val$value) || val$value <= thresholds[n]) {
        cps <- cps[-i]
        deleted <- TRUE
      } else {
        i <- i + 1L
      }
    }
    if (!deleted) break
  }
  cps[-c(1L, length(cps))]
}
