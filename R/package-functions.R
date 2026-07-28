minimumSegmentLength <- function(test) {
  if (!test %in% c("lepage", "mw", "mood", "cvm", "baumgartner")) {
    stop(sprintf("Error: unsupported method '%s'", test))
  }
  10L
}

thresholdCalibrationMaxN <- function() {
  10000L
}

prepareThresholds <- function(thresholds, yLength, method = "lepage") {
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

scanWBSBatch <- function(y, starts, ends, test, d = 2) {
  switch(test,
    lepage = cpp_scan_lepage(y, starts, ends, d),
    mw = cpp_scan_mw(y, starts, ends, d),
    mood = cpp_scan_mood(y, starts, ends, d),
    cvm = cpp_scan_cvm(y, starts, ends, d),
    baumgartner = cpp_scan_baumgartner(y, starts, ends, d),
    stop(sprintf("Error: unsupported method '%s'", test))
  )
}

WBS <- function(y, sims = 10000, test = "lepage", d = 2) {
  starts <- integer(sims)
  ends <- integer(sims)
  minlength <- minimumSegmentLength(test)
  for (i in seq_len(sims)) {
    interval <- sampleWBSInterval(length(y), minlength)
    starts[i] <- interval[["start"]]
    ends[i] <- interval[["end"]]
  }
  scanWBSBatch(y, starts, ends, test, d)
}

selectThresholds <- function(method, alpha) {
  if (method == "lepage") {
    if (!alpha %in% c(0.05, 0.01)) {
      stop("Error: method='lepage' supports only alpha=0.05 or alpha=0.01")
    }
    return(if (alpha == 0.05) LPthresholds05 else LPthresholds01)
  }
  if (alpha != 0.05) {
    stop(sprintf("Error: method='%s' supports only alpha=0.05", method))
  }
  switch(method,
    mw = MWthresholds05,
    mood = Moodthresholds05,
    cvm = CVMthresholds05,
    baumgartner = BaumgartnerThresholds05
  )
}

detectChanges <- function(y, alpha = 0.05, prune = TRUE, M = 10000, d = 2,
                          displayOutput = FALSE, method = "lepage", breakTies = TRUE) {
  if (!is.numeric(y) || length(y) == 0L || any(!is.finite(y))) {
    stop("Error: y must be a non-empty finite numeric vector")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("lepage", "mw", "mood", "cvm", "baumgartner")) {
    stop("Error: method must be exactly one of 'lepage', 'mw', 'mood', 'cvm', or 'baumgartner'")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || !is.finite(alpha)) {
    stop("Error: alpha must be a single supported numeric value")
  }
  if (!is.numeric(M) || length(M) != 1L || is.na(M) || !is.finite(M) || M != 10000) {
    stop("Error: bundled thresholds support only M=10000")
  }
  if (!is.numeric(d) || length(d) != 1L || is.na(d) || !is.finite(d) || d != 2) {
    stop("Error: bundled thresholds support only d=2")
  }
  if (!is.logical(prune) || length(prune) != 1L || is.na(prune)) {
    stop("Error: prune must be TRUE or FALSE")
  }
  if (!is.logical(breakTies) || length(breakTies) != 1L || is.na(breakTies)) {
    stop("Error: breakTies must be TRUE or FALSE")
  }

  thresholds <- prepareThresholds(selectThresholds(method, alpha), length(y), method)
  if (anyDuplicated(y)) {
    if (!breakTies) stop("Error: ties are present and breakTies=FALSE")
    warning("Ties detected; applying one random strict ordering to the complete input.")
    y <- rank(y, ties.method = "random")
  }

  cps <- detectChanges_aux(y, 1, length(y), method, thresholds, M, d, displayOutput)
  if (prune) cps <- prunecps(y, cps, method, thresholds, M = M, d = d)
  cps
}

detectChanges_aux <- function(y, start, end, test, thresholds, M = 10000, d = 2,
                              displayOutput = FALSE) {
  n <- end - start + 1L
  if (n < minimumSegmentLength(test)) return(numeric())
  val <- WBS(y[start:end], M, test, d)
  if (!is.finite(val$value) || val$value <= thresholds[n]) return(numeric())
  k <- start + val$split - 1L
  if (displayOutput) {
    print(sprintf("found change point at %s on [%s,%s] with statistic %s",
                  k, start, end, val$value))
  }
  c(
    detectChanges_aux(y, start, k, test, thresholds, M, d),
    k,
    detectChanges_aux(y, k + 1L, end, test, thresholds, M, d)
  )
}

prunecps <- function(y, cps, test, thresholds, M = 10000, d = 2) {
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
        next
      }
      val <- WBS(y[(cps[i - 1L] + 1L):cps[i + 1L]], M, test, d)
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
