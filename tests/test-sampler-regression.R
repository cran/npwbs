library(npwbs)

ns <- asNamespace("npwbs")
get_ns <- function(name) get(name, envir = ns)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

legacy_reference <- function(n, minlength, sims) {
  starts <- integer(sims)
  ends <- integer(sims)
  for (i in seq_len(sims)) {
    max_start <- n - minlength + 1L
    starts[i] <- sample.int(max_start, 1L)
    min_end <- starts[i] + minlength - 1L
    ends[i] <- sample.int(n - min_end + 1L, 1L) + min_end - 1L
  }
  list(starts = starts, ends = ends)
}

for (n in c(10L, 11L, 100L, 1000L)) {
  set.seed(5100L + n)
  expected <- legacy_reference(n, 10L, 1000L)
  set.seed(5100L + n)
  actual <- get_ns("sampleWBSIntervals")(n, 10L, 1000L)
  stopifnot(identical(actual, expected))
}

mann_whitney_reference <- function(y, d = 2L) {
  ranks <- rank(y)
  n <- length(ranks)
  values <- numeric(n)
  rank_sum <- if (d > 1L) sum(ranks[seq_len(d - 1L)]) else 0
  for (i in seq.int(d, n - d + 1L)) {
    n1 <- i
    n2 <- n - n1
    rank_sum <- rank_sum + ranks[i]
    u <- rank_sum - n1 * (n1 + 1) / 2
    mean_u <- n1 * n2 / 2
    sd_u <- sqrt(n1 * n2 * (n + 1) / 12)
    values[i] <- abs((u - mean_u) / sd_u)
  }
  if (d > 1L) values[seq_len(d - 1L)] <- -Inf
  values[(n - d + 1L):n] <- -Inf
  values
}

mood_reference <- function(y, d = 2L) {
  ranks <- (rank(y) - (length(y) + 1) / 2)^2
  n <- length(ranks)
  values <- numeric(n)
  rank_sum <- if (d > 1L) sum(ranks[seq_len(d - 1L)]) else 0
  for (i in seq.int(d, n - d + 1L)) {
    n1 <- i
    n2 <- n - n1
    rank_sum <- rank_sum + ranks[i]
    mean_mood <- n1 * (n^2 - 1) / 12
    sd_mood <- sqrt(n1 * n2 * (n + 1) * (n^2 - 4) / 180)
    values[i] <- abs((rank_sum - mean_mood) / sd_mood)
  }
  if (d > 1L) values[seq_len(d - 1L)] <- -Inf
  values[(n - d + 1L):n] <- -Inf
  values
}

lepage_reference_interval <- function(y, start, end, d = 2L) {
  local_y <- y[start:end]
  values <- mann_whitney_reference(local_y, d)^2 + mood_reference(local_y, d)^2
  if (d > 1L) values[seq_len(d - 1L)] <- -Inf
  values[(length(local_y) - d + 1L):length(local_y)] <- -Inf
  best <- which.max(values)
  list(value = values[best], split = as.integer(start + best - 1L))
}

lepage_reference_batch <- function(y, starts, ends, d = 2L) {
  results <- Map(
    function(start, end) lepage_reference_interval(y, start, end, d),
    starts, ends
  )
  values <- vapply(results, `[[`, numeric(1), "value")
  best <- which.max(values)
  list(value = values[best], split = results[[best]]$split)
}

expect_scanner_equal <- function(y, starts, ends, tolerance = 1e-10) {
  expected <- lepage_reference_batch(y, starts, ends)
  actual <- get_ns("cpp_scan_lepage")(y, as.integer(starts), as.integer(ends), 2L)
  scale <- max(1, abs(expected$value))
  stopifnot(
    abs(actual$value - expected$value) <= tolerance * scale,
    identical(as.integer(actual$split), as.integer(expected$split))
  )
  invisible(actual)
}

for (n in c(10:100, 125L, 250L, 500L, 750L, 1000L)) {
  set.seed(6000L + n)
  expect_scanner_equal(rnorm(n), 1L, n)
}

set.seed(7001)
y_collection <- rnorm(1000)
starts <- c(1L, 991L, 1L, 251L, 101L, 400L, 73L, 500L)
ends <- c(10L, 1000L, 1000L, 750L, 110L, 409L, 872L, 999L)
expect_scanner_equal(y_collection, starts, ends)

set.seed(7002)
y_left <- c(rnorm(2, 12, 0.01), rnorm(98))
left_result <- expect_scanner_equal(y_left, 1L, 100L)
stopifnot(identical(as.integer(left_result$split), 2L))

set.seed(7003)
y_central <- c(rnorm(50), rnorm(50, 8))
central_result <- expect_scanner_equal(y_central, 1L, 100L)
stopifnot(identical(as.integer(central_result$split), 50L))

set.seed(7004)
y_right <- c(rnorm(98), rnorm(2, 12, 0.01))
right_result <- expect_scanner_equal(y_right, 1L, 100L)
stopifnot(identical(as.integer(right_result$split), 98L))

set.seed(7101)
y_invariant <- rnorm(1000)
invariant_starts <- c(1L, 12L, 101L, 401L, 700L)
invariant_ends <- c(1000L, 311L, 600L, 900L, 1000L)
base_result <- expect_scanner_equal(y_invariant, invariant_starts, invariant_ends)
increasing_result <- expect_scanner_equal(
  exp(y_invariant / 4), invariant_starts, invariant_ends
)
reversed_result <- expect_scanner_equal(-y_invariant, invariant_starts, invariant_ends)
stopifnot(
  identical(as.integer(base_result$split), as.integer(increasing_result$split)),
  identical(as.integer(base_result$split), as.integer(reversed_result$split)),
  abs(base_result$value - increasing_result$value) <= 1e-12,
  abs(base_result$value - reversed_result$value) <= 1e-10
)

prunecps_buggy_reference <- function(y, cps, test, thresholds, M = 10000, d = 2) {
  minlength <- get_ns("minimumSegmentLength")(test)
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
      value <- get_ns("WBS")(
        y[(cps[i - 1L] + 1L):cps[i + 1L]], M, test, d
      )
      if (!is.finite(value$value) || value$value <= thresholds[n]) {
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

set.seed(7201)
y_prune <- rnorm(20)
prune_thresholds <- rep(-Inf, 20L)
prune_thresholds[13L] <- Inf
set.seed(7202)
buggy_result <- prunecps_buggy_reference(
  y_prune, c(4L, 10L, 13L), "lepage", prune_thresholds, M = 20L
)
set.seed(7202)
fixed_result <- get_ns("prunecps")(
  y_prune, c(4L, 10L, 13L), "lepage", prune_thresholds, M = 20L
)
stopifnot(
  identical(as.integer(buggy_result), c(4L, 13L)),
  identical(as.integer(fixed_result), 13L)
)
