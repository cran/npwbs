library(npwbs)
ns <- asNamespace("npwbs")
get_ns <- function(name) get(name, envir = ns, inherits = FALSE)

expect_error_without_rng <- function(expression, pattern) {
  set.seed(19401L)
  before <- .Random.seed
  message <- tryCatch({ force(expression); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE), identical(.Random.seed, before))
}

stopifnot(
  as.character(utils::packageVersion("npwbs")) == "1.0",
  identical(formals(detectChanges)$M, 1000),
  identical(formals(detectChanges)$combination, "sum")
)

frozen_hashes <- c(
  LPthresholds05 = "ff1fbd837122d568f5f8b77359c0603251a00076e2501e2de80050522fb6e81b",
  LPthresholds01 = "871e6cb2763ac0b38a9557fdf85dde2f3d04d13b2d675e09e6fc5611cfbd8ee9",
  MWthresholds05 = "7240470c3126a9927fb332cbbcbb4066dd38cd163eef01e2dafcd149a73e1258",
  Moodthresholds05 = "a5192c806366df7c78b1bf1af3400e6c654688ff6cb02536114d879440750032",
  CVMthresholds05 = "5cbf0363b37334b4b39cf7e2e1f9abf6ee83683405a0d8fe744c58f39705e0c0",
  BaumgartnerThresholds05 = "ca7dc80c7192ff64c9ca458b3b8aac1f375336e7aec73e4e6e194eae4b368fc4",
  ZhangThresholds05 = "dc5778f9955d215b80d1bfd7f1d9896b49efcc270675e48452828cdfa15034e3"
)
stopifnot(all(vapply(names(frozen_hashes), function(name) {
  identical(digest::digest(get_ns(name), algo = "sha256", serialize = TRUE),
            unname(frozen_hashes[[name]]))
}, logical(1L))))

registry <- get_ns(".thresholdRegistry")()
required_keys <- c(
  "lepage|sum|1000|0.05|2", "lepage|sum|1000|0.01|2",
  "lepage|max|1000|0.05|2", "lepage|max|10000|0.05|2",
  "lepage|sum|10000|0.05|2", "lepage|sum|10000|0.01|2",
  "mw|scalar|1000|0.05|2", "mood|scalar|1000|0.05|2",
  "cvm|scalar|1000|0.05|2", "baumgartner|scalar|1000|0.05|2",
  "zhang|scalar|1000|0.05|4", "anderson_darling|scalar|1000|0.05|2",
  "mw|scalar|10000|0.05|2", "mood|scalar|10000|0.05|2",
  "cvm|scalar|10000|0.05|2", "baumgartner|scalar|10000|0.05|2",
  "zhang|scalar|10000|0.05|4"
)
stopifnot(identical(sort(names(registry)), sort(required_keys)))

continuous <- c(8, 1, 14, 4, 11, 2, 16, 6, 10, 3, 15, 7, 13, 5, 12, 9)
expect_error_without_rng(detectChanges(continuous, M = 999), "M must be exactly")
for (nearby_alpha in c(0.05 - 1e-12, 0.05 + 1e-12,
                       0.01 - 1e-12, 0.01 + 1e-12)) {
  expect_error_without_rng(
    detectChanges(continuous, alpha = nearby_alpha),
    "alpha must be exactly 0.05 or 0.01"
  )
}
expect_error_without_rng(detectChanges(matrix(continuous, ncol = 1L)),
                         "finite numeric vector")
expect_error_without_rng(detectChanges(continuous, displayOutput = NA),
                         "displayOutput must be TRUE or FALSE")
expect_error_without_rng(detectChanges(continuous, displayOutput = c(TRUE, FALSE)),
                         "displayOutput must be TRUE or FALSE")
expect_error_without_rng(detectChanges(continuous, method = "mw", alpha = 0.01),
                         "no built-in threshold calibration")
expect_error_without_rng(detectChanges(continuous, method = "mw", combination = "max"),
                         "supported only for method='lepage'")
expect_error_without_rng(detectChanges(continuous, method = "lepage", combination = "max", alpha = 0.01),
                         "no built-in threshold calibration")
for (method in c("lepage", "mw", "mood", "cvm", "baumgartner")) {
  expect_error_without_rng(detectChanges(continuous, method = method, d = 3),
                           "d=2")
}
expect_error_without_rng(detectChanges(continuous, method = "zhang", d = 2),
                         "supports only d=4")

local({
  old_options <- options(digits = 1)
  on.exit(options(old_options), add = TRUE)
  stopifnot(
    identical(get_ns(".thresholdKey")("lepage", "sum", 1000, 0.05, 2),
              "lepage|sum|1000|0.05|2"),
    identical(get_ns(".thresholdKey")("lepage", "sum", 1000, 0.01, 2),
              "lepage|sum|1000|0.01|2")
  )
  set.seed(19402L)
  digits_result <- detectChanges(continuous, alpha = 0.05, prune = FALSE)
  options(old_options)
  set.seed(19402L)
  default_digits_result <- detectChanges(continuous, alpha = 0.05, prune = FALSE)
  stopifnot(identical(digits_result, default_digits_result))
})

mw_reference <- function(y, d = 2L) {
  ranks <- rank(y)
  n <- length(y)
  out <- rep(-Inf, n)
  for (split in d:(n - d)) {
    n1 <- split; n2 <- n - split
    u <- sum(ranks[seq_len(split)]) - n1 * (n1 + 1) / 2
    out[[split]] <- (u - n1 * n2 / 2)^2 / (n1 * n2 * (n + 1) / 12)
  }
  out
}
mood_reference <- function(y, d = 2L) {
  ranks <- (rank(y) - (length(y) + 1) / 2)^2
  n <- length(y)
  out <- rep(-Inf, n)
  for (split in d:(n - d)) {
    n1 <- split; n2 <- n - split
    z <- sum(ranks[seq_len(split)]) - n1 * (n^2 - 1) / 12
    out[[split]] <- z^2 / (n1 * n2 * (n + 1) * (n^2 - 4) / 180)
  }
  out
}
cvm_reference <- function(y, d = 2L) {
  ranks <- rank(y)
  n <- length(y)
  out <- rep(-Inf, n)
  for (split in d:(n - d)) {
    m <- split; q <- n - m
    left_cdf <- cumsum(tabulate(ranks[seq_len(m)], nbins = n)) / m
    right_cdf <- cumsum(tabulate(ranks[-seq_len(m)], nbins = n)) / q
    w <- m * q / n^2 * sum((left_cdf - right_cdf)^2)
    mean <- (n + 1) / (6 * n)
    variance <- (n + 1) *
      (4 * m * q * n - 3 * m^2 - 3 * q^2 - 2 * m * q) /
      (180 * m * q * n^2)
    out[[split]] <- (w - mean) / sqrt(variance)
  }
  out
}
baumgartner_reference <- function(y, d = 2L) {
  ranks <- rank(y)
  n <- length(y)
  component <- function(group_ranks, other_size) {
    group_ranks <- sort(group_ranks)
    size <- length(group_ranks)
    i <- seq_len(size)
    p <- i / (size + 1)
    expected <- (n + 1) * i / (size + 1)
    variance_scale <- other_size * (n + 1) / (size + 2)
    mean((group_ranks - expected)^2 / (p * (1 - p) * variance_scale))
  }
  out <- rep(-Inf, n)
  for (split in d:(n - d)) {
    out[[split]] <- 0.5 * (
      component(ranks[seq_len(split)], n - split) +
      component(ranks[-seq_len(split)], split)
    )
  }
  out
}
ad_variance_reference <- function(n, split) {
  harmonic <- c(0, cumsum(1 / seq_len(n - 1L)))
  g <- 0
  if (n > 2L) {
    for (current in 2:(n - 1L)) {
      g <- g + 1 / current^2 +
        2 * (harmonic[[current]] - 1) / (current * (current + 1))
    }
  }
  m <- split; q <- n - m; H <- 1 / m + 1 / q
  h <- harmonic[[n]]
  a <- 4 * g - 6 + (10 - 6 * g) * H
  b <- 12 * g + 8 * h - 22 + (2 * g - 14 * h - 4) * H
  cc <- 36 * h + 4 + (2 * h - 6) * H
  (a * n^3 + b * n^2 + cc * n + 24) /
    ((n - 1) * (n - 2) * (n - 3))
}
anderson_darling_reference <- function(y, d = 2L) {
  ranks <- rank(y)
  n <- length(y)
  out <- rep(-Inf, n)
  for (split in d:(n - d)) {
    cumulative_left <- cumsum(tabulate(ranks[seq_len(split)], nbins = n))
    j <- seq_len(n - 1L)
    raw <- sum((n * cumulative_left[j] - split * j)^2 /
                 (j * (n - j))) / (split * (n - split))
    out[[split]] <- (raw - 1) / sqrt(ad_variance_reference(n, split))
  }
  out
}
zhang_reference <- function(y, d = 4L) {
  ranks <- rank(y)
  n <- length(y)
  artifact <- get_ns(".zhang_moment_artifacts")(n)$bundled
  offsets <- artifact$data$offsets
  out <- rep(-Inf, n)
  weights <- function(size) log(size / (seq_len(size) - 0.5) - 1)
  pooled_weights <- weights(n)
  for (split in d:(n - d)) {
    left <- sort(ranks[seq_len(split)])
    right <- sort(ranks[-seq_len(split)])
    zc <- (sum(weights(split) * pooled_weights[left]) +
             sum(weights(n - split) * pooled_weights[right])) / n
    m_star <- min(split, n - split)
    position <- offsets[[n - artifact$metadata$min_N + 1L]] + m_star - 4L + 1L
    mean <- artifact$data$mean_zplus[[position]]
    sd <- artifact$data$sd_zplus[[position]]
    out[[split]] <- ((4 - zc) - mean) / sd
  }
  out
}

expect_reference_scan <- function(scanner, reference, y, d,
                                  start = 1L, end = length(y)) {
  local_y <- y[start:end]
  expected_values <- reference(local_y, d)
  expected_local_split <- which.max(expected_values)
  expected_split <- as.integer(start + expected_local_split - 1L)
  arguments <- list(y, as.integer(start), as.integer(end), as.integer(d))
  if (identical(scanner, "cpp_scan_zhang")) {
    moments <- get_ns(".zhang_moment_artifacts")(length(local_y))
    arguments <- c(arguments, list(moments$bundled, moments$extended))
  }
  observed <- do.call(get_ns(scanner), arguments)
  stopifnot(
    isTRUE(all.equal(observed$value, expected_values[[expected_local_split]],
                     tolerance = 1e-10)),
    identical(as.integer(observed$split), as.integer(expected_split))
  )
}

numerical_cases <- list(
  c(7, 9, 4, 8, 6, 3, 1, 5, 0, 2),
  c(8, 1, 14, 4, 11, 2, 16, 6, 10, 3, 15, 7, 13, 5, 12, 9),
  seq_len(16L),
  rev(seq_len(16L))
)
for (case in numerical_cases) {
  expect_reference_scan("cpp_scan_mw", mw_reference, case, 2L)
  expect_reference_scan("cpp_scan_mood", mood_reference, case, 2L)
  expect_reference_scan("cpp_scan_cvm", cvm_reference, case, 2L)
  expect_reference_scan("cpp_scan_baumgartner", baumgartner_reference, case, 2L)
  expect_reference_scan("cpp_scan_anderson_darling", anderson_darling_reference, case, 2L)
  if (length(case) >= 8L) {
    expect_reference_scan("cpp_scan_zhang", zhang_reference, case, 4L)
  }
  if (length(case) >= 16L) {
    expect_reference_scan("cpp_scan_mw", mw_reference, case, 2L, 3L, 14L)
    expect_reference_scan("cpp_scan_mood", mood_reference, case, 2L, 3L, 14L)
    expect_reference_scan("cpp_scan_cvm", cvm_reference, case, 2L, 3L, 14L)
    expect_reference_scan("cpp_scan_baumgartner", baumgartner_reference,
                          case, 2L, 3L, 14L)
    expect_reference_scan("cpp_scan_anderson_darling",
                          anderson_darling_reference, case, 2L, 3L, 14L)
    expect_reference_scan("cpp_scan_zhang", zhang_reference,
                          case, 4L, 3L, 14L)
  }
}

starts <- c(1L, 2L, 5L)
ends <- c(16L, 13L, 16L)
for (combination in c("sum", "max")) {
  references <- Map(function(start, end) {
    local <- continuous[start:end]
    mw <- mw_reference(local)
    mood <- mood_reference(local)
    values <- if (combination == "sum") mw + mood else pmax(mw, mood)
    split <- which.max(values)
    list(value = values[[split]], split = as.integer(start + split - 1L))
  }, starts, ends)
  best <- which.max(vapply(references, `[[`, numeric(1L), "value"))
  observed <- get_ns("cpp_scan_lepage")(continuous, starts, ends, 2L, combination)
  stopifnot(
    isTRUE(all.equal(observed$value, references[[best]]$value, tolerance = 1e-12)),
    identical(as.integer(observed$split), references[[best]]$split)
  )
}

set.seed(19501L)
sum_default <- detectChanges(continuous, M = 1000, prune = FALSE)
set.seed(19501L)
sum_explicit <- detectChanges(continuous, M = 1000, prune = FALSE, combination = "sum")
stopifnot(identical(sum_default, sum_explicit))

set.seed(19502L)
recursive_output <- capture.output(get_ns("detectChanges_aux")(
  rep(seq_len(20L), 2L), 1L, 40L, "lepage", rep(-Inf, 40L),
  M = 1L, d = 2L, displayOutput = TRUE
))
stopifnot(length(recursive_output) > 1L)

ad_case <- c(7, 9, 4, 8, 6, 3, 1, 5, 0, 2)
ad_scan <- get_ns("cpp_scan_anderson_darling")(ad_case, 1L, 10L, 2L)
ad_threshold <- get_ns("AndersonDarlingThresholds05M1000")[[10L]]
set.seed(19503L)
ad_raw <- detectChanges(
  ad_case, method = "anderson_darling", M = 1000,
  alpha = 0.05, d = 2, prune = FALSE
)
set.seed(19503L)
ad_pruned <- detectChanges(
  ad_case, method = "anderson_darling", M = 1000,
  alpha = 0.05, d = 2, prune = TRUE
)
ad_expected <- if (ad_scan$value > ad_threshold) ad_scan$split else numeric()
stopifnot(
  identical(length(ad_raw), length(ad_expected)),
  identical(as.integer(ad_raw), as.integer(ad_expected)),
  identical(length(ad_pruned), length(ad_expected)),
  identical(as.integer(ad_pruned), as.integer(ad_expected))
)

boundary_thresholds <- rep(Inf, 10L)
boundary_thresholds[[10L]] <- ad_scan$value
set.seed(19504L)
stopifnot(length(get_ns("detectChanges_aux")(
  ad_case, 1L, 10L, "anderson_darling", boundary_thresholds,
  M = 1L, d = 2L
)) == 0L)
boundary_thresholds[[10L]] <- ad_scan$value -
  .Machine$double.eps * max(1, abs(ad_scan$value))
set.seed(19504L)
stopifnot(identical(as.integer(get_ns("detectChanges_aux")(
  ad_case, 1L, 10L, "anderson_darling", boundary_thresholds,
  M = 1L, d = 2L
)), as.integer(ad_scan$split)))

expect_error_without_rng(
  detectChanges(ad_case, method = "anderson_darling", M = 10000),
  "no built-in threshold calibration"
)
expect_error_without_rng(
  detectChanges(ad_case, method = "anderson_darling", alpha = 0.01),
  "no built-in threshold calibration"
)
expect_error_without_rng(
  detectChanges(ad_case, method = "anderson_darling", d = 3),
  "d=2"
)
expect_error_without_rng(
  detectChanges(seq_len(3001L), method = "anderson_darling"),
  "supports sequences of length at most 3000"
)

exports <- getNamespaceExports("npwbs")
stopifnot(
  identical(sort(exports), sort(c("detectChanges", "download_zhang_moments"))),
  !any(grepl("^bs|uniform", exports, ignore.case = TRUE))
)
