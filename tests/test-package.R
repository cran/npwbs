library(npwbs)
ns <- asNamespace("npwbs")
get_ns <- function(name) get(name, envir = ns)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

expected_formals <- c("y", "alpha", "prune", "M", "d", "displayOutput", "method", "breakTies", "combination")
stopifnot(identical(names(formals(detectChanges)), expected_formals))
stopifnot(identical(formals(detectChanges)$method, "lepage"))
stopifnot(identical(formals(detectChanges)$breakTies, TRUE))
stopifnot(identical(formals(detectChanges)$combination, "sum"))
stopifnot(identical(formals(detectChanges)$M, 1000))
stopifnot(
  identical(names(formals(download_zhang_moments)), "overwrite"),
  identical(formals(download_zhang_moments)$overwrite, FALSE)
)

for (name in c("MWthresholds05", "Moodthresholds05", "CVMthresholds05",
               "BaumgartnerThresholds05")) {
  value <- get_ns(name)
  stopifnot(is.numeric(value), length(value) == 10000L, all(is.finite(value)))
}
stopifnot(identical(
  get_ns("BaumgartnerThresholds05")[c(3000L, 3001L, 4000L, 5000L, 7500L, 10000L)],
  c(14.26525503482025, 14.265406878639029, 14.394027935713112,
    14.490872517956545, 14.660522679273727, 14.776279053219584)
))
stopifnot(length(get_ns("LPthresholds05")) == 10000L)
stopifnot(length(get_ns("LPthresholds01")) == 10000L)
stopifnot(
  length(get_ns("ZhangThresholds05")) == 3000L,
  isTRUE(all.equal(
    get_ns("ZhangThresholds05")[3000L],
    14.056959675710399,
    tolerance = 1e-15
  ))
)

warning_count <- 0L
warning_messages <- character()
tail_thresholds <- withCallingHandlers(
  get_ns("prepareThresholds")(get_ns("selectThresholdEntry")(
    "baumgartner", "sum", 10000, 0.05, 2
  ), 10002L),
  warning = function(w) {
    warning_count <<- warning_count + 1L
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(warning_count == 1L)
stopifnot(identical(
  warning_messages,
  "Thresholds are provided through n = 10000; the n = 10000 threshold will be used for longer segments."
))
stopifnot(identical(tail_thresholds[10001:10002], rep(tail_thresholds[10000], 2L)))

y <- c(3, 1, 4, 1.5, 5, 9, 2, 6, 8, 7, 10, 0)
starts <- c(1L, 1L)
ends <- c(12L, 12L)
for (scanner in c("cpp_scan_lepage", "cpp_scan_mw", "cpp_scan_mood", "cpp_scan_cvm", "cpp_scan_baumgartner")) {
  result <- get_ns(scanner)(y, starts, ends, 2L)
  stopifnot(identical(names(result), c("value", "split")), is.finite(result$value), !is.na(result$split))
}

set.seed(1002)
y_location <- c(rnorm(30), rnorm(30, 4))
set.seed(2002)
stopifnot(identical(as.integer(detectChanges(y_location)), 30L))

set.seed(1005)
y_pruning <- c(rep(0, 30), rep(3, 30), rep(0, 40)) + rnorm(100, sd = 0.1)
set.seed(2005)
stopifnot(identical(as.integer(detectChanges(y_pruning)), c(30L, 60L)))

for (method in c("mw", "mood", "cvm", "baumgartner")) {
  set.seed(3100 + match(method, c("mw", "mood", "cvm", "baumgartner")))
  result <- detectChanges(rnorm(40), method = method, prune = FALSE)
  stopifnot(is.numeric(result))
}

old_user_data <- Sys.getenv("R_USER_DATA_DIR", unset = NA_character_)
zhang_user_data <- tempfile("npwbs-zhang-user-data-")
dir.create(zhang_user_data)
Sys.setenv(R_USER_DATA_DIR = zhang_user_data)
zhang_cache <- get_ns(".zhang_moment_cache")
rm(list = ls(zhang_cache, all.names = TRUE), envir = zhang_cache)

set.seed(4301)
zhang_boundary_data <- rnorm(1000)
set.seed(4302)
zhang_boundary <- detectChanges(
  zhang_boundary_data, method = "zhang", prune = FALSE
)
stopifnot(is.numeric(zhang_boundary))
stopifnot(
  exists("bundled", envir = zhang_cache, inherits = FALSE),
  !exists("extended", envir = zhang_cache, inherits = FALSE)
)

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}
expect_exact_error <- function(expr, expected) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), identical(message, expected))
  message
}
missing_zhang_message <- paste0(
  "The Zhang method requires the extended exact moment table for sequences ",
  "longer than 1000. Install it once by running:\n\n",
  "    npwbs::download_zhang_moments()\n\n",
  "Then rerun detectChanges()."
)
invalid_zhang_message <- paste0(
  "The installed Zhang moment table is invalid or incompatible. Reinstall ",
  "it by running:\n\n",
  "    npwbs::download_zhang_moments(overwrite = TRUE)"
)
unsupported_zhang_message <-
  "The Zhang method currently supports sequences of length at most 3000."
zhang_error_messages <- expect_exact_error(
  detectChanges(seq_len(1001), method = "zhang", prune = FALSE),
  missing_zhang_message
)
zhang_error_messages <- c(zhang_error_messages, expect_exact_error(
  detectChanges(seq_len(3001), method = "zhang", prune = FALSE),
  unsupported_zhang_message
))
expect_error(
  detectChanges(seq_len(40), method = "zhang", d = 2),
  "method='zhang' supports only d=4"
)

local_extension <- Sys.getenv("NPWBS_LOCAL_ZHANG_EXTENSION")
if (nzchar(local_extension) && file.exists(local_extension)) {
  installed_path <- get_ns(".zhang_extension_path")()
  dir.create(dirname(installed_path), recursive = TRUE)
  stopifnot(file.copy(local_extension, installed_path))
  if (exists("extended", envir = zhang_cache, inherits = FALSE)) {
    rm(list = "extended", envir = zhang_cache)
  }

  set.seed(4304)
  zhang_long_data <- rnorm(1001)
  set.seed(4305)
  zhang_long <- detectChanges(
    zhang_long_data, method = "zhang", prune = FALSE
  )
  stopifnot(is.numeric(zhang_long))

  already_installed <- capture.output(
    returned_path <- download_zhang_moments(),
    type = "message"
  )
  stopifnot(
    identical(returned_path, installed_path),
    any(grepl("already installed", already_installed, fixed = TRUE))
  )

  altered <- tempfile(fileext = ".rds")
  stopifnot(file.copy(local_extension, altered))
  connection <- file(altered, open = "r+b")
  seek(connection, where = 100L)
  byte <- readBin(connection, what = "raw", n = 1L)
  seek(connection, where = 100L)
  writeBin(as.raw(bitwXor(as.integer(byte), 1L)), connection)
  close(connection)
  stopifnot(file.copy(altered, installed_path, overwrite = TRUE))
  rm(list = "extended", envir = zhang_cache)
  zhang_error_messages <- c(
    zhang_error_messages,
    expect_exact_error(
      detectChanges(seq_len(1001), method = "zhang", prune = FALSE),
      invalid_zhang_message
    )
  )
}

stopifnot(!any(grepl(
  "Extended Zhang moments are not loaded",
  zhang_error_messages,
  fixed = TRUE
)))

if (identical(Sys.getenv("NPWBS_RUN_ONLINE_TESTS"), "true")) {
  online <- tempfile(fileext = ".rds")
  utils::download.file(
    "http://gordonjross.co.uk/zhang_zc_moments_d4_n1001_3000.rds",
    online,
    mode = "wb",
    quiet = FALSE
  )
  get_ns(".zhang_read_moment_file")(
    online, 1001L, 3000L, get_ns(".zhang_extension_sha256")
  )
}

if (is.na(old_user_data)) {
  Sys.unsetenv("R_USER_DATA_DIR")
} else {
  Sys.setenv(R_USER_DATA_DIR = old_user_data)
}

continuous <- rnorm(40)
set.seed(4101)
named <- detectChanges(continuous, prune = FALSE, M = 10000)
set.seed(4101)
positional <- detectChanges(continuous, 0.05, FALSE, 10000, 2, FALSE)
stopifnot(identical(named, positional))

tied <- rep(1:4, each = 8L)
warnings <- character()
set.seed(4201)
first <- withCallingHandlers(
  detectChanges(tied, prune = FALSE),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
set.seed(4201)
second <- suppressWarnings(detectChanges(tied, prune = FALSE))
stopifnot(length(warnings) == 1L, identical(first, second))

expect_error(detectChanges(continuous, method = "MW"), "method must be exactly")
expect_error(detectChanges(continuous, method = "mw", alpha = 0.01), "no built-in threshold calibration")
expect_error(detectChanges(continuous, method = "baumgartner", alpha = 0.01),
             "no built-in threshold calibration")
expect_error(detectChanges(continuous, method = "baumgartner", M = 100), "M must be exactly 1000 or 10000")
expect_error(detectChanges(continuous, method = "baumgartner", d = 3), "d=2")
expect_error(detectChanges(tied, breakTies = FALSE), "ties are present")
