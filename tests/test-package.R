library(npwbs)
ns <- asNamespace("npwbs")
get_ns <- function(name) get(name, envir = ns)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

expected_formals <- c("y", "alpha", "prune", "M", "d", "displayOutput", "method", "breakTies")
stopifnot(identical(names(formals(detectChanges)), expected_formals))
stopifnot(identical(formals(detectChanges)$method, "lepage"))
stopifnot(identical(formals(detectChanges)$breakTies, TRUE))

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

warning_count <- 0L
warning_messages <- character()
tail_thresholds <- withCallingHandlers(
  get_ns("prepareThresholds")(get_ns("BaumgartnerThresholds05"), 10002L,
                              "baumgartner"),
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

continuous <- rnorm(40)
set.seed(4101)
named <- detectChanges(continuous, prune = FALSE)
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

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}
expect_error(detectChanges(continuous, method = "MW"), "method must be exactly")
expect_error(detectChanges(continuous, method = "mw", alpha = 0.01), "supports only alpha=0.05")
expect_error(detectChanges(continuous, method = "baumgartner", alpha = 0.01),
             "supports only alpha=0.05")
expect_error(detectChanges(continuous, method = "baumgartner", M = 100), "M=10000")
expect_error(detectChanges(continuous, method = "baumgartner", d = 3), "d=2")
expect_error(detectChanges(tied, breakTies = FALSE), "ties are present")
