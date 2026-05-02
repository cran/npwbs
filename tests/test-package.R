library(npwbs)

namespace <- asNamespace("npwbs")
LPthresholds05 <- get("LPthresholds05", envir = namespace)
LPthresholds01 <- get("LPthresholds01", envir = namespace)
prepareThresholds <- get("prepareThresholds", envir = namespace)

stopifnot(identical(formals(npwbs::detectChanges)$alpha, 0.05))
stopifnot(identical(formals(npwbs::detectChanges)$prune, TRUE))
stopifnot(identical(formals(npwbs::detectChanges)$M, 10000))
stopifnot(identical(formals(npwbs::detectChanges)$d, 2))

stopifnot(is.numeric(LPthresholds05), is.numeric(LPthresholds01))
stopifnot(all(is.finite(LPthresholds05[c(10, 100, 1000, 10000)])))
stopifnot(all(is.finite(LPthresholds01[c(10, 100, 1000, 10000)])))
stopifnot(all(LPthresholds01[10:10000] >= LPthresholds05[10:10000]))

msg <- NULL
thresholdsHigh <- withCallingHandlers(
  prepareThresholds(LPthresholds05, 13000),
  warning = function(w) {
    msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(!is.null(msg))
stopifnot(grepl("calibrated up to n = 10000", msg, fixed = TRUE))
stopifnot(identical(thresholdsHigh[10001], thresholdsHigh[10000]))
stopifnot(identical(thresholdsHigh[13000], thresholdsHigh[10000]))

msg <- NULL
thresholdsOrdinary <- withCallingHandlers(
  prepareThresholds(LPthresholds05, 10000),
  warning = function(w) {
    msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(is.null(msg))
stopifnot(identical(length(thresholdsOrdinary), length(LPthresholds05)))

errorMsg <- tryCatch(
  {
    npwbs::detectChanges(rnorm(20), alpha = 0.02)
    NA_character_
  },
  error = function(err) err$message
)
stopifnot(grepl("0.05 and 0.01", errorMsg, fixed = TRUE))

errorMsg <- tryCatch(
  {
    npwbs::detectChanges(rnorm(20), M = 9999)
    NA_character_
  },
  error = function(err) err$message
)
stopifnot(grepl("M=10000", errorMsg, fixed = TRUE))

errorMsg <- tryCatch(
  {
    npwbs::detectChanges(rnorm(20), d = 3)
    NA_character_
  },
  error = function(err) err$message
)
stopifnot(grepl("d=2", errorMsg, fixed = TRUE))

set.seed(1)
y <- c(rep(0, 50), rep(10, 50)) + rnorm(100, 0, 1e-4)
cps <- npwbs::detectChanges(y)
stopifnot(identical(as.integer(cps), 50L))

set.seed(2)
y <- c(rep(0, 40), rep(10, 20), rep(0, 40)) + rnorm(100, 0, 1e-4)
cps <- sort(as.integer(npwbs::detectChanges(y)))
stopifnot(identical(cps, c(40L, 60L)))
stopifnot(all(cps > 0L & cps < length(y)))

set.seed(3)
y <- c(rep(0, 30), rep(3, 30), rep(0, 40)) + rnorm(100, 0, 0.1)
cpsPruned <- npwbs::detectChanges(y, prune = TRUE)
cpsUnpruned <- npwbs::detectChanges(y, prune = FALSE)
stopifnot(length(cpsPruned) <= length(cpsUnpruned))
