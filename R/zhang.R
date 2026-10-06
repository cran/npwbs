.zhang_moment_cache <- new.env(parent = emptyenv())
.zhang_extension_filename <- "zhang_zc_moments_d4_n1001_3000.rds"
.zhang_extension_sha256 <- "5fdf35c3649b0182241c4d29962267eea66acb8821b292874f754583877bbe55"

.zhang_extension_path <- function() {
  file.path(tools::R_user_dir("npwbs", "data"), .zhang_extension_filename)
}

.zhang_missing_moments <- function() {
  stop(
    paste0(
      "The Zhang method requires the extended exact moment table for sequences ",
      "longer than 1000. Install it once by running:\n\n",
      "    npwbs::download_zhang_moments()\n\n",
      "Then rerun detectChanges()."
    ),
    call. = FALSE
  )
}

.zhang_invalid_moments <- function() {
  stop(
    paste0(
      "The installed Zhang moment table is invalid or incompatible. Reinstall ",
      "it by running:\n\n",
      "    npwbs::download_zhang_moments(overwrite = TRUE)"
    ),
    call. = FALSE
  )
}

.zhang_read_moment_file <- function(path, min_N, max_N, sha256 = NULL) {
  if (!file.exists(path)) stop("Zhang moment table does not exist: ", path)
  if (!is.null(sha256)) {
    observed <- digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
    if (!identical(observed, sha256)) {
      stop("Zhang moment table checksum validation failed.")
    }
  }
  artifact <- readRDS(path)
  metadata <- artifact$metadata
  data <- artifact$data
  expected_N <- seq.int(as.integer(min_N), as.integer(max_N))
  expected_offsets <- as.integer(c(
    0L, cumsum(floor(expected_N / 2L) - 3L)
  ))
  expected_size <- expected_offsets[[length(expected_offsets)]]
  valid_metadata <- is.list(metadata) &&
    identical(metadata$schema_version, 1L) &&
    identical(metadata$statistic, "zhang_zc") &&
    identical(metadata$standardised_method, "zhang_zc_z") &&
    identical(metadata$d, 4L) &&
    identical(metadata$min_N, as.integer(min_N)) &&
    identical(metadata$max_N, as.integer(max_N)) &&
    identical(
      metadata$moment_definition_version,
      "zhang_zc_finite_permutation_dp_v1"
    ) &&
    identical(
      metadata$storage,
      "flat non-redundant m=4:floor(N/2) vectors with zero-based offsets"
    )
  valid_data <- is.list(data) &&
    identical(names(data), c("N", "offsets", "mean_zplus", "sd_zplus")) &&
    identical(data$N, expected_N) &&
    identical(data$offsets, expected_offsets) &&
    is.numeric(data$mean_zplus) &&
    is.numeric(data$sd_zplus) &&
    length(data$mean_zplus) == expected_size &&
    length(data$sd_zplus) == expected_size &&
    all(is.finite(data$mean_zplus)) &&
    all(is.finite(data$sd_zplus)) &&
    all(data$sd_zplus > 0)
  if (!valid_metadata || !valid_data) {
    stop("Zhang moment table metadata validation failed.")
  }
  artifact
}

.zhang_bundled_moments <- function() {
  if (exists("bundled", envir = .zhang_moment_cache, inherits = FALSE)) {
    return(get("bundled", envir = .zhang_moment_cache, inherits = FALSE))
  }
  path <- system.file(
    "extdata", "zhang_zc_moments_d4_n0010_1000.rds", package = "npwbs"
  )
  artifact <- .zhang_read_moment_file(path, 10L, 1000L)
  assign("bundled", artifact, envir = .zhang_moment_cache)
  artifact
}

.zhang_extended_moments <- function() {
  if (exists("extended", envir = .zhang_moment_cache, inherits = FALSE)) {
    return(get("extended", envir = .zhang_moment_cache, inherits = FALSE))
  }
  path <- .zhang_extension_path()
  if (!file.exists(path)) .zhang_missing_moments()
  artifact <- tryCatch(
    .zhang_read_moment_file(
      path, 1001L, 3000L, .zhang_extension_sha256
    ),
    error = function(e) .zhang_invalid_moments()
  )
  assign("extended", artifact, envir = .zhang_moment_cache)
  artifact
}

.zhang_moment_artifacts <- function(n) {
  if (n > 3000L) {
    stop(
      "The Zhang method currently supports sequences of length at most 3000.",
      call. = FALSE
    )
  }
  list(
    bundled = .zhang_bundled_moments(),
    extended = if (n > 1000L) .zhang_extended_moments() else NULL
  )
}

.zhang_moment_lookup <- function(N, m) {
  artifacts <- .zhang_moment_artifacts(N)
  artifact <- if (N <= 1000L) artifacts$bundled else artifacts$extended
  position <- N - artifact$metadata$min_N + 1L
  m_star <- min(m, N - m)
  index <- artifact$data$offsets[position] + m_star - 3L
  list(
    mean_zplus = artifact$data$mean_zplus[index],
    sd_zplus = artifact$data$sd_zplus[index]
  )
}

download_zhang_moments <- function(overwrite = FALSE) {
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE")
  }
  path <- .zhang_extension_path()
  if (file.exists(path) && !overwrite) {
    artifact <- tryCatch(
      .zhang_read_moment_file(
        path, 1001L, 3000L, .zhang_extension_sha256
      ),
      error = function(e) NULL
    )
    if (!is.null(artifact)) {
      assign("extended", artifact, envir = .zhang_moment_cache)
      message("Zhang moment table is already installed.")
      return(invisible(path))
    }
  }

  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("zhang-moments-", tmpdir = directory, fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)

  message("Downloading Zhang moment table (approximately 19.3 MB)...")
  utils::download.file(
    "https://gordonjross.co.uk/zhang_zc_moments_d4_n1001_3000.rds",
    temporary,
    mode = "wb",
    quiet = FALSE
  )
  artifact <- .zhang_read_moment_file(
    temporary, 1001L, 3000L, .zhang_extension_sha256
  )
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace the existing Zhang moment table.")
  }
  if (!file.rename(temporary, path)) {
    stop("Could not install the Zhang moment table.")
  }
  assign("extended", artifact, envir = .zhang_moment_cache)
  message("Zhang moment table installed.")
  invisible(path)
}
