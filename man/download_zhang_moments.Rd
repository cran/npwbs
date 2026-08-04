\name{download_zhang_moments}
\alias{download_zhang_moments}
\title{Install the extended exact Zhang moment table}
\description{
Downloads and verifies the fixed exact-moment artifact required by the
standardised Zhang detector for sequence lengths above 1000.
}
\usage{
download_zhang_moments(overwrite = FALSE)
}
\arguments{
\item{overwrite}{If \code{TRUE}, replace an existing installed artifact.}
}
\details{
Exact moments through sequence length 1000 are bundled with \pkg{npwbs}.
The \code{method = "zhang"} detector requires this one-time download for
lengths 1001 through 3000. The artifact is installed in the package-specific
user data directory returned by \code{tools::R_user_dir("npwbs", "data")}.
Use \code{overwrite = TRUE} to reinstall an invalid or incompatible artifact.
Sequence lengths above 3000 are unsupported.
}
\value{
The installed file path, returned invisibly.
}
