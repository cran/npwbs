\name{detectChanges}
     \alias{detectChanges}
     \title{Nonparametric detection of multiple change points using Wild Binary Segmentation}
     \author{Gordon J. Ross \email{gordon.ross@ed.ac.uk}}
     \description{
       Returns the estimated number and locations of the change points in a sequence of univariate observations. For full details of how this procedure works, please see G. J. Ross (2026) - "Nonparametric Detection of Multiple Location-Scale Change Points via Wild Binary Segmentation" at \url{https://arxiv.org/abs/2107.01742}
    }
    
    \usage{
     detectChanges(y, alpha=0.05, prune=TRUE, M=1000, d=2,
                   displayOutput=FALSE, method="lepage", breakTies=TRUE,
                   combination="sum")
     }
     \arguments{
       \item{y}{A non-empty finite numeric vector to test for change points. Matrices and other dimensioned inputs are rejected.}
       \item{alpha}{Required Type I error rate. Lepage with \code{combination = "sum"} supports 0.05 and 0.01; Lepage with \code{combination = "max"} and all other methods support 0.05 only.}
       \item{prune}{Whether to prune potential excess change points via post-processing. Most likely should be left as TRUE.}
       \item{M}{Number of subsequences to sample during WBS and pruning. The default is \code{M = 1000}. Only exact tuples represented in the private calibration registry are accepted.}
       \item{d}{Endpoint trimming within each sampled interval. The existing methods, including Anderson--Darling, use \code{d = 2}; \code{method = "zhang"} uses \code{d = 4} automatically and accepts only an explicit value of 4.}
       \item{displayOutput}{A single non-missing logical. If TRUE, prints information throughout the recursive search for change points.}
       \item{method}{Rank statistic: exactly one of \code{"lepage"}, \code{"mw"}, \code{"mood"}, \code{"cvm"}, \code{"baumgartner"}, \code{"zhang"}, or \code{"anderson_darling"}.}
       \item{breakTies}{If TRUE, tied inputs receive one random strict ordering before detection; if FALSE, ties are an error.}
       \item{combination}{For \code{method = "lepage"}, exactly \code{"sum"} (the default, \eqn{Z_{MW}^2 + Z_{Mood}^2}) or \code{"max"} (\eqn{\max(Z_{MW}^2,Z_{Mood}^2)}). Other methods reject \code{"max"}.}
     }
     
     \details{
       Verified \code{M = 1000} tables are used by default. Exact frozen \code{M = 10000} tables remain available for all original methods. Unsupported method/combination/M/alpha/d tuples fail rather than reuse a mismatched table. The sampled interval minimum length is fixed at 10. Thresholds for Lepage, Mann--Whitney, Mood, Cramér--von Mises, and modified Baumgartner are provided through segment length \code{n = 10000}; longer inputs retain the package's existing behaviour of reusing the \code{n = 10000} threshold with one warning. Zhang and Anderson--Darling support lengths only through \code{n = 3000}.

       The \code{"zhang"} method is the exact mean/standard-deviation-standardised Zhang \eqn{Z_C} detector used in the paper. It retains \code{d = 4}, \code{alpha = 0.05}, strict threshold exceedance, and registered length-specific thresholds for \code{M = 1000} and \code{M = 10000}. Exact moments through \code{n = 1000} are bundled. Sequence lengths \code{n = 1001} through 3000 require one prior call to \code{\link{download_zhang_moments}}; an invalid installed artifact can be replaced with \code{download_zhang_moments(overwrite = TRUE)}. Lengths above 3000 are unsupported. Anderson--Darling supports only \code{M = 1000}, \code{alpha = 0.05}, \code{d = 2}, and lengths through 3000.
       
       If \code{prune = TRUE}, a post-processing pruning step retests merged neighbouring segments and removes changepoints that are no longer supported.
     }
     
     \value{
       A vector containing the detected changepoint locations. A returned value \code{k} denotes a split between observations \code{k} and \code{k + 1}, i.e. the change occurs after observation \code{k}.
     }
     \examples{
     \donttest{
set.seed(100)
y <- c(rnorm(30,0,1),rnorm(30,3,1), rnorm(30,0,1),rnorm(30,0,3))
detectChanges(y)
}
}
