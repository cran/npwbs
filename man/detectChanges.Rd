\name{detectChanges}
     \alias{detectChanges}
     \title{Nonparametric detection of multiple change points using Wild Binary Segmentation}
     \author{Gordon J. Ross \email{gordon.ross@ed.ac.uk}}
     \description{
       Returns the estimated number and locations of the change points in a sequence of univariate observations. For full details of how this procedure works, please see G. J. Ross (2026) - "Nonparametric Detection of Multiple Location-Scale Change Points via Wild Binary Segmentation" at \url{https://arxiv.org/abs/2107.01742}
    }
    
    \usage{
     detectChanges(y, alpha=0.05, prune=TRUE, M=10000, d=2,
                   displayOutput=FALSE, method="lepage", breakTies=TRUE)
     }
     \arguments{
       \item{y}{The sequence to test for change points}
       \item{alpha}{Required Type I error rate. Lepage supports 0.05 and 0.01; the other methods support 0.05.}
       \item{prune}{Whether to prune potential excess change points via post-processing. Most likely should be left as TRUE.}
       \item{M}{Number of subsequences to sample during WBS and pruning. The built-in thresholds are calibrated for \code{M = 10000}, so other values are not supported by the packaged calibration.}
       \item{d}{Endpoint trimming within each sampled interval. The built-in thresholds support only \code{d = 2}.}
       \item{displayOutput}{If TRUE then will print some information while searching for change points}
       \item{method}{Rank statistic: exactly one of \code{"lepage"}, \code{"mw"}, \code{"mood"}, \code{"cvm"}, or \code{"baumgartner"}.}
       \item{breakTies}{If TRUE, tied inputs receive one random strict ordering before detection; if FALSE, ties are an error.}
     }
     
     \details{
       The built-in thresholds are provided for Lepage, Mann--Whitney, Mood, Cramér--von Mises, and modified Baumgartner statistics using \code{d = 2} and \code{M = 10000}. Lepage at \code{alpha = 0.05} uses the final registered book-chapter calibration, while Lepage at \code{alpha = 0.01} retains the established package calibration. Mann--Whitney, Mood, Cramér--von Mises, and modified Baumgartner support \code{alpha = 0.05} only. The sampled interval minimum length is fixed at 10. Thresholds are provided through segment length \code{n = 10000}; longer inputs reuse the \code{n = 10000} threshold with one warning. Modified-Baumgartner values above \code{n = 3000} use the fitted extension accepted for the package rather than direct Monte-Carlo calibration.
       
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
