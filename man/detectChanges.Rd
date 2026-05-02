\name{detectChanges}
     \alias{detectChanges}
     \title{Nonparametric detection of multiple change points using Wild Binary Segmentation}
     \author{Gordon J. Ross \email{gordon@gordonjross.co.uk}}
     \description{
       Returns the estimated number and locations of the change points in a sequence of univariate observations. For full details of how this procedure works, please see G. J. Ross (2021) - "Nonparametric Detection of Multiple Location-Scale Change Points via Wild Binary Segmentation" at \url{https://arxiv.org/abs/2107.01742}
    }
    
    \usage{
     detectChanges(y,alpha=0.05,prune=TRUE,M=10000,d=2,displayOutput=FALSE)
     }
     \arguments{
       \item{y}{The sequence to test for change points}
       \item{alpha}{Required Type I error (i.e. false positive) rate. Can be set to either 0.05 or 0.01}
       \item{prune}{Whether to prune potential excess change points via post-processing. Most likely should be left as TRUE.}
       \item{M}{Number of subsequences to sample during WBS and pruning. The built-in thresholds are calibrated for \code{M = 10000}, so other values are not supported by the packaged calibration.}
       \item{d}{Minimum number of observations between change points. The built-in thresholds currently support only \code{d = 2}.}
       \item{displayOutput}{If TRUE then will print some information while searching for change points}
     }
     
     \details{
       The built-in thresholds are precomputed for the Lepage statistic using \code{d = 2}, \code{M = 10000}, and \code{alpha} in \code{\{0.05, 0.01\}}. Thresholds are calibrated up to segment length \code{n = 10000}. For \code{n > 10000}, the implementation uses the \code{n = 10000} threshold and issues a warning. For \code{n > 1000}, the built-in thresholds use a smooth tail approximation fitted to additional Monte Carlo calibration points.
       
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
