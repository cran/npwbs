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
       \item{M}{Number of subsequences to sample during WBS. Should be left as M=10000}
       \item{d}{Minimum number of observations between change points. Should be left as d=2.}
       \item{displayOutput}{If TRUE then will print some information while searching for change points}
     }
     
     \value{
       A vector containing the location of the detected change points
     }
     \examples{
     \donttest{
set.seed(100)
y <- c(rnorm(30,0,1),rnorm(30,3,1), rnorm(30,0,1),rnorm(30,0,3))
detectChanges(y)
}
}


