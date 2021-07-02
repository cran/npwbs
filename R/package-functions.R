
mysample <- function(x,size,replace=FALSE,prob=NULL) {
  if (length(x)==1) {x <- c(x,x)}
  sample(x,size,replace,prob)
}

mannWhitneySegment <- function(y,d=2,doabs=TRUE) {
  ranks <- rank(y)
  n <- length(y)
  vals <- numeric(n)
  r1 <- 0
  if (d > 1) {
    r1 <- sum(ranks[1:(d-1)])
  }

  for (i in d:(n-d+1)) {
    n1 <- i;
    n2 <- n-n1
    r1 <- r1 + ranks[i]
    U<- r1- n1*(n1+1)/2
    mu <- n1*n2 / 2
    sigma <- sqrt( (n1*n2*(n1+n2+1))/12)

    vals[i] <- (U-mu)/sigma
    if (doabs==TRUE) {
      vals[i] <- abs(vals[i])
    }
  }
  if (doabs==TRUE) {
    if (d > 1) {vals[1:(d-1)] <- -Inf}
    vals[(n-d+1):n] <- -Inf
    } else {
    if (d > 1) {vals[1:(d-1)] <- Inf}
    vals[(n-d+1):n] <- Inf
  }
  vals
}

moodSegment <- function(y,d=2,doabs=TRUE) {
  ranks <- rank(y)
  n <- length(y)
  ranks <- (ranks - (n+1)/2)^2
  vals <- numeric(n)
  r1 <- 0
  if (d > 1) {
    r1 <- sum(ranks[1:(d-1)])
  }
  for (i in d:(n-d+1)) {
    n1 <- i; n2 <- n-n1
    r1 <- r1 + ranks[i]
    U <- r1
    mu <- n1*(n^2-1)/12
    sigma <- sqrt( (n1*n2*(n+1)*(n^2-4))/180)

    vals[i] <- (U-mu)/sigma
    if (doabs==TRUE) {
      vals[i] <- abs(vals[i])
    }
  }
  if (doabs==TRUE) {
    if (d > 1) {vals[1:(d-1)] <- -Inf}
    vals[(n-d+1):n] <- -Inf
  } else {
    if (d > 1) {vals[1:(d-1)] <- Inf}
    vals[(n-d+1):n] <- Inf
  }
  vals
}

lepageSegment <- function(y,d=2) {
   vals <- (mannWhitneySegment(y,d,doabs=TRUE)^2 + moodSegment(y,d,doabs=TRUE)^2)
   #mu <- 2; sigma <- sqrt(2*2) #this is chi squared, so rescale it
  # vals <- (vals-mu)/sigma
   n <- length(y)
   if (d > 1) {vals[1:(d-1)] <- -Inf}
   vals[(n-d+1):n] <- -Inf
   vals
}

WBS <- function(y,sims=10000,test="MW",d=2) {
  vals <- numeric(sims)
  starts <- numeric(sims)
  ends <- numeric(sims)
  n <- length(y)
  minlength <- 8
  if (test=="LP") {minlength=10}
  for (i in 1:sims) {
    start <- mysample(1:(length(y)-minlength),1)
    end <- mysample( (start+minlength-1):n,1)
    if (test=="MW") {
      vals[i] <- max(mannWhitneySegment(y[start:end],d=2))
    } else if (test=="LP") {
      vals[i] <- max(lepageSegment(y[start:end],d=2))
    }
    starts[i] <- start
    ends[i] <- end
  }
  list(statistics=vals,starts=starts,ends=ends)
}


detectChanges <- function(y,alpha=0.05,prune=TRUE,M=10000,d=2,displayOutput=FALSE) {
  test <- "LP"
  if (alpha != 0.05 && alpha != 0.01) {
    stop("Error: only supported values for alpha are 0.05 and 0.01")
  }

  thresholds <- NULL
  if (alpha == 0.05) {
    thresholds <- LPthresholds05
  } else if (alpha == 0.01) {
    thresholds <- LPthresholds01
  }

  if (length(thresholds) < length(y)) {
    warning("Warning: length(y) is longer than maximum precomputed threshold n. Procedure will still work, but false positives may be slightly inflated")
    z <-length(thresholds)
    thresholds <- c(thresholds, rep(thresholds[z], length(y)-z))
  }

  cps <-  detectChanges_aux(y,1,length(y),test,thresholds,M,d,displayOutput)
  if (prune==TRUE) {
    cps <- prunecps(y,cps,test,thresholds,M=10000,d=2)
  }
  return(cps)
}



detectChanges_aux <- function(y,start,end,test,thresholds,M=10000,d=2,displayOutput=FALSE) {
  minlength <- NA
  if (test=="MW") {minlength=8}
  if (test=="LP") {minlength=10}

  n <- end-start+1

  if (n < minlength) {return(numeric())}

  val <-WBS(y[start:end],M,test,d)

  threshold <- thresholds[n]
  if (max(val$statistics) <= threshold) { #need absolutely greater than
    return(numeric())
  } else {
    ind <- which.max(val$statistics)
    thisstart <- start  + val$starts[ind] -1
    thisend <- start + val$ends[ind] -1
    z <- y[thisstart:thisend]
    k <- NA
    if (test=="MW") {k <- which.max(mannWhitneySegment(z))}
    if (test=="LP") {k <- which.max(lepageSegment(z))}
    k <- k + thisstart  -1 #check indexing
    #k <- k + start-1

    if (displayOutput==TRUE) {
      print(sprintf("found change point at %s on [%s,%s] with statistic %s",k,start,end, max(val$statistics)))
    }
    return(c(detectChanges_aux(y,start,k,test,thresholds,M,d), k, detectChanges_aux(y,k+1,end,test,thresholds,M,d)))
  }
}

prunecps <- function(y,cps,test,thresholds,M=10000,d=2) {
  cps <- c(0,cps,length(y))

  while (TRUE) {
    deleted <- FALSE
    if (length(cps) < 3) {break}

    i <- 2
    while(TRUE) {
     # print(i)
      if (i >= length(cps)) {break}
      n <- cps[i+1] - cps[i-1]

      if ( (test=="LP" && n < 10) || (test=="LP" && n < 8)) {
      	cps <- cps[-i]; next #sequence is too short to find any chnages
      }

      val <- WBS(y[(cps[i-1]+1):cps[i+1] ],M,test,d)
      if (max(val$statistics) <= thresholds[n]) {
        #print(sprintf("removed %s",cps[i]))
        cps <- cps[-i]
        deleted <- TRUE
      } else {
        i <- i+1
      }
    }
    if (deleted==FALSE) {
      break
    }
  }
  return(cps[-c(1,length(cps))])
}
