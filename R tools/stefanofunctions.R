######################
## Stefano's functions
##
####################

required.packages <- c("clue")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]

if(length(missing.packages)) install.packages(missing.packages)


## Quantile Normalization
quantileNormalization <- function (wd, distribution) 
{
  n <- nrow(wd)
  m <- ncol(wd)
  if (!missing(distribution)) 
    if (m != length(distribution)) 
      stop("The reference distribution has length different from the col dimension of the data matrix.")
  else distribution <- sort(distribution)
  o <- matrix(0, n, m)
  for (i in 1:n) o[i, ] <- order(wd[i, ])
  j <- 1
  tmp <- rep(0, n)
  while (j <= m) {
    if (missing(distribution)) {
      for (i in 1:n) tmp[i] <- wd[i, o[i, j]]
      value <- mean(tmp)
    }
    else value <- distribution[j]
    for (i in 1:n) wd[i, o[i, j]] <- value
    j <- j + 1
  }
  return(wd)
}


## Calinsky's function to determine optimal k value
calinsky <- function (hhc, ddist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 
                                                                    10))) 
{
  ans <- rep(0, gMax)
  uDist <- as.cl_ultrametric(hhc)
  uDist <- as.matrix(uDist)
  totalMedoid <- which.min(apply(uDist, 1, sum))
  totalSum <- sum(uDist[totalMedoid, ])
  n <- length(hhc$order)
  for (g in 2:gMax) {
    cclust <- cutree(hhc, k = g)
    groupMedoids <- rep(0, g)
    withinSum <- 0
    for (k in 1:g) {
      if (sum(cclust == k) == 1) 
        next
      tmp <- which.min(apply(uDist[cclust == k, cclust == 
                                     k], 1, sum))
      groupMedoids[k] <- which(names(cclust) == names(tmp))
      withinSum <- withinSum + sum(uDist[groupMedoids[k], 
                                         cclust == k])
    }
    betweenSum <- totalSum - withinSum
    betweenSum <- betweenSum/(g - 1)
    withinSum <- withinSum/(n - g)
    ans[g] <- betweenSum/withinSum
  }
  class(ans) <- "calinsky"
  return(ans)
}
