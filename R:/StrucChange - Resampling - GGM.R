# rm(list = ls())
library("readr")
library("networktree")
library("MASS")
library("strucchange")
library("modelr")
library("tidyverse")
library("parallel")

# ----------------------------------------------------------------------

# resampling the dataset 
permutationData <- function(origdata){
  n <- nrow(origdata)
  index <- sample(1:n, n, replace = FALSE)
  newdata <- origdata
  newdata[, 2] <- origdata[index, 2] 
  return(newdata)
}


# ----------------------------------------------------------------------
# Computing structural change test

strucchangeGGM <- function(process, splitvar, model = "all"){
  
  k <- NCOL(process)
  n <- NROW(process)
  nobs <- n 
  
  ## scale process
  process <- process/sqrt(nobs)
  meat <- crossprod(process)
  J12 <- root.matrix(chol2inv(chol(meat)))
  
  process <- t(J12 %*% t(process))  
  
  ## Order along splitting variable
  zi <- as.matrix(splitvar)
  oi <- order(zi)
  proci <- process[oi, , drop = FALSE]
  
  ## order partitioning variable
  zi <- zi[oi]
  # re-apply factor() added to drop unused levels
  zi <- factor(zi, levels = unique(zi))
  
  segweights <- table(zi)
  segweights <- as.vector(segweights)/nobs
  
  proci <- apply(proci, 2L, cumsum)
  tt0 <- head(cumsum(table(zi)), -1L)
  tt <- head(cumsum(segweights), -1L)
  
  proci <- rowSums(proci^2)
  stat <- max(proci[tt0] / (tt * (1-tt)))
  
  return(stat)
}

## Test function
# Generate a dataset
# origdata <- data.frame(id = 1:500, GGMSimulationSplit(p = 5, n = 500, prob = .2, delta_interaction = 0))
# nodevars <- as.data.frame(origdata[, 3:7])
# splitvar <- as.data.frame(origdata[, 2])
# splitvars <- as.data.frame(as.factor(origdata[, 2]))
# 
# strucchangeGGM(nodevars, splitvar, model = "correlation")


# ------------------------------------------------------------
# Obtain p-value trough permutation

permutationPval <- function(i, p, n, deltacor = 0, permiter){ 
  
  estim <- FALSE
  while(estim == FALSE){
    # Generate Original Dataset
    
    origdata <- data.frame(id = 1:n, 
                           GGMSimulationSplit(p = p, n = n, prob = .2, 
                                              delta_interaction = deltacor))
    nodevars <- as.data.frame(origdata[, 3:(p+2)])
    process <- networktree::mvnfit(nodevars, estfun = TRUE)$estfun
    
    splitvar <- as.data.frame(origdata[, 2])
    
    orig_stat <- try(strucchangeGGM(process, splitvar), silent = TRUE)
    if(!inherits(orig_stat, "try-error")) estim <- TRUE
  }
  # Resample from that original dataset
  out_stat <- numeric(permiter)
  
  for(i in 1:permiter) {
    run <- FALSE
    while(run == FALSE){
      permdata <- permutationData(origdata)
      permsplitvar <- as.data.frame(permdata[, 2])
      
      out <- try(strucchangeGGM(process, splitvar = permsplitvar), silent = TRUE)
      if(!inherits(out, "try-error")) run <- TRUE
    }
    out_stat[i] <- out
  }
  
  # Determine p-value
  pval <- mean(out_stat > orig_stat)
  
  res <- list(p = p, n = n, cor = deltacor, pval = pval)
  
  return(res)
}

#permutationPval(1, p = 5, n = 250, bootiter  = 100)

#----------------------------------------------------------------------
# Simulations Regenerate the p-value using resampling

n <- c(200, 500, 2000) # sample size
p <- c(5, 10, 15) # size of network
repiter <- 5000 # number of repetitions 
permiter <- 1000 # number of samples for permutation approach
cor <- 0 # measurement invariance violation (set to 0 to simulate data under the null - no measurement invariance violation)
numCores <- detectCores()
cntr <- 0
res_GGMPerm <- list()

for (pi in 1:length(p)) {
  for (ni in 1:length(n)){
    for (ci in 1:length(cor)) {
      cntr <- cntr + 1 
      pit <- p[pi]
      nit <- n[ni]
      cit <- cor[ci]
      res_GGMPerm[[cntr]] <- mclapply(1:repiter, permutationPval, p = pit, n = nit,
                                   deltacor = cit,
                                   permiter = permiter, mc.cores = numCores)
    }
  }
}

out <- as.numeric(unlist(res_GGMPerm))
res <- na.omit(out)
pvalPermutationGGM <- as.data.frame(matrix(res, ncol = 4, byrow = T))
colnames(pvalPermutationGGM) <- c("p", "n", "cor", "pval")
write.csv(pvalPermutationGGM, "pval-Permutation-GGM-cont.csv")

