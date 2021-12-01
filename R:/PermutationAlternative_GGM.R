#rm(list = ls())
library("readr")
library("networktree")
library("MASS")
library("strucchange")
library("modelr")
library("tidyverse")
library("parallel")
library("qgraph")
data(big5)
# ----------------------------------------------------------------------
# Load helper function
# resampling the dataset 
permutationData <- function(origdata){
  n <- nrow(origdata)
  index <- sample(1:n, n, replace = FALSE)
  newdata <- origdata
  newdata[, 1] <- origdata[index, 1] 
  return(newdata)
}

# ------------------------------------------------------------
# Obtain p-value trough permutation

permutationPval <- function(nodevars, splitvar, permiter = 1e4){ 
  
  estim <- FALSE
  while(estim == FALSE){
    nodevars <- as.data.frame(nodevars)
    splitvars <- as.data.frame(splitvar)
    origdata <- as.data.frame(cbind(splitvar, nodevars))
    #MOB direct
    p <- ncol(nodevars)
    n <- nrow(nodevars)
    minsplit <- p*(p-1)/2 + 1
    if(minsplit > n/2) {minsplit <- (n/2 -10)}
    ## MOB
    out <-  try(networktree::networktree(nodevars = nodevars, splitvars = splitvars, method = c("mob"), model = "correlation",
                                        verbose = F, minsplit = minsplit, maxdepth = 2))

    if(!inherits(out, "try-error")) {
      estim <- TRUE
      ## MOB
      orig_stat <- unlist(out)$node.info.test1
    } else {
      warning("MOB for GGM could not be calculated. Check the input nodevars and splitvar.")
    }
  }
  # Resample from that original dataset
  out_stat <- numeric(permiter)
  
  for(f in 1:permiter) {
    run <- FALSE
    while(run == FALSE){
      permdata <- permutationData(origdata)
      permsplitvar <- permdata[, 1]
      permsplitvars <- as.data.frame(permsplitvar)
      
      out <-  try(networktree::networktree(nodevars = nodevars, splitvars = permsplitvars, method = c("mob"), model = "correlation",
                                           verbose = F, minsplit = minsplit, maxdepth = 2))
      if(!inherits(out, "try-error")) run <- TRUE
    }
    #out_stat[f] <- out
    out_stat[f] <- unlist(out)$node.info.test1
  }
  
  # Determine p-value
  pval <- mean(out_stat >= orig_stat)
  
  if(pval == 0){
    pval <- 1/permiter
  }
  
  res <- list(pval = pval, testStat = orig_stat)
  
  return(res)
}

# For user: Specify the variables of your network (nodevars) and 
# the auxiliary variable to test against (splitvar)
permutationPval(nodevars = big5[, 1:5], splitvar = big5[, 6], permiter = 1e4)
