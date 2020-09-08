# rm(list = ls())
library("readr")
library("networktree")
library("MASS")
library("strucchange")
library("modelr")
library("tidyverse")
library("parallel")
source("~/Documents/PhD/03_Programming/MOBSimulation/R/DataSimulation.R")

# ----------------------------------------------------------------------
# Bootstrap Function
bootstrapData <- function(coef, obs, par){
  
  mean <- coef[1:par] 
  sigma <- matrix(0, ncol = par, nrow = par)
  diag(sigma) <- coef[(par+1) : (par+par)] 
  sigma[lower.tri(sigma)] <- coef[(2*par+1):length(coef)] 
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)] 
  data <- rmvnorm(n, mean, sigma,
                     method=c("chol"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)
  z1 <- sort(sample(18:75, size = n, replace = TRUE)) 
  newdata <- cbind(as.data.frame(z1), data)
  return(newdata)
}

#bootstrapData(coef = coef, obs = 200, par = 5)

permutationData <- function(origdata){
  n <- nrow(origdata)
  index <- sample(1:n, n, replace = TRUE)
  newdata <- origdata
  newdata[, 2] <- origdata[index, 2] 
  return(newdata)
}


# ----------------------------------------------------------------------
# Structural Change Function

strucchangeGGM <- function(process, splitvar, model = "all"){
  
  k <- NCOL(process)
  n <- NROW(process)
  nobs <- n 
  
  ## scale process
  process <- process/sqrt(nobs)
  meat <- crossprod(process)
  J12 <- root.matrix(chol2inv(chol(meat)))
  
  process <- t(J12 %*% t(process))  
  
  if(model == "correlation"){
    ## select parameters to test
    process <- process[, 11:20]
    k <- NCOL(process)
  }
  
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
# 
# strucchangeGGM(nodevars, splitvar, model = "correlation")

# Check against networktree MOB
#networktree(nodevars, splitvars, method = c("mob"), model = "correlation", verbose = TRUE)

# ------------------------------------------------------------
# Obtain p-value trough bootstrap

bootPval <- function(i, p, n, deltacor = 0, bootiter){ 
  # Generate Original Dataset
  origdata <- data.frame(id = 1:n, 
                          GGMSimulationSplit(p = p, n = n, prob = .2, 
                                            delta_interaction = deltacor))
  nodevars <- as.data.frame(origdata[, 3:(p+2)])
  splitvar <- as.data.frame(origdata[, 2])
  
  fit <- mvnfit(nodevars, estfun = TRUE)
  orig_stat <- strucchangeGGM(process = fit$estfun, splitvar)
  
  # Bootstrap from that original dataset
  out_stat <- numeric(bootiter)
  
  for(i in 1:bootiter) {
    run <- FALSE
    while(run == FALSE){
      bootdata <- bootstrapData(coef = fit$coefficients, obs = n, par = p)
      bootnodevars <- as.data.frame(bootdata[, 2:(p+1)])
      bootsplitvar <- as.data.frame(bootdata[, 1])
      
      bootfit <- mvnfit(bootnodevars, estfun = TRUE)
      
      out <- try(strucchangeGGM(bootfit$estfun, bootsplitvar), silent = TRUE)
      if(!inherits(out, "try-error")) run <- TRUE
    }
    out_stat[i] <- out
  }
  out_stat <- out_stat[order(out_stat)]
  closestValue <- which.min(abs(out_stat - orig_stat))
  
  pval <- (1- closestValue/bootiter)
  
  res <- list(p = p, n = n, cor = deltacor, pval = pval)
  
  return(res)
}

#bootPval(1, p = 5, n = 100, bootiter  = 200)

# ------------------------------------------------------------
# Obtain p-value trough permutation

permutationPval <- function(i, p, n, deltacor = 0, bootiter){ 
  
  estim <- FALSE
  while(estim == FALSE){
    # Generate Original Dataset
    
    origdata <- data.frame(id = 1:n, 
                           GGMSimulationSplit(p = p, n = n, prob = .2, 
                                              delta_interaction = deltacor))
    nodevars <- as.data.frame(origdata[, 3:(p+2)])
    process <- mvnfit(nodevars, estfun = TRUE)$estfun
    
    splitvar <- as.data.frame(origdata[, 2])
    
    orig_stat <- try(strucchangeGGM(process, splitvar), silent = TRUE)
    if(!inherits(orig_stat, "try-error")) estim <- TRUE
  }
  # Resample from that original dataset
  out_stat <- numeric(bootiter)
  
  for(i in 1:bootiter) {
    run <- FALSE
    while(run == FALSE){
      permdata <- permutationData(origdata)
      permsplitvar <- as.data.frame(permdata[, 2])
      
      out <- try(strucchangeGGM(process, splitvar = permsplitvar), silent = TRUE)
      if(!inherits(out, "try-error")) run <- TRUE
    }
    out_stat[i] <- out
  }
  
  out_stat <- out_stat[order(out_stat)]
  closestValue <- which.min(abs(out_stat - orig_stat))
  
  pval <- (1 - closestValue/bootiter)
  
  res <- list(p = p, n = n, cor = deltacor, pval = pval)
  
  return(res)
}

#permutationPval(1, p = 5, n = 250, bootiter  = 100)

#----------------------------------------------------------------------
# Simulations Regenerate the p-value using resampling

n <- c(200, 500, 2000)
p <- c(5, 10, 15)
repiter <- 5000
bootiter <- 1000 
cor <- 0 #c(0, 0.05, 0.1, 0.2)
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
                                   bootiter = bootiter, mc.cores = numCores)
    }
  }
}

out <- as.numeric(unlist(res_GGMPerm))
res <- na.omit(out)
pvalPermutationGGM <- as.data.frame(matrix(res, ncol = 4, byrow = T))
colnames(pvalPermutationGGM) <- c("p", "n", "cor", "pval")
write.csv(pvalPermutationGGM, "pval-Permutation-GGM-cont.csv")

