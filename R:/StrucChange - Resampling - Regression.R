rm(list = ls())
library("readr")
library("networktree")
library("MASS")
library("strucchange")
library("modelr")
library("tidyverse")
library("parallel")

# ----------------------------------------------------------------------
# Bootstrap Function
bootstrapData <- function(origdata){
  n <- nrow(origdata)
  index <- sample(1:n, n, replace = TRUE)
  newdata <- origdata[index, ]
  return(newdata)
}

permutationData <- function(origdata){
  n <- nrow(origdata)
  index <- sample(1:n, n, replace = TRUE)
  newdata <- origdata
  newdata[, 2] <- origdata[index, 2] 
  return(newdata)
}

# ------------------------------------------------------------
# Obtain p-value trough bootstrap

bootPval <- function(i, p, n, deltacor = 0, bootiter){ 
  # Generate Original Dataset
  origdata <- data.frame(id = 1:n, 
                         GGMSimulationSplit(p = p, n = n, prob = .2, 
                                            delta_interaction = deltacor))
  nodevars <- as.data.frame(origdata[, 3:(p+2)])
  
  # Obtain original test statistic for Structural change test 
  form <- paste("y1", "~", paste(colnames(nodevars[,2:p]), collapse=" + "))
  form <- as.formula(form)
  resorig <- gefp(form, fit = lm, scores = estfun, data = nodevars, order.by = origdata[, 2])
  orig_stat <- sctest(resorig, functional = supLM())$statistic
  
  # Bootstrap from that original dataset
  out_stat <- numeric(bootiter)
  
  for(i in 1:bootiter) {
    run <- FALSE
    while(run == FALSE){
      bootdata <- bootstrapData(origdata)
      bootnodevars <- as.data.frame(bootdata[, 3:(p+2)])
      
      form <- paste("y1", "~", paste(colnames(bootnodevars[,2:p]), collapse=" + "))
      form <- as.formula(form)
      resboot <- gefp(form, fit = lm, scores = estfun, data = bootnodevars, order.by = bootdata[, 2])
      out <- try(sctest(resboot, functional = supLM()), silent = TRUE)
      if(!inherits(out, "try-error")) run <- TRUE
    }
    out_stat[i] <- out$statistic
  }
  out_stat <- out_stat[order(out_stat)]
  closestValue <- which.min(abs(out_stat - orig_stat))
  
  pval <- (1- closestValue/bootiter)
  
  res <- list(p = p, n = n, cor = deltacor, pval = pval)
  
  return(res)
}

#bootPval(1, p = 5, n = 100, bootiter  = 100)

# ------------------------------------------------------------
# Obtain p-value trough permutation

permutationPval <- function(i, p, n, deltacor = 0, bootiter){ 
  
  # Generate Original Dataset
  
  origdata <- data.frame(id = 1:n, 
                         GGMSimulationSplit(p = p, n = n, prob = .2, 
                                            delta_interaction = deltacor))
  nodevars <- as.data.frame(origdata[, 3:(p+2)])
  
  # Obtain original test statistic for Structural change test 
  form <- paste("y1", "~", paste(colnames(nodevars[,2:p]), collapse=" + "))
  form <- as.formula(form)
  resorig <- gefp(form, fit = lm, scores = estfun, data = nodevars, order.by = origdata[, 2])
  orig_stat <- sctest(resorig, functional = supLM())$statistic
  
  # Resample from that original dataset
  out_stat <- numeric(bootiter)
  
  for(i in 1:bootiter) {
    run <- FALSE
    while(run == FALSE){
      permdata <- permutationData(origdata)
      permnodevars <- as.data.frame(permdata[, 3:(p+2)])
      
      form <- paste("y1", "~", paste(colnames(permnodevars[,2:p]), collapse=" + "))
      form <- as.formula(form)
      resperm <- gefp(form, fit = lm, scores = estfun, data = permnodevars, order.by = permdata[, 2])
      out <- try(sctest(resperm, functional = supLM()), silent = TRUE)
      if(!inherits(out, "try-error")) run <- TRUE
    }
    out_stat[i] <- out$statistic
  }
  
  out_stat <- out_stat[order(out_stat)]
  closestValue <- which.min(abs(out_stat - orig_stat))
  
  pval <- (1 - closestValue/bootiter)
  
  res <- list(p = p, n = n, cor = deltacor, pval = pval)
  
  return(res)
}

#permutationPval(1, p = 5, n = 250, bootiter  = 1000)

#----------------------------------------------------------------------
# Simulations Regenerate the p-value 

n <- c(50, 200, 1000)
p <- c(3, 5, 9)
repiter <- 5000
bootiter <- 1000 
cor <- c(0, 0.05, 0.1, 0.2)
numCores <- detectCores()
cntr <- 0
res_RegPerm <- list()

for (pi in 1:length(p)) {
  for (ni in 1:length(n)){
    for (ci in 1:length(cor)) {
      cntr <- cntr + 1 
      pit <- p[pi]
      nit <- n[ni]
      cit <- cor[ci]
      res_RegPerm[[cntr]] <- mclapply(1:repiter, permutationPval, p = pit, n = nit,
                                      deltacor = cit,
                                      bootiter = bootiter, mc.cores = numCores)
    }
  }
}

out <- as.numeric(unlist(res_RegPerm))
res <- na.omit(out)
pvalPermutationReg <- as.data.frame(matrix(res, ncol = 4, byrow = T))
colnames(pvalPermutationReg) <- c("p", "n", "cor", "pval")
write.csv(pvalPermutationReg, "pval-Permutation-Reg-cont.csv")
