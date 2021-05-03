#rm(list = ls())
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

# ------------------------------------------------------------
# Obtain p-value trough permutation

permutationPval <- function(i, p, n, deltacor = 0, permiter){ 
  
  estim <- FALSE
  while(estim == FALSE){
    # Generate Original Dataset
    
    origdata <- data.frame(id = 1:n, 
                           GGMSimulationSplit(p = p, n = n, prob = .5, 
                                              delta_interaction = deltacor))
    nodevars <- as.data.frame(origdata[, 3:(p+2)])
    #process <- networktree::mvnfit(nodevars, estfun = TRUE)$estfun
    
    splitvar <- origdata[, 2]
    splitvars <- as.data.frame(splitvar)
    #MOB direct
    minsplit <- p*(p-1)/2 + 1
    if(minsplit > n/2) {minsplit <- (n/2 -10)}
    ## MOB
    out <-  try(networktree::networktree(nodevars = nodevars, splitvars = splitvars, method = c("mob"), model = "correlation",
                                        verbose = F, minsplit = minsplit, maxdepth = 2))

    ## own Function 
    #out <- try(strucchangeGGM(process, splitvar, minsplit), silent = TRUE)
    if(!inherits(out, "try-error")) {
      estim <- TRUE
      ## MOB
      orig_stat <- unlist(out)$node.info.test1
      ## own function
      #orig_stat <- out
    }
  }
  # Resample from that original dataset
  out_stat <- numeric(permiter)
  
  for(f in 1:permiter) {
    run <- FALSE
    while(run == FALSE){
      permdata <- permutationData(origdata)
      permsplitvar <- permdata[, 2]
      permsplitvars <- as.data.frame(permsplitvar)
      
      out <-  try(networktree::networktree(nodevars = nodevars, splitvars = permsplitvars, method = c("mob"), model = "correlation",
                                           verbose = F, minsplit = minsplit, maxdepth = 2))
      #out <- try(strucchangeGGM(process, splitvar = permsplitvar, minsplit), silent = TRUE)
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
  
  res <- list(p = p, n = n, cor = deltacor, pval = pval, testStat = orig_stat)
  
  return(res)
}

# permutationPval(1, p = 10, n = 200, deltacor = 0.1,  permiter  = 1000)

#----------------------------------------------------------------------
# Simulations Regenerate the p-value using resampling

n <- c(200, 500, 2000) # sample size
p <- c(5, 10, 15) # size of network
repiter <- 1000 # number of repetitions 
permiter <- 5000 # number of samples for permutation approach
cor <- c(0.1, 0.3, 0.5) # measurement invariance violation (set to 0 to simulate data under the null - no measurement invariance violation)
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
      
      ## Printing intermediate results
      out <- as.numeric(unlist(res_GGMPerm))
      res <- na.omit(out)
      pvalPermutationGGM <- as.data.frame(matrix(res, ncol = 5, byrow = T))
      colnames(pvalPermutationGGM) <- c("p", "n", "cor", "pval", "testStat")
      intermediate_res <- pvalPermutationGGM %>%
        filter(cor != 0) %>%
        mutate(sig = ifelse(pval < 0.05, 1L, 0L)) %>% 
        group_by(n, p, cor) %>% 
        na.omit() %>%
        dplyr::summarize(count_all = n(),
                         power = mean(sig, na.rm = TRUE))
      print(intermediate_res)
    }
    print(nit)
  }
  print(pit)
}

out <- as.numeric(unlist(res_GGMPerm))
res <- na.omit(out)
pvalPermutationGGM <- as.data.frame(matrix(res, ncol = 5, byrow = T))
colnames(pvalPermutationGGM) <- c("p", "n", "cor", "pval", "testStat")
write.csv(pvalPermutationGGM, "power-Permutation-GGM-cont.csv")
head(pvalPermutationGGM)
