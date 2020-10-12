library("readr")
library("networktree")
library("MASS")
library("strucchange")
library("modelr")
library("tidyverse")
library("parallel")

### TO DO'S:
# adapt to also allow for other test statistics but maxLM, e.g., CvM and DM

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

permutationPval <- function(i, p, n, deltacor = 0, bootiter){ 
  
  # Generate Original Dataset
  
  origdata <- data.frame(id = 1:n, 
                         GGMSimulationSplit(p = p, n = n, prob = .2, 
                                            delta_interaction = deltacor))
  nodevars <- as.data.frame(origdata[, 3:(p+2)])
  
  # Obtain original test statistic for Structural change test 
  form <- paste("y1", "~", paste(colnames(nodevars[,2:p]), collapse=" + "))
  form <- as.formula(form)
  resorig <- lm(form, data = nodevars)
  orig_stat <- strucchange::sctest(resorig, functional = "maxLM", order.by = origdata[, 2], 
                                   vcov = "info")$statistic
  
  # Resample from that original dataset
  out_stat <- numeric(bootiter)
  
  for(i in 1:bootiter) {
    run <- FALSE
    while(run == FALSE){
      permdata <- permutationData(origdata)
      permnodevars <- as.data.frame(permdata[, 3:(p+2)])
      
      resperm <- lm(form, data = permnodevars)
      out <- try(strucchange::sctest(resperm, functional = "maxLM", order.by = permdata[, 2], vcov = "info"), silent = TRUE)
      if(!inherits(out, "try-error")) run <- TRUE
    }
    out_stat[i] <- out$statistic
  }
  
  # Determine p-value
  pval <- mean(out_stat > orig_stat)
  
  res <- list(p = p, n = n, cor = deltacor, pval = pval)
  
  return(res)
}

#permutationPval(1, p = 5, n = 250, bootiter  = 1000)

#----------------------------------------------------------------------
# Simulations Regenerate the p-value 

n <- c(50, 200, 1000) # sample size
p <- c(3, 5, 9) # model size
repiter <- 10#5000 # number of repetitions 
bootiter <- 10#1000 # number of samples for permutation approach
cor <- c(0) # measurement invariance violation (set to 0 to simulate data under the null - no measurement invariance violation)
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

