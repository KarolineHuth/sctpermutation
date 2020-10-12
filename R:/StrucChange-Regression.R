library("networktree")
library("MASS")
library("strucchange")
library("parallel")

# -------------------------------------------

# Function to simulate data and compute structural change test for linear regression model

strucChangeReg <- function(i, p, n, deltacor = 0){ 
  
  run <- FALSE
  while(run == FALSE) {
    data <- GGMSimulationSplit(p, n, prob = 0.2, delta_interaction = deltacor)
    nodevars <- as.data.frame(data[, 2:(p+1)])
    
    # Structural change test 
    form <- paste("y1", "~", paste(colnames(nodevars[,2:p]), collapse=" + "))
    form <- as.formula(form)
    
    res <- lm(form, data = nodevars)
    out_DM <- try(strucchange::sctest(res, functional = "DM", order.by = data[, 1], vcov = "info"), 
                  silent = TRUE)
    if(!inherits(out_DM, "try-error")) run <- TRUE
  }
  
  # Check for all test statistics
  out_CvM <- strucchange::sctest(res, functional = "CvM", order.by = data[, 1], vcov = "info")
  out_maxLM <- strucchange::sctest(res, functional = "maxLM", order.by = data[, 1], vcov = "info")
  
  res <- list(p = p, n = n, cor = deltacor, 
              DMpval = out_DM$p.value, DMstat = out_DM$statistic, 
              CvMpval = out_CvM$p.value, CvMstat = out_CvM$statistic, 
              maxLMpval = out_maxLM$p.value, maxLMstat = out_maxLM$statistic)
  return(res)
}

#strucChangeReg(1, 5, 100)

# -------------------------------------------

# Simulation of SCT under the Null-Hypothesis

n <- c(50, 200, 1000) # sample size
p <- c(3, 5, 9) # model size
cor <- c(0) # invariance violation (set to 0 to simulate data under the null - no measurement invariance violation)
repiter <- 5000 # number of repetitions 
numCores <- detectCores()
res_reg <- list()
cntr <- 0

for (pi in 1:length(p)) {
  for (ni in 1:length(n)) {
    for (ci in 1:length(cor)) {
      cntr <- cntr + 1 
      pit <- p[pi]
      nit <- n[ni]
      cit <- cor[ci]
      res_reg[[cntr]] <- mclapply(1:repiter, strucChangeReg, p = pit, n = nit, 
                                  deltacor = cit, mc.cores = numCores)
    }
  }
}

out <- unlist(res_reg)
res <- na.omit(as.numeric(out))
resReg <- as.data.frame(matrix(res, ncol = 11, byrow = T))
colnames(resReg) <- c("p", "n", "cor", "DMpval", "DMstat",
                      "CvMpval", "CvMstat", 
                      "maxLMpval", "maxLMstat")
resReg %>% group_by(n, p, cor) %>% count()
write.csv(resReg, "pval-StrucChange-Reg-cont.csv")


