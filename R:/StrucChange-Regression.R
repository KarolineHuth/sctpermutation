library("networktree")
library("MASS")
library("strucchange")
library("parallel")

# -------------------------------------------

strucChangeReg <- function(i, p, n, deltacor = 0){ 
  
  run <- FALSE
  while(run == FALSE) {
    data <- GGMSimulationSplit(p, n, prob = 0.2, delta_interaction = deltacor)
    nodevars <- as.data.frame(data[, 2:(p+1)])
    
    # Structural change test 
    form <- paste("y1", "~", paste(colnames(nodevars[,2:p]), collapse=" + "))
    form <- as.formula(form)
    
    #res <- gefp(form, fit = lm, scores = estfun, data = nodevars, order.by = as.factor(data[, 1]))
    res <- lm(form, data = nodevars)
    ## FOR ORDINAL AND BINARY SV
    out_DM <- try(sctest(res, functional = "DM", order.by = data[, 1], vcov = "info"), 
                  silent = TRUE)
    if(!inherits(out_DM, "try-error")) run <- TRUE
  }
  # out_DM <- sctest(res, functional = "DM", order.by = as.factor(data[, 1])
  # out_CvM <- sctest(res, functional = "CvM", order.by = as.factor(data[, 1]), vcov = "info")
  # out_maxLM <- sctest(res, functional = "maxLM", order.by = as.factor(data[, 1]), vcov = "info")
  # out_LMuo <- sctest(res, functional = "LMuo", order.by = as.factor(data[, 1]), vcov = "info")
  ## FOR CONTINUOUS SV
  #out_DM <- sctest(res, functional = "DM", order.by = data[, 1], vcov = "info")
  out_CvM <- sctest(res, functional = "CvM", order.by = data[, 1], vcov = "info")
  out_maxLM <- sctest(res, functional = "maxLM", order.by = data[, 1], vcov = "info")
  out_LMuo <- sctest(res, functional = "LMuo", order.by = data[, 1], vcov = "info")
  
  res <- list(p = p, n = n, cor = deltacor, 
              DMpval = out_DM$p.value, DMstat = out_DM$statistic, 
              CvMpval = out_CvM$p.value, CvMstat = out_CvM$statistic, 
              maxLMpval = out_maxLM$p.value, maxLMstat = out_maxLM$statistic, 
              LMuopval = out_LMuo$p.value, LMuostat = out_LMuo$statistic)
  return(res)
}

#strucChangeReg(1, 5, 100)

n <- c(50, 200, 1000) #c(2000, 4000) 
p <- c(3, 5, 9)
cor <- c(0, 0.05, 0.1, 0.2)
repiter <- 5000
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
                      "maxLMpval", "maxLMstat", 
                      "LMuopval", "LMuostat"
)
head(resReg)
resReg %>% group_by(n, p, cor) %>% count()
write.csv(resReg, "pval-StrucChange-Reg-cont.csv")


