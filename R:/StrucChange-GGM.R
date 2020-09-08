rm(list = ls())

library("networktree")
library("MASS")
library("strucchange")
library("parallel")
library("tidyverse")

### TO DO'S:
# change reading in txt-file to reading in using string 

# -------------------------------------------

strucChangeGGM <- function(i, p, n, deltacor = 0){ 
  n.nodes <- p
  n.obs <- n
  minsplit <- p*(p-1)/2
  if(minsplit > n/2) {minsplit <- (n/2 -10)} # change minsplit size in case it is smaller than the number of parameters
  run <- FALSE
  while(run == FALSE) {
    data <- GGMSimulationSplit(p = n.nodes, n = n.obs, prob = 0.2, 
                                             delta_interaction = deltacor)
    nodevars <- as.data.frame(data[, 2:(p+1)])
    splitvar <- data[, 1]
    splitvars <- as.data.frame(splitvar)
    
    out <- try(networktree::networktree(nodevars, splitvars, method = c("mob"), model = "correlation", 
                    verbose = TRUE, minsplit = minsplit))
    if(!inherits(out, "try-error")) run <- TRUE
  }
  
  # Structural change test 
  capture.output(networktree::networktree(nodevars, splitvars, method = c("mob"), model = "correlation", 
                             verbose = TRUE, minsplit = minsplit), file = paste0("out", i, ".txt"))
  out <- readr::read_delim(paste0("out", i, ".txt"), delim = "-")
  
  # Very complicated but necessary way to get data out
  out <- strsplit(as.character(out$X1), ' ')
  outunl <- unlist(out)
  p.in <- which(grepl("p.value", outunl))
  s.in <- which(grepl("statistic", outunl))
  best <- which(grepl("Best", outunl))
  #p-value
  p.val <- outunl[p.in:(best-1)]
  p.val <- p.val %>% as.numeric() 
  indexp <- which(!is.na(p.val))
  pval <- p.val[indexp]
  #statistic
  statistic <- outunl[s.in:(p.in-1)]
  statistic <- statistic %>% as.numeric() 
  indexs <- which(!is.na(statistic))
  stat <- statistic[indexs]
  
  
  res <- list(p = p, n = n, cor = deltacor, 
              pval = pval, stat = stat)
  return(res)
}

#strucChangeGGM(1, p = 15, deltacor = 0, n = 200)

n <- c(200, 500, 2000)
p <- c(5, 10, 15)
cor <- c(0, 0.05, 0.1, 0.2)
repiter <- 5000
numCores <- detectCores()
res_GGM <- list()
cntr <- 0

for (pi in 1:length(p)) {
  for (ni in 1:length(n)){
    for (ci in 1:length(cor)) {
      cntr <- cntr + 1 
      pit <- p[pi]
      nit <- n[ni]
      cit <- cor[ci]
      res_GGM[[cntr]] <- mclapply(1:repiter, strucChangeGGM, p = pit, n = nit, 
                                  deltacor = cit, mc.cores = numCores)
    }
  }
}

out <- unlist(res_GGM)
res <- na.omit(as.numeric(out))
resGGM <- as.data.frame(matrix(res, ncol = 5, byrow = T))
colnames(resGGM) <- c("p", "n", "cor", "pval", "stat")
trest <- resGGM %>% group_by(n, p, cor) %>% count()
write.csv(resGGM, "pval-StrucChange-GGM-cont.csv")
