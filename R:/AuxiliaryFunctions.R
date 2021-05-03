library(corpcor)
GGMSimulationSplit <- function(p, n, prob, delta_interaction){
  
  nedges <- p*(p-1)/2 # number of edges
  
  posdef <- FALSE
  while(posdef == FALSE) {
    # generate random network with x nodes, simulate pcor
    graph <- GeneNet::ggm.simulate.pcor(p, etaA = prob)
    sigma <- corpcor::pcor2cor(graph)
    
    if (delta_interaction == 0) {
      
      # Simulate only one big dataset in case of no correlation difference
      mean <- rep(0, p)
      
      d <- as.data.frame(mvtnorm::rmvnorm(n, mean = mean, sigma = sigma))
      posdef <- TRUE
    } else {
      
      sigma2 <- sigma
      # to avoid correlation above |1| 
      modVals <- c(sigma2[1, 2] + delta_interaction, sigma2[1, 2] - delta_interaction)
      modVal <- modVals[which.min(abs(modVals))]
      sigma2[1, 2] <- sigma2[2, 1] <- modVal
      
      
      #if(all(eigen(graph)$values > 0) & all(eigen(graph2)$values > 0))
      # if (!corpcor::is.positive.definite(sigma2, tol=1e-8)){
      #   posdef <- FALSE
      # } else {
        mean <- rep(0, p)
        sigma2 <- as.matrix(Matrix::nearPD(sigma2)$mat)
        # Sample Data 
        data1 <- try(mvtnorm::rmvnorm(n/2, mean = mean,
                                      sigma = sigma))
        data2 <- try( mvtnorm::rmvnorm(n/2, mean = mean,
                                       sigma = sigma2))
        
        
        if(!inherits(data1, "try-error") & !inherits(data2, "try-error")) {
          posdef <- TRUE
          d <- as.data.frame(rbind(data1, data2))
        }
      # }
    }
  }
  
  # Generate a splitting variable/moderator variable
  #z1 <- sort(sample(18:75, size = n, replace = TRUE)) # continuous splitting variable
  z1 <- sort(rnorm(n)) 
  #z1 <- c(rep(0, n/2), rep(1, n/2)) # binary splitting variable
  
  data <- cbind(as.data.frame(z1) , d)
  
  colnames(data)[2:(1+p)] <- paste0("y", 1:p)
  
  return(as.data.frame(data))
  #return(list(data, graph, graph2))
}
#data <- GGMSimulationSplit(p = 5, n = 500, prob = 0.5, delta_interaction = .5)

