#### GGM Data Simulation ####

GGMSimulationSplit <- function(p, n, prob, delta_interaction){
  
  nedges <- p*(p-1)/2 # number of edges
  
  posdef <- FALSE
  
  while(posdef == FALSE) {
    # Step A:
    # generate random network with x nodes, simulate pcor
    graph <- GeneNet::ggm.simulate.pcor(p, etaA=prob)
    diag(graph) <- 1
    
    if (delta_interaction == 0) {
      
      # Simulate only one big dataset in case of no correlation difference
      pcor.inv <- corpcor::pseudoinverse(graph)
      mean <- rep(0, p)
      
      d <- mvtnorm::rmvnorm(n, mean = mean, sigma = pcor.inv)
      posdef <- TRUE
    } else {
      
      # If in variable correlation setting, modify the partial correlation matrix
      change <- sample(0:1, size = nedges, prob = c(1-prob, prob), replace  = TRUE)
      change[change == 1] <- sample(c(+ delta_interaction, - delta_interaction), size = sum(change), replace = TRUE)
      
      #Interaction Matrix
      graph2 <- matrix(0, p, p)
      graph2[lower.tri(graph2)] <- change
      graph2[upper.tri(graph2)] <- t(graph2)[upper.tri(graph2)] 
      graph2 <- graph + graph2
      
      
      if (!corpcor::is.positive.definite(graph2, tol=1e-8)){
        posdef <- FALSE
      } else {
        posdef <- TRUE
      }
      
      pcor.inv <- corpcor::pseudoinverse(graph)
      pcor.inv2 <- corpcor::pseudoinverse(graph2)
      mean <- rep(0, p)
      
      # Sample Data 
      d <- rbind(
        mvtnorm::rmvnorm(n/2, mean = mean,
                         sigma = pcor.inv),
        mvtnorm::rmvnorm(n/2, mean = mean,
                         sigma = pcor.inv2)
      )
    }
    
  }
  
  # Generate a splitting variable/moderator variable
  z1 <- sort(sample(18:75, size = n, replace = TRUE)) # continuous splitting variable
  #z1 <- c(rep(0, n/2), rep(1, n/2)) # binary splitting variable

  
  
  data <- cbind(as.data.frame(z1) , d)
  
  colnames(data)[2:(1+p)] <- paste0("y", 1:p)
  
  return(data)
}
