non_centrality <- function(gamma, pa, N){
  p_plus <- (gamma * pa)/((gamma - 1) * pa + 1)
  p_minus <- pa
  p_mean <- (p_plus + p_minus) / 2
  lambda <- (p_plus - p_minus) / (sqrt(2 * p_mean *(1 - p_mean)))
  non_centrality_parameter = lambda * sqrt(N)
  return(non_centrality_parameter)
}

power <- function(gamma, pa, N, alpha=0.05){
  return(pnorm(qnorm(alpha/2) + non_centrality(gamma, pa, N)) + 1 -
         pnorm((-qnorm(alpha/2) + non_centrality(gamma, pa , N))))
}

number_individuals <- function(gamma, pa){
  individuals <- seq(1:2000)
  ind <- numeric()
  close_to_80 <- numeric() 
  for (i in individuals)
  {
    if (abs(power(gamma,pa,i) - .80) <= .005) {
      close_to_80 <- c(close_to_80, abs(power(gamma,pa,i) - .80))
      ind <- c(ind, i)
    }
  }
  return (ind[which.min(close_to_80)])
}

gamma <- c(1.5, 2.0, 3.0)
pa <- c(0.05, 0.2, 0.4)
N <- c(500, 1000)

## Non Centrality Parameter
for (n in N) {
  for (g in gamma) {
    for (p in pa){
      cat("Gamma:", g, " PA:", p, " N:", n, " NCP = ", non_centrality(g,p,n), "\n")
      cat("Gamma:", g, " PA:", p, " N:", n, " POWER = ", power(g,p,n), "\n")
      x <- number_individuals(g, p)
      cat("Gamma:", g, " PA:", p, " N:", x, " Power = ", power(g,p,x), "\n")  
      cat("\n")
    }
  }
}
