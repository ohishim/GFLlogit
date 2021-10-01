
#This is required to execute "GFLlogit.R"

#argument: 
#  m: natural value
#  y: natural value (y < m)
#  lambda: positive value
#  w: vector of positive values
#  z: vector which are all distinct values

#output:
#  x.hat: the minimizer

library(magrittr)

GFLlogit1 <- function(m, y, lambda, w, z){
  
  r <- length(z)
  
  judge0 <- sapply(1:r, function(j){
    zj <- z[j]
    
    a1 <- (m*exp(zj) / (1+exp(zj))) - y + 2*lambda*sum(w*sign(zj - z))
    a2 <- 2*lambda*w[j]
    
    return(c(
      a1 - a2, a1 + a2
    ))
  })
  
  judge1 <- apply(judge0, 2, function(x){
    c(x[1] <= 0, 0 <= x[2]) %>% sum
  })
  
  s <- which(judge1 == 2)
  
  if(length(s) == 1)
  {
    x.hat <- z[s]
  } else if(length(s) == 0)
  {
    judge2 <- sign(judge0[, order(z), drop=F]) %>% apply(2, sum)
    t <- sort(z)
    
    if(all(judge2 == 2))
    {
      wa <- -sum(w)
    } else if(all(judge2 == -2))
    {
      wa <- sum(w)
    } else
    {
      rl <- which(judge2 == -2) %>% max %>% t[.]
      
      idx <- which(z <= rl)
      wa <- sum(w[idx]) - sum(w[-idx])
    }
    
    c <- 2*lambda*wa - y
    x.hat <- log(-c / (m + c))
    
  } else
  {
    f <- function(X){
      sapply(X, function(x){
        m*log(1+exp(x)) - y*x + 2*lambda*sum(w*abs(x-z))
      })
    }
    
    x0 <- z[s]
    x.hat <- f(x0) %>% which.min %>% x0[.]
  }
  
  return(x.hat)
}

