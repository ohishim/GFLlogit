
#Coordinate optimization for GFL logistic regression

#arguments:
#  m: vector of the numbers of trials
#  y: vector of the numbers of successes
#  adjCD: list of adjacent information
#  lambda: tuning parameter; "NULL" or vector/scalar
#  lambda.type: option when is.null(lambda)=F
#    value: lambda is searching points
#    rate: lam.max*lambda is searching points
#  weight: penalty weights; "NULL" or list
#  conv: "coef" or "pll"
#    A convergence judgement is based on "conv"
#  progress: If TRUE, progress is displayed

#output:
#  mu.hat: vector of estimates
#  lam.hat: the optimal tuning parameter
#  pll: the minimum of the objective function
#  time: runtime (s)

library(magrittr)

GFLlogit <- function(m, y, adjCD, lambda=NULL, lambda.type=NULL, weight=NULL, conv="coef", progress=FALSE){
  
  ##############################################################################
  ###   preparation
  ##############################################################################
  
  n <- length(y)
  r <- sapply(adjCD, length)
  
  Mu.MLE <- log(y / (m - y))
  mu.max <- log(sum(y) / sum(m-y))
  
  if(is.null(weight))
  {
    Wli <- lapply(1:n, function(i){
      mui <- Mu.MLE[i]
      Mui <- Mu.MLE[adjCD[[i]]]
      
      return(1 / abs(mui - Mui))
    })
  } else
  {
    Wli <- weight
  }
  
  if(is.null(lambda))
  {
    Lam.max0 <- sapply(1:n, function(i){ 
      max(
        (y[i] - (m[i]-y[i])*exp(mu.max)) / (2*(1+exp(mu.max))*sum(Wli[[i]])),
        ((m[i]-y[i])*exp(mu.max) - y[i]) / (2*(1+exp(mu.max))*sum(Wli[[i]]))
      )
    })
    lam.max <- Lam.max0 %>% max
    
    Lambda <- rev(lam.max * ((3/4)^((1:100)-1)))
  } else if(is.null(lambda.type))
  {
    stop("Please select lambda.type")
  } else if(lambda.type == "rate")
  {
    Lam.max0 <- sapply(1:n, function(i){
      max(
        (y[i] - (m[i]-y[i])*exp(mu.max)) / (2*(1+exp(mu.max))*sum(Wli[[i]])),
        ((m[i]-y[i])*exp(mu.max) - y[i]) / (2*(1+exp(mu.max))*sum(Wli[[i]]))
      )
    })
    lam.max <- Lam.max0 %>% max
    
    Lambda <- lam.max * lambda
  } else if(lambda.type == "value")
  {
    Lambda <- lambda
  } else
  {
    stop("The lambda.type is missing")
  }
  
  lam.n <- length(Lambda)
  
  PLL <- function(Mu, lambda, y, m, adjCD, n, Wli){
    L <- -sum(y*Mu) + sum(m*log(1+exp(Mu)))
    P <- lapply(1:n, function(i){
      wi <- Wli[[i]]
      mui <- Mu[i]
      Mu.adj <- Mu[adjCD[[i]]]
      return(sum(wi * abs(mui - Mu.adj)))
    }) %>% unlist %>% sum
    
    return(L + lambda*P)
  }
  
  ##############################################################################
  ###   main
  ##############################################################################
  
  Mu.aft <- Mu.MLE
  
  if(lam.n > 1){MU <- matrix(0, n, lam.n)}
  
  t.s <- proc.time()[3]
  for(lam.i in 1:lam.n)
  {
    lambda <- Lambda[lam.i]
    if(conv == "pll"){pll.aft <- PLL(Mu.aft, lambda, y, m, adjCD, n, Wli)}
    
    dif <- 1
    iter <- 0
    
    while(dif > thres)
    {
      iter <- iter + 1
      Mu.bef <- Mu.aft
      if(conv == "pll"){pll.bef <- pll.aft}
      
      ##########################################################################
      ###   descent cycle
      ##########################################################################
      
      for(i in 1:n)
      {
        Di <- adjCD[[i]]
        Mu.adj <- Mu.aft[Di]  
        wi <- Wli[[i]]
        ri <- r[i]
        
        if(length(unique(Mu.adj)) < ri)
        {
          
          Mu.adj0 <- Mu.adj
          Mu.adj <- unique(Mu.adj0)
          labi <- match(Mu.adj0, Mu.adj)
          wi <- split(wi, labi) %>% sapply(., sum)
        }
        
        Mu.aft[i] <- GFLlogit1(m[i], y[i], lambda, wi, Mu.adj)
        
      } #end for i
      
      ##########################################################################
      ###   fusion cycle
      ##########################################################################
      
      Xi <- unique(Mu.aft)
      gr.labs <- match(Mu.aft, Xi)
      b <- unique(gr.labs) %>% length
      
      if(b < n)
      {
        if(b == 1)
        {
          Mu.aft <- rep(mu.max, n)
        } else
        {
          E <- split(1:n, gr.labs)
          idx.f <- (sapply(E, length) > 1) %>% which
          
          for(l in idx.f)
          {
            El <- E[[l]]
            
            Fl <- adjCD[El] %>% unlist %>% unique %>% sort %>% setdiff(., El)
            
            if(length(Fl) > 0)
            {
              Dl <- sapply((1:b)[-l], function(s){
                if(length(intersect(E[[s]], Fl)) > 0){return(s)}
              }) %>% unlist
              
              Jl <- lapply(Dl, function(i){
                Ei <- E[[i]]
                
                Out <- lapply(El, function(j){
                  Dj <- adjCD[j][[1]]
                  if(length(intersect(Ei, Dj)) > 0)
                  {
                    return(cbind(j, intersect(Ei, Dj)))
                  }
                }) %>% do.call(rbind, .)
              })
              
              wl <- lapply(Jl, function(x){
                lapply(1:nrow(x), function(i){
                  Wli[[x[i,1]]][which(x[i,2] == adjCD[x[i,1]][[1]])]
                }) %>% unlist %>% sum
              }) %>% unlist
              
              Xi.adj <- Xi[Dl]
              rl <- length(Dl)
              
              if(length(unique(Xi.adj)) < rl)
              {
                Xi.adj0 <- Xi.adj
                Xi.adj <- unique(Xi.adj0)
                labl <- match(Xi.adj0, Xi.adj)
                wl <- split(wl, labl) %>% sapply(., sum)
              }
              
              Mu.aft[El] <- Xi[l] <- GFLlogit1(sum(m[El]), sum(y[El]), lambda, wl, Xi.adj)
              
            } else 
            {
              Mu.aft[El] <- Xi[l] <- log(sum(y[El]) / (sum(m[El] - y[El])))
            } #end if
          } #end for l
        } #end if
      } #end if
      
      ##############################################################################
      ###   convergence judgement
      ##############################################################################
      
      if(conv == "pll")
      {
        if(b == 1)
        {
          dif <- 0
        } else
        {
          pll.aft <- PLL(Mu.aft, lambda, y, m, adjCD, n, Wli)
          dif <- abs(pll.aft - pll.bef)
        }
      } else if(conv == "coef")
      {
        if(b == 1)
        {
          dif <- 0
        } else
        {
          dif <- max((Mu.aft - Mu.bef)^2) / max(Mu.bef^2)
        }
      }
      
      if(progress){print(paste(lam.i, iter, dif, sep="_"))}
      
    } #end while dif
    
    if(lam.n > 1){MU[,lam.i] <- Mu.aft}
    
  } #end for lam.i
  t.f <- proc.time()[3]
  
  ##############################################################################
  ###   select tuning parameter
  ##############################################################################
  
  if(lam.n == 1){
    Mu.hat <- Mu.aft
    lam.hat <- lambda
  } else
  {
    LL <- function(Mu){-sum(y*Mu) + sum(m*log(1+exp(Mu)))}
    BIC <- function(Mu){2*LL(Mu) + log(n)*length(unique(Mu))}
    
    opt <- apply(MU, 2, BIC) %>% which.min
    Mu.hat <- MU[,opt]
    lam.hat <- Lambda[opt]
  }
  
  pll.hat <- PLL(Mu.hat, lambda, y, m, adjCD, n, Wli)
  
  ##############################################################################
  ###   output
  ##############################################################################
  
  return(list(
    mu.hat = Mu.hat,
    lam.hat = lam.hat,
    pll = pll.hat,
    time = t.f-t.s
  ))
}