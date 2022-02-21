## code to prepare `DATASET` dataset goes here

library(magrittr)

rm(list=ls(all.names = T));

adj <- read.table("inst/adj.txt", header = T)

set.seed(22)

#---   setting
p1 <- 10  #地域数
p2 <- 20  #年数
n <- p1 * p2
m <- sample(100:1000, n, replace=T)

#---   adjacency
D <- lapply(1:p2, function(j){
  Dj <- split(
    adj$adjNO + p1*(j-1),
    adj$areaNO + p1*(j-1)
  )

  if(j == 1)
  {
    Dj <- lapply(1:p1, function(l){
      c(Dj[[l]], l + p1)
    })
    names(Dj) <- 1:p1
  } else if(j == p2)
  {
    Dj <- lapply(1:p1, function(l){
      c(l + p1*(j-2), Dj[[l]])
    })
    names(Dj) <- (1:p1) + p1*(j-1)
  } else
  {
    Dj <- lapply(1:p1, function(l){
      c(l + p1*(j-2), Dj[[l]], l + p1*j)
    })
    names(Dj) <- (1:p1) + p1*(j-1)
  }

  return(Dj)
}) %>% do.call("c", .)

#---  true adjacency
adjD <- lapply(1:n, function(i){
  cbind(lab=i, adj=D[[i]])
}) %>% do.call(rbind, .)

b.sta <- 80
base <- sample(1:n, b.sta) %>% sort
Tgroup <- numeric(n)
Tgroup[base] <- paste0("g", 1:b.sta)
idx <- which(Tgroup == "0")

while("0" %in% Tgroup)
{
  i <- idx[1]

  cand <- D[[i]][!D[[i]] %in% idx]

  if(length(cand) == 0)
  {
    join <- F
  } else
  {
    join <- rbinom(1, 1, 0.5) == 1
  }

  if(join)
  {
    Tgroup[i] <- sample(Tgroup[cand] %>% unique, 1)
    idx <- idx[-1]
  } else
  {
    idx <- c(idx[-1], idx[1])
  }
} #end while

#---   simulation data
Mu.sta0 <- rnorm(b.sta)
Mu.sta <- numeric(n)

for(j in 1:b.sta)
{
  Mu.sta[Tgroup == paste0("g", j)] <- Mu.sta0[j]
}

Pi <- exp(Mu.sta) / (1 + exp(Mu.sta))
y <- mapply(rbinom, 1, m, Pi)

exdata <- data.frame(
  y = y,
  m = m,
  area = rep(1:p1, p2),
  year = rep(1:p2, each=p1),
  lab = 1:n
)

#---   output
adj <- data.frame(adjD)
usethis::use_data(exdata, overwrite = TRUE)
usethis::use_data(adj, overwrite = TRUE)
