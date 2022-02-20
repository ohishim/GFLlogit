
exdata <- read.table("data/exdata.txt", header = T)
adj <- read.table("data/adj.txt", header = T)

# devtools::use_data(exdata)
# devtools::use_data(adj)

save(exdata, file="data/exdata.rda")
save(adj, file="data/adj.rda")