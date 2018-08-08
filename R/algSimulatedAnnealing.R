library(ggplot2)
library(profvis)
source('./R/toolbox.R')
source('./R/greedySoler.R')
source('./R/simulatedAnnealing.R')
NodesA <- 30
NodesB <- 30
Edges <- 120
mtx <- greedy_solve(NodesA, NodesB, Edges)
print(mtx)
print(nestedness_NODF(mtx))
#print(sum(mtx))
mtx <- sim_anneal_optimisation(mtx)
print(nestedness_NODF(mtx))

