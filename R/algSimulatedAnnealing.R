source('./R/toolbox.R')
source('./R/greedySoler.R')
source('./R/simulatedAnnealing.R')
library(profvis)
library(ggplot2)
library(profvis)
NodesA <- 14
NodesB <- 13
Edges <- 52
greedy_solve(NodesA, NodesB, Edges)
# print("Constructing a greedy solution:")
# Rprof(greedy_solve(NodesA, NodesB, Edges), gc.profiling = TRUE, line.profiling = TRUE)
#pd <- readProfileData(filename = "Rprof.out")
#head(funSummary(pd), 10)
# mtx <- sim_anneal_optimisation(mtx)
# print(nestedness_NODF(mtx))

