source('./R/toolbox.R')
source('./R/greedySoler.R')
source('./R/simulatedAnnealing.R')

maxnodf <- function(NodesA, NodesB, Edges, quality = 0){
    valid_call <- 0
    if(quality == 0){
        # Use a basic greedy algorithm to find the optimum
        valid_call <- 1
    }else if(quality == 1){
        # Use a greedy algorithm with one round of hill climbing to
        # improve the results
        valid_call <- 1
    }else if(quality == 2){
        # Use a full round of simulated annealing to find a highly
        # improved solution.
        valid_call <- 1
    }else{
        print("Please chose a valid quality parameter. Options: ")
        print("quality = 0 -> Use a very fast greedy algorithm.")
        print("quality = 1 -> Use a slower greedy algorithm with better results.")
        print("quality = 2 -> Use a simulated annealing algorith. Slowest but best resutls.")
    }
    if(valid_call == 1){
        command <- "python3 ./py3_src/maxnodf.py"
        command <- paste(command, NodesA,NodesB,Edges,quality, sep = " ")
        system(command)
    }
    # Read output matrix from csv file:
    temp <- read.csv("./outputmtx.csv", header = FALSE)
    mtx <- as.matrix(temp)
    return(mtx)

}
mtx <- maxnodf(14, 13, 52, 2)
print(mtx)
print(nestedness_NODF(mtx))
