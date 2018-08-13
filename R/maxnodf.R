#' @useDynLib maxnodf
#' @importFrom Rcpp sourceCPP
NULL

maxnodf <- function(NodesA, NodesB, Edges, quality = 0){
    if(quality == 0){
        # Use a basic greedy algorithm to find the optimum
        mtx <- greedySolve2(NodesA, NodesB, Edges)
    }else if(quality == 1){
        # Use a greedy algorithm with one round of hill climbing to
        # improve the results
        mtx <- greedySolve2(NodesA, NodesB, Edges)
        mtx <- full_hill_climb_cpp(NodesA, NodesB, Edges)
    }else if(quality == 2){
        # Use a full round of simulated annealing to find a highly
        # improved solution.
        mtx <- greedySolve2(NodesA, NodesB, Edges)
        mtx <- sim_anneal_opt_cpp(mtx)
    }else{
        print("Please chose a valid quality parameter. Options: ")
        print("quality = 0 -> Use a very fast greedy algorithm.")
        print("quality = 1 -> Improved result using hillclimbing in combination with greedy.")
        print("quality = 2 -> Use a simulated annealing algorith. Best results but requires the most computation time")
    }
    return(mtx)
}
