#' Compute NODF-maximising graph
#'
#' @param NodesA    The number of rows in the desired network.
#' @param NodesB    The number of columns in the desired network.
#' @param Edges     The number of Edges in the desired network.
#' @param quality   An optional quality parameter to control the
#'                  tradeoff between computation time and result quality.
#' @return Approximation of the NODF maximising graph.
#' @examples
#' maxnodf(14, 13, 52)
#' maxnodf(14, 13, 52, 2)
maxnodf <- function(NodesA, NodesB, Edges, quality = 0){
    #' @useDynLib maxnodf
    #' @import Rcpp
    #' @export
    if(quality == 0){
        # Use a basic greedy algorithm to find the optimum
        mtx <- greedy_solve2(NodesA, NodesB, Edges)
    }else if(quality == 1){
        # Use a greedy algorithm with one round of hill climbing to
        # improve the results
        mtx <- greedy_solve2(NodesA, NodesB, Edges)
        mtx <- full_hill_climb_cpp(mtx)
    }else if(quality == 2){
        # Use a full round of simulated annealing to find a highly
        # improved solution.
        mtx <- greedy_solve2(NodesA, NodesB, Edges)
        mtx <- sim_anneal_opt_cpp(mtx)
    }else{
        print("Please chose a valid quality parameter. Options: ")
        print("quality = 0 -> Use a very fast greedy algorithm.")
        print("quality = 1 -> Improved result using hillclimbing in combination with greedy.")
        print("quality = 2 -> Use a simulated annealing algorith. Best results but requires the most computation time")
    }
    return(mtx)
}

