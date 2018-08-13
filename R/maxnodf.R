#' Compute NODF-maximising graph
#'
#' @param quality   An optional quality parameter to control the
#'                  tradeoff between computation time and result quality.
#' @return Approximation of the NODF maximising graph.
#' @examples
#' maxnodf(14, 13, 52)
#' maxnodf(14, 13, 52, 2)
maxnodf <- function(web, quality = 0){
    #' @useDynLib maxnodf
    #' @import Rcpp
    #' @export
    web[web>0] <- 1
    if(any(web < 0)){
        stop("Invalid network. Ensure all elements of web >= 0.")
    }
    if( !quality %in% 0:2){
        stop("Please chose a valid quality parameter. Options: \n\tquality = 0 -> Use a very fast greedy algorithm.\n\tquality = 1 -> Improved result using hillclimbing in combination with greedy.\n\tquality = 2 -> Use a simulated annealing algorith. Best results but requires the most computation time")
    }
    mt_0  <- computeMT0(web)
    mt_t  <- computeMT0(web)
    if(any(mt_0 < 1)){
        stop("Invalid network. Ensure all marginal totals are >= 1.")
    }
    if(any(mt_t < 1)){
        stop("Invalid network. Ensure all marginal totals are >= 1.")
    }
    NodesA <- nrow(web)
    NodesB <- ncol(web)
    Edges <- sum(web)
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
    }
    return(list(nodf_cpp(mtx), mtx))
}

