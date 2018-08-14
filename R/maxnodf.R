#' Compute NODF-maximising network
#'
#' @param quality   An optional quality parameter to control the
#'                  tradeoff between computation time and result quality.
#' @return Nestedness and graph-matrix of the NODF maximising network.
#' @examples
#' maxnodf(web)
#' maxnodf(web, 2)
maxnodf <- function(web, quality = 0){
    #' @useDynLib maxnodf
    #' @import Rcpp
    #' @export
    config <- sanity_check(web, quality)
    NodesA <- config[[1]]
    NodesB <- config[[2]]
    Edges <- config[[3]]
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
    print(" ")
    return(list(nodf_cpp(mtx), mtx))
}

sanity_check <- function(web, quality){
    NodesA <- -1
    NodesB <- -1
    Edges <- -1
    if(is.matrix(web)){
        if(all(is.numeric(web))){
            web[web>0] <- 1
            if(any(web < 0)){
                stop("Invalid network. Ensure all elements of web >= 0.")
            }
            mt_0  <- computeMT0(web)
            mt_t  <- computeMTt(web)
            if(any(mt_0 < 1)){
                stop("Invalid network. Ensure all marginal totals are >= 1.")
            }
            if(any(mt_t < 1)){
                stop("Invalid network. Ensure all marginal totals are >= 1.")
            }
            NodesA <- nrow(web)
            NodesB <- ncol(web)
            Edges <- sum(web)
        }
        else{
            stop("Parameter 'web' is expected to be a numeric matrix or a numeric vector.")
        }
    }else if(is.vector(web)){
        if(length(web) == 3){
            NodesA <- web[[1]]
            NodesB <- web[[2]]
            Edges <- web[[3]]
        }else{
            stop("The vector 'web' is expected to have three entries.")
        }
    }else{
        stop("Parameter 'web' is expected to either be a matrix or a vector containing the matrix dimentions and number of links.")
    }
    if(Edges <= NodesA + NodesB){
        stop("Number of links needs to satisfy 'Links > nrow(web) + ncol(web).")
    }
    if(Edges > NodesA * NodesB){
        stop("Number of links needs to satisfy 'Links <= nrow(web) * ncol(web).")
    }
    if( !quality %in% 0:2){
        stop("Please chose a valid quality parameter. Options: \n\tquality = 0 -> Use a very fast greedy algorithm.\n\tquality = 1 -> Improved result using hillclimbing in combination with greedy.\n\tquality = 2 -> Use a simulated annealing algorith. Best results but requires the most computation time")
    }
    return(list(NodesA, NodesB, Edges, quality))
}
