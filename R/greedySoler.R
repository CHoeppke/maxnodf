source('./R/toolbox.R')

greedy_add_link <- function(mtx, support_data){
    # Find the most promising zeros:
    zPosList <- websearch_NODF_fast(mtx)
    opt_nodf <- -100.0
    opt_pos <- c(-1,-1)
    for(idx in 1:nrow(zPosList)){
        zPos <- zPosList[idx, ]
        my_res <- nodf_one_link_added(mtx, zPos, support_data)
        new_nodf <- my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
        if(new_nodf > opt_nodf){
            opt_nodf <- new_nodf
            opt_pos <- zPos
        }
        # Revert the change:
        my_res <- nodf_one_link_removed(mtx, zPos, support_data)
        new_nodf <- my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
    }
    # actually perform the update
    my_res <- nodf_one_link_added(mtx, opt_pos, support_data)
    nodf <- my_res[[1]]
    mtx <- my_res[[2]]

    return(my_res)
}

greedy_add_link_smart <- function(mtx, support_data){
    # Find the most promising zeros:
    NodesA = nrow(mtx)
    NodesB = ncol(mtx)
    zPosList <- websearch_NODF_fast(mtx)
    opt_nodf <- -100.0
    opt_pos <- c(-1,-1)
    for(idx in 1:nrow(zPosList)){
        zPos <- zPosList[idx, ]
        my_res <- nodf_one_link_added(mtx, zPos, support_data)
        new_nodf <- my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
        if(new_nodf > opt_nodf){
            opt_nodf <- new_nodf
            opt_pos <- zPos
        }
        if(new_nodf == opt_nodf){
            score_old <- (opt_pos[1] / NodesA) + (opt_pos[2] / NodesB)
            score_new <- (zPos[1] / NodesA) + (zPos[2] / NodesB)
            if(score_new < score_old){
                opt_nodf <- new_nodf
                opt_pos <- zPos
            }
        }
        # Revert the change:
        my_res <- nodf_one_link_removed(mtx, zPos, support_data)
        new_nodf <- my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
    }
    # actually perform the update
    my_res <- nodf_one_link_added(mtx, opt_pos, support_data)
    nodf <- my_res[[1]]
    mtx <- my_res[[2]]

    return(my_res)
}

greedy_solve <- function(NodesA, NodesB, Edges){
    mtx <- matrix(0.0, nrow=NodesA, ncol=NodesB)
    mtx[1,] = 1.0
    mtx[,1] = 1.0
    mtx[2,2] = 1.0
    support_data <- init_nodf(mtx)
    nodf <- nestedness_NODF(mtx)
    tp <- txtProgressBar(min = sum(mtx), max = Edges, style = 3)
    while(sum(mtx) < Edges){
        my_res <- greedy_add_link_smart(mtx, support_data)
        nodf <-my_res[[1]]
        mtx <- my_res[[2]]
        support_data <- my_res[[3]]
        setTxtProgressBar(tp, sum(mtx))
    }
    print(sum(mtx))
    return(mtx)
}

