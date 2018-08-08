source('./R/toolbox.R')

hill_climb_step <- function(mtx, R){
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    oPosList <- get_valid_ones(mtx)
    support_data <- init_nodf(mtx)
    opt_mtx <- mtx
    opt_nodf <- nestedness_NODF(mtx)
    for(idx in 1:nrow(oPosList)){
        opos <- oPosList[idx,]
        for(xshift in -R:R){
            for(yshift in -R:R){
                newx <- opos[1] + xshift
                newy <- opos[2] + yshift
                if(newx>= 1 & newx <= NodesA & newy >= 1 & newy <= NodesB){
                    if(mtx[newx, newy] == 0){
                        zpos <- c(newx, newy)
                        my_res <- nodf_one_link_removed(mtx, opos, support_data)
                        mtx <- my_res[[2]]
                        support_data <- my_res[[3]]
                        my_res <- nodf_one_link_added(mtx, zpos, support_data)
                        nodf <- my_res[[1]]
                        if(nodf > opt_nodf){
                            opt_mtx <- my_res[[2]]
                            opt_nodf <- nodf
                        }
                        my_res <- nodf_one_link_added(mtx, opos, support_data)
                        mtx <- my_res[[2]]
                        support_data <- my_res[[3]]
                    }
                }
            }
        }
    }
    return(opt_mtx)
}

full_hill_climb <- function(mtx, R=1){
    old_nodf <- -100.0
    count <- 0
    while(old_nodf < nestedness_NODF(mtx)){
        count <- count + 1
        old_nodf <- nestedness_NODF(mtx)
        mtx <- hill_climb_step(mtx, R)
    }
    return(mtx)
}

