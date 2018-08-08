library(nnet)

#function that loads the plant-pollinator network
#input: number = label of a network
#output: web = mutualsitic network
load_data <- function(){
    name <- paste('network.csv',sep='')
    d <- read.csv(file=name,header=FALSE)
    web <- as.matrix(d)
    web[web > 0] = 1
    return(web)
}

#computes the raw NODF
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
    web[web > 0] = 1
    SA <- nrow(web)
    SP <- ncol(web)
    N <- t(web) %*% web # Standard matrix multiplication
    num <- N
    num[lower.tri(num,diag=TRUE)]=1
    den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
    dele <- den - t(den)
    dele[lower.tri(dele,diag=TRUE)] <- 1
    num[dele == 0] <- 0
    den <- pmin(den,t(den))
    den[lower.tri(den,diag=TRUE)] = 1
    nes <- num/den
    nes[lower.tri(nes,diag=TRUE)] = 0
    nes[is.na(nes)] <- 0
    n1 <- sum(nes)
    N <- web %*% t(web)
    num <- N
    num[lower.tri(num,diag=TRUE)]=1
    den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
    dele <- den - t(den)
    dele[lower.tri(dele,diag=TRUE)] <- 1
    num[dele ==0 ] <- 0
    den <- pmin(den,t(den))
    den[lower.tri(den,diag=TRUE)]=1
    nes <- num/den
    nes[lower.tri(nes,diag=TRUE)] = 0
    nes[is.na(nes)] <- 0
    n2 <- sum(nes)
    out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
    return(out)
}

#finds the maximizum raw value of NODF with given connectance and community size
#input: web = mutualistic network
#output: the maximun row value of NODF
max_nest2 <- function(NodesA, NodesB, Edges){
    #binarize the interaction matrix
    SA <- NodesA
    SP <- NodesB
    SI <- Edges
    #initialize the interaction matrix with minimum requirements
    web_opt <- matrix(0, nrow=SA, ncol=SP)
    web_opt[1,] <- 1
    web_opt[,1] <- 1
    web_opt[2,2] <- 1
    #counting the number of
    SI_left <- SI-SP-SA
    if(SI_left>0){
        #search the best possible location
        for(j in 1:SI_left){
            print("\n---------------------\n")
            myedges <- SI_left - j
            #compare all possible locations and the maximum one
            # Print the potential next position to see what exactly
            # is happening here:
            position_potential <- websearch_NODF(web_opt)
            nest_poten <- c()
            for(i in 1:nrow(position_potential)) {
                web_poten <- web_opt
                web_poten[position_potential[i,1],position_potential[i,2]] <- 1
                nest_poten[i] <- nestedness_NODF(web_poten)
                print(c(position_potential[i,1], position_potential[i,2], nest_poten[i]))
            }
            position_the <- which.is.max(nest_poten)
            print(c("--->", position_potential[position_the,1], position_potential[position_the,2]))
            web_opt[position_potential[position_the,1],position_potential[position_the,2]] <- 1
        }
        return (web_opt)
    }
    #this is to prevent the trivial case
    else{
        return(-1)
    }
}

# finds the maximizum raw value of NODF with given
# connectance and community size
# input: web = mutualistic network
# output: the maximun row value of NODF
max_nest <- function(web){
    #binarize the interaction matrix
    web_binary <- web
    web_binary[web_binary > 0] = 1
    #compute the number of pollinators, plants and interactions
    SA <- nrow(web_binary)
    SP <- ncol(web_binary)
    SI <- floor(sum(web_binary))
    #initialize the interaction matrix with minimum requirements
    web_opt <- matrix(0, nrow=SA, ncol=SP)
    web_opt[1,] <- 1
    web_opt[,1] <- 1
    web_opt[2,2] <- 1
    #counting the number of
    SI_left <- SI-SP-SA
    if(SI_left>0){
        #search the best possible location
        for(j in 1:SI_left){
            #compare all possible locations and the maximum one
            # Print the potential next position to see what exactly
            # is happening here:
            position_potential <- websearch_NODF(web_opt)
            nest_poten <- c()
            for(i in 1:nrow(position_potential)) {
                web_poten <- web_opt
                web_poten[position_potential[i,1],position_potential[i,2]] <- 1
                nest_poten[i] <- nestedness_NODF(web_poten)
            }
            position_the <- which.is.max(nest_poten)
            web_opt[position_potential[position_the,1],position_potential[position_the,2]] <- 1
        }
        return(nestedness_NODF(web_opt))
    }
    #this is to prevent the trivial case
    else{
        return(-1)
    }
}

#finds all possible positions to add an interaction
#input: web = mutualistic network
#output: all achievable positions of adding an interaction to the network
websearch_NODF <- function(web){
    SA <- nrow(web)
    SP <- ncol(web)
    domain <- web
    position <- which(domain == 1, arr.ind=T)
    position <- subset(position, position[,2] != 1)
    position <- subset(position, position[,1] != 1)
    boundary <- matrix(0, nrow=2*nrow(position),ncol=2)
    j=1
    #choose boundary points
    for(i in 1:nrow(position)){
        if(position[i,1]<nrow(domain)&&position[i,2]<ncol(domain))
            if(domain[position[i,1]+1,position[i,2]]+
               domain[position[i,1]-1,position[i,2]]+
               domain[position[i,1],position[i,2]+1]+
               domain[position[i,1],position[i,2]-1]<=3){
                boundary[j,1] <- position[i,1]+1
                boundary[j,2] <- position[i,2]
                boundary[j+1,1] <- position[i,1]
                boundary[j+1,2] <- position[i,2]+1
                j <- j+2
            }
    }
    #delete those with zero entries which entered as auxiliary in the first place
    keep <- c()
    for(i in 1:nrow(boundary)){
        if(boundary[i,1]+boundary[i,2]>0) keep <- append(keep,i)
    }
    boundary <- boundary[keep,]
    #choose true boundary points
    stay <- c()
    for(i in 1:nrow(boundary)){
        if(boundary[i,1]<SA&&boundary[i,2]<SP){
            if(domain[boundary[i,1]+1,boundary[i,2]]+
               domain[boundary[i,1]-1,boundary[i,2]]+
               domain[boundary[i,1],boundary[i,2]+1]+
               domain[boundary[i,1],boundary[i,2]-1]==2
           && domain[boundary[i,1],boundary[i,2]]==0){
                stay <- append(stay,i)
            }
        }
    }
    boundary <- boundary[stay,]
    return(boundary)
}

websearch_NODF_fast <- function(mtx){
    k = matrix(0,3,3)
    k[1,2] = 1
    k[2,1] = 1
    k[2,3] = 1
    k[3,2] = 1

    NodesA = nrow(mtx)
    NodesB = ncol(mtx)
    mtx2 <- 0.0*mtx
    mtx2[1:NodesA - 1,] <- mtx2[1:NodesA -1,] + mtx[2:NodesA,]
    mtx2[2:NodesA,] <- mtx2[2:NodesA,] + mtx[1:NodesA-1,]
    mtx2[,1:NodesB-1] <- mtx2[,1:NodesB-1] + mtx[,2:NodesB]
    mtx2[,2:NodesB] <- mtx2[,2:NodesB] + mtx[,1:NodesB-1]
    mtx3 <- matrix(as.numeric((mtx2 <= 4.0) * (mtx >= 0.5)), nrow(mtx), ncol(mtx))
    mtx4 <- 0.0*mtx3
    mtx4[2:NodesA,] <- mtx4[2:NodesA,] + mtx3[1:NodesA-1,]
    mtx4[,2:NodesB] <- mtx4[,2:NodesB] + mtx3[,1:NodesB-1]
    mtx5 = (mtx4 >= 0.5)
    mtx6 = matrix(as.numeric(mtx5 * (mtx == 0.0) * (mtx4 == 2.0)), nrow(mtx), ncol(mtx))
    posList = which(mtx6 == 1, arr.ind = T)
    return(posList)
}

# calculates the combined NODF statistic
# inputs: web = mutualistic network, raw NODF, maximum raw NODF
# output: the combined NODF statistic
comb_nest <- function(web,NODF,max_NODF){
    C <- sum(web)/(ncol(web)*nrow(web))
    S <- sqrt(ncol(web) * nrow(web) )
    out <- NODF / (max_NODF * C * log10(S))
    return(out)
}

# Computes the row marginal totals
marginal_total0 <- function(mtx){
    mt_0 <- rowSums(mtx, na.rm = FALSE, dims = 1)
    return(mt_0)
}

# Computes the column marginal totals
marginal_totalt <- function(mtx){
    mt_t <- colSums(mtx, na.rm = FALSE, dims = 1)
    return(mt_t)
}

# Assembles the list of marginal totals
compute_marginal_totals <- function(mtx){
    mt_0  <- marginal_total0(mtx)
    mt_t  <- marginal_totalt(mtx)
    return(list(mt_0, mt_t))
}

# Computes a list containing both Fill matrices
compute_fill_factors <- function(mtx){
    F0 <- mtx %*% t(mtx)
    Ft <- t(mtx) %*% mtx
    return(list(F0, Ft))
}

# Compute a list containing both degree matrices
compute_deg_mtx <- function (mtx) {
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    mt_0 <- marginal_total0(mtx)
    mt_t <- marginal_totalt(mtx)
    deg_mtx0 <- matrix(mt_0, nrow=length(mt_0),ncol=length(mt_0),byrow=TRUE)
    deg_mtxt <- matrix(mt_t, nrow=length(mt_t),ncol=length(mt_t),byrow=TRUE)
    return(list(deg_mtx0, deg_mtxt))
}

# Compute a list containing both degree minima matrices
compute_deg_minima <- function(mtx){
    my_ans <- compute_deg_mtx(mtx)
    DM0 <- my_ans[[1]]
    DMt <- my_ans[[2]]
    deg_min0 <- pmin(DM0, t(DM0))
    deg_mint <- pmin(DMt, t(DMt))
    return(list(deg_min0, deg_mint))
}

# Compute a list containing both negative delta matrices
compute_neg_deltas <- function(mtx){
    mt_0 <- marginal_total0(mtx)
    mt_t <- marginal_totalt(mtx)
    neg_delta0 = outer(mt_0, mt_0, FUN = ">")
    neg_deltat = outer(mt_t, mt_t, FUN = ">")
    return(list(neg_delta0, neg_deltat))
}

# Compute a list containing both the row and column sum values
compute_sums <- function(mtx){
    my_res <- compute_fill_factors(mtx)
    F0 <- my_res[[1]]
    Ft <- my_res[[2]]
    my_res <- compute_deg_minima(mtx)
    DM0 <- my_res[[1]]
    DMt <- my_res[[2]]
    my_res <- compute_neg_deltas(mtx)
    ND0 <- my_res[[1]]
    NDt <- my_res[[2]]
    n_paris0 = F0[ND0] / (DM0[ND0])
    n_parist = Ft[NDt] / (DMt[NDt])
    sum0 = sum(n_paris0)
    sumt = sum(n_parist)
    return(list(sum0, sumt))
}

# Assembles the list containing all the additional data required
# for fast_nodf computations
init_nodf <- function(mtx){
    MT <- compute_marginal_totals(mtx)
    Fill <- compute_fill_factors(mtx)
    DM <- compute_deg_minima(mtx)
    ND <- compute_neg_deltas(mtx)
    S <- compute_sums(mtx)
    return(list(MT, Fill, DM, ND, S))
}

get_contributions <- function(Fill, ND, DM, idx){
    A1 <- Fill[idx,][ND[idx,]] / (DM[idx,][ND[idx,]])
    A2 <- Fill[,idx][ND[,idx]] / (DM[,idx][ND[,idx]])
    return(sum(A1) + sum(A2))
}

# Efficient way to compute the nodf value of a matrix where the
# link at position pos = list(xpos, ypos) is removed.
nodf_one_link_removed <- function(mtx, pos, support_data){
    # Unpack and update all the support data:
    xpos <- pos[[1]]
    ypos <- pos[[2]]
    MT <- support_data[[1]]
    Fill <- support_data[[2]]
    DM <- support_data[[3]]
    ND <- support_data[[4]]
    S <- support_data[[5]]
    # print("B")
    # Unpack even futher to get a set a variables that we can work with:
    mt_0 <- MT[[1]]
    mt_t <- MT[[2]]
    # print("B1")
    F0 <- Fill[[1]]
    Ft <- Fill[[2]]
    # print("B2")
    DM0 <- DM[[1]]
    DMt <- DM[[2]]
    # print("B3")
    ND0 <- ND[[1]]
    NDt <- ND[[2]]
    # print("B4")
    S0 <- S[[1]]
    St <- S[[2]]
    # print("B5")
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    my_norm <- 0.5*((NodesA *(NodesA -1)) + (NodesB*(NodesB-1)))
    # print("B6")
    # Modify the matrix appropriately (remove a link):
    mtx[xpos, ypos] <- 0.0
    # print("C")

    # Compute old contributions
    old_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    old_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)

    # modify the marginal totals
    mt_0[xpos] <- mt_0[xpos] - 1
    mt_t[ypos] <- mt_t[ypos] - 1

    # Update the degree matrix:
    m0 <- matrix(mt_0[xpos], nrow=length(mt_0), ncol=1)
    mt <- matrix(mt_t[ypos], nrow=length(mt_t), ncol=1)

    # Update the degree minima:
    DM0[xpos,] <- pmin(m0, mt_0)
    DM0[,xpos] <- pmin(m0, mt_0)
    DMt[ypos,] <- pmin(mt, mt_t)
    DMt[,ypos] <- pmin(mt, mt_t)

    # Update negative deltas:
    ND0[xpos, ] <- (m0 > mt_0)
    ND0[, xpos] <- (m0 < mt_0)
    NDt[ypos, ] <- (mt > mt_t)
    NDt[, ypos] <- (mt < mt_t)

    # Update the fill factors
    F0[,xpos] <- F0[,xpos] - mtx[,ypos]
    F0[xpos,] <- F0[xpos,] - t(mtx[,ypos])
    F0[xpos,xpos] <- F0[xpos,xpos] - 1
    #
    Ft[, ypos] <- Ft[,ypos] - mtx[xpos,]
    Ft[ypos,] <- Ft[ypos,] - mtx[xpos,]
    Ft[ypos,ypos] <- Ft[ypos,ypos] - 1

    # Compute the new contributions:
    new_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    new_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)

    # Update the sums
    S0 <- S0 - old_contrib_0 + new_contrib_0
    St <- St - old_contrib_t + new_contrib_t

    # Compute the new nodf:
    nodf <- (S0+St) / my_norm

    # Pack up the results
    # Unpack even futher to get a set a variables that we can work with:
    MT <- list(mt_0, mt_t)
    Fill <- list(F0, Ft)
    DM <- list(DM0, DMt)
    ND <- list(ND0, NDt)
    S <- list(S0, St)
    # Unpack all the support data:
    support_data <- list(MT, Fill, DM, ND, S)
    return(list(nodf, mtx, support_data))
}

# Efficient way to compute the nodf value of a matrix where the
# link at position pos = list(xpos, ypos) is added:
nodf_one_link_added <- function(mtx, pos, support_data){
    # Unpack and update all the support data:
    xpos <- pos[[1]]
    ypos <- pos[[2]]
    MT <- support_data[[1]]
    Fill <- support_data[[2]]
    DM <- support_data[[3]]
    ND <- support_data[[4]]
    S <- support_data[[5]]
    # print("B")
    # Unpack even futher to get a set a variables that we can work with:
    mt_0 <- MT[[1]]
    mt_t <- MT[[2]]
    # print("B1")
    F0 <- Fill[[1]]
    Ft <- Fill[[2]]
    # print("B2")
    DM0 <- DM[[1]]
    DMt <- DM[[2]]
    # print("B3")
    ND0 <- ND[[1]]
    NDt <- ND[[2]]
    # print("B4")
    S0 <- S[[1]]
    St <- S[[2]]
    # print("B5")
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    my_norm <- 0.5*((NodesA *(NodesA -1)) + (NodesB*(NodesB-1)))
    # print("B6")
    # Modify the matrix appropriately (remove a link):
    mtx[xpos, ypos] <- 1.0
    # print("C")

    # Compute old contributions
    old_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    old_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)

    # modify the marginal totals
    mt_0[xpos] <- mt_0[xpos] + 1
    mt_t[ypos] <- mt_t[ypos] + 1

    # Update the degree matrix:
    m0 <- matrix(mt_0[xpos], nrow=length(mt_0), ncol=1)
    mt <- matrix(mt_t[ypos], nrow=length(mt_t), ncol=1)

    # Update the degree minima:
    DM0[xpos,] <- pmin(m0, mt_0)
    DM0[,xpos] <- pmin(m0, mt_0)
    DMt[ypos,] <- pmin(mt, mt_t)
    DMt[,ypos] <- pmin(mt, mt_t)

    # Update negative deltas:
    ND0[xpos, ] <- (m0 > mt_0)
    ND0[, xpos] <- (m0 < mt_0)
    NDt[ypos, ] <- (mt > mt_t)
    NDt[, ypos] <- (mt < mt_t)

    # Update the fill factors
    F0[,xpos] <- F0[,xpos] + mtx[,ypos]
    F0[xpos,] <- F0[xpos,] + t(mtx[,ypos])
    F0[xpos,xpos] <- F0[xpos,xpos] - 1
    #
    Ft[,ypos] <- Ft[,ypos] + t(mtx[xpos,])
    Ft[ypos,] <- Ft[ypos,] + mtx[xpos,]
    Ft[ypos,ypos] <- Ft[ypos,ypos] - 1

    # Compute the new contributions:
    new_contrib_0 <- get_contributions(F0, ND0, DM0, xpos)
    new_contrib_t <- get_contributions(Ft, NDt, DMt, ypos)
    #print(c(S0, old_contrib_0, new_contrib_0))
    #print(c(St, old_contrib_t, new_contrib_t))

    # Update the sums
    S0 <- S0 - old_contrib_0 + new_contrib_0
    St <- St - old_contrib_t + new_contrib_t
    #print(c(S0, St))

    # Compute the new nodf:
    nodf <- (S0+St) / my_norm

    # Pack up the results
    # Unpack even futher to get a set a variables that we can work with:
    MT <- list(mt_0, mt_t)
    Fill <- list(F0, Ft)
    DM <- list(DM0, DMt)
    ND <- list(ND0, NDt)
    S <- list(S0, St)
    # Unpack all the support data:
    support_data <- list(MT, Fill, DM, ND, S)
    return(list(nodf, mtx, support_data))
}

# Efficient way to compute the nodf value of a neighbor graph
nodf_neighbor <- function(mtx, oPos, zPos, support_data){
    my_res <- nodf_one_link_removed(mtx, oPos, support_data)
    mtx <- my_res[[2]]
    support_data <- my_res[[3]]
    my_res  <- nodf_one_link_added(mtx, zPos, support_data)
    # my_res = list(nodf, mtx, support_data)
    return(my_res)
}

get_valid_ones <- function(mtx){
    NodesA <- nrow(mtx)
    NodesB <- ncol(mtx)
    sub_mtx <- mtx[2:NodesA,2:NodesB]
    one_pos <- which(sub_mtx == 1, arr.ind = T)
    shift_mtx <- matrix(1, nrow = nrow(one_pos), ncol = 2)
    return(one_pos + shift_mtx)
}

get_zeros <- function(mtx){
    zero_pos <- which(mtx == 0, arr.ind = T)
    return(zero_pos)
}

accept_probability <- function(new_cost, old_cost, temp){
    if(new_cost < old_cost){
        result <- 1.0
    }
    else{
        a <- -1.0*(new_cost- old_cost) / temp
        result <- exp(a)
    }
    return(result)
}
