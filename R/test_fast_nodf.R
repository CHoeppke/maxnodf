source('./R/toolbox.R')
source('./R/hill_climb.R')
source('./R/greedySoler.R')

A <- matrix(0.0, 7, 5)
A[1,] <- 1.0
A[2,0:3] <- 1.0
A[3,0:3] <- 1.0
A[4,0:1] <- 1.0
A[5,0:1] <- 1.0
A[6,0:1] <- 1.0
A[7,0:1] <- 1.0
pos = list(7,5)
support_data <- init_nodf(A)
my_res <- nodf_one_link_added(A, pos, support_data)
# Update the matrix
nodf <- my_res[[1]]
A <- my_res[[2]]
support_data <- my_res[[3]]
# Test if the nodf value is correct:
nodf_2 <- nestedness_NODF(A)
# test if the support_data is still correct:
support_data_2 <- init_nodf(A)
#print(support_data[[5]][[1]] - support_data_2[[5]][[1]])
#print(support_data[[5]][[2]] - support_data_2[[5]][[2]])
#print(support_data[[4]][[1]] - support_data_2[[4]][[1]])
#print(support_data[[4]][[2]] - support_data_2[[4]][[2]])
#print(support_data[[3]][[1]] - support_data_2[[3]][[1]])
#print(support_data[[3]][[2]] - support_data_2[[3]][[2]])
#print(support_data[[2]][[1]] - support_data_2[[2]][[1]])
#print(support_data[[2]][[2]] - support_data_2[[2]][[2]])
#print(support_data[[1]][[1]] - support_data_2[[1]][[1]])
#print(support_data[[1]][[2]] - support_data_2[[1]][[2]])

# Actually test fast nodf greedy solving:
B <- greedy_solve(14, 13, 52)
print(B)
print(nestedness_NODF(B))
