context("Cpp implementation of nodf")
library(maxnodf)

test_that("testing nodf_cpp",{
    A <- matrix(0.0, 7, 5)
    A[1,] <- 1.0
    A[2,1:3] <- 1.0
    A[3,1:3] <- 1.0
    A[4,1:1] <- 1.0
    A[5,1:1] <- 1.0
    A[6,1:1] <- 1.0
    A[7,1:1] <- 1.0
    expect_equal(nodf_cpp(A), 0.709677419354839)
    B <- maxnodf(14, 13, 52, 0)
    expect_equal(nodf_cpp(B), 0.875739644970414)
})

test_that("testing nodf link addition"){
    A <- matrix(0.0, 7, 5)
    A[1,] <- 1.0
    A[2,1:3] <- 1.0
    A[3,1:3] <- 1.0
    A[4,1:1] <- 1.0
    A[5,1:1] <- 1.0
    A[6,1:1] <- 1.0
    A[7,1:1] <- 1.0

}
