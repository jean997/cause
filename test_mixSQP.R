library(mixsqp)
library(Rmosek)
library(ashr)

#library(dplyr)
#library(cause)
#set.seed(1)
#X1 <- read_tsv("test_data/giant_height_summary_statistics.tsv.gz")
#X2 <- read_tsv("test_data/vanderHarst_cad_summary_statistics.tsv.gz")

#X <- gwas_format_cause(X1, X2)

#varlist <- sample(X$snp, size=1e5, replace=FALSE)

#X <- filter(X, snp %in% varlist) %>% cause:::new_cause_data()
#saveRDS(X, file="X.RDS")
X <- readRDS("X.RDS")
system.time(
fit1 <- with(X, ash(betahat = beta_hat_1, sebetahat = seb1,
                    mixcompdist = "normal", prior="nullbiased",
                    optmethod="mixIP"))
)

system.time(
fit2 <- with(X, ash(betahat = beta_hat_1, sebetahat = seb1,
                    mixcompdist = "normal", prior="nullbiased",
                    optmethod="mixSQP"))
)
all(fit1$fitted_g$sd == fit2$fitted_g$sd)
#[1] TRUE

cbind(fit1$fitted_g$pi, fit2$fitted_g$pi)
# [,1]         [,2]
# [1,] 1.774469e-01 2.656528e-12
# [2,] 1.158056e-10 0.000000e+00
# [3,] 1.174241e-10 0.000000e+00
# [4,] 1.208374e-10 0.000000e+00
# [5,] 1.282829e-10 0.000000e+00
# [6,] 1.457089e-10 0.000000e+00
# [7,] 1.932258e-10 0.000000e+00
# [8,] 3.802649e-10 9.352648e-01
# [9,] 5.018178e-02 0.000000e+00
# [10,] 5.574481e-01 0.000000e+00
# [11,] 4.853750e-10 0.000000e+00
# [12,] 1.309805e-01 0.000000e+00
# [13,] 4.751587e-02 0.000000e+00
# [14,] 2.727578e-02 0.000000e+00
# [15,] 5.528548e-03 0.000000e+00
# [16,] 3.622427e-03 0.000000e+00
# [17,] 1.407412e-09 6.464525e-02
# [18,] 6.951089e-10 0.000000e+00
# [19,] 5.527904e-12 0.000000e+00
# [20,] 5.427889e-13 0.000000e+00
# [21,] 1.754670e-13 0.000000e+00
# [22,] 8.957492e-14 0.000000e+00
# [23,] 5.699670e-14 0.000000e+00

max(abs(fit1$fitted_g$pi - fit2$fitted_g$pi))

#mix_grid <- variance_pair_candidates(X, optmethod="mixIP")
#K <- nrow(mix_grid)
#rho <- 0
#matrix_llik1 <- loglik_mat(rho, 0, 0, 0,
#                           mix_grid$S1, mix_grid$S2,
#                           X$beta_hat_1,
#                           X$beta_hat_2,
#                           X$seb1,
#                           X$seb2)
#matrix_llik = matrix_llik1 - apply(matrix_llik1, 1, max)
#matrix_lik = exp(matrix_llik)
#saveRDS(matrix_lik, file="matrix_lik.RDS")


matrix_lik <- readRDS("matrix_lik.RDS")
K <- ncol(matrix_lik)
null_wt <- 10
system.time(
  w_res1 <- ashr:::mixIP(matrix_lik =matrix_lik,
                     prior=c(null_wt, rep(1, K-1)),
                     weights=rep(1, nrow(matrix_lik)))
)
# user  system elapsed
# 27.285   1.088   8.006

system.time(
  w_res2 <- ashr:::mixSQP(matrix_lik =matrix_lik,
                      prior=c(null_wt, rep(1, K-1)),
                      weights=rep(1, nrow(matrix_lik)))
)
# Error in mixsqp_rcpp(L, w, x0, convtol.sqp, convtol.activeset, zero.threshold.solution,  :
#                        Step size is too small; consider relaxing convergence criteria
#                      Timing stopped at: 215.5 5.083 131.1


matrix_lik <- matrix_lik[1:50000,]
system.time(
  w_res1 <- ashr:::mixIP(matrix_lik =matrix_lik,
                         prior=c(null_wt, rep(1, K-1)),
                         weights=rep(1, nrow(matrix_lik)))
)
# user  system elapsed
# 12.842   0.571   3.758
system.time(
  w_res2 <- ashr:::mixSQP(matrix_lik =matrix_lik,
                          prior=c(null_wt, rep(1, K-1)),
                          weights=rep(1, nrow(matrix_lik)))
)
# user  system elapsed
# 0.945   0.023   0.750

w_res1$pihat
# [1] 2.154525e-01 3.944225e-10 2.274752e-01 3.372971e-02 3.191436e-12 8.426919e-13 5.594566e-13
# [8] 3.046829e-12 9.408704e-12 8.706004e-11 1.812455e-01 1.224553e-01 4.778707e-12 1.237769e-12
# [15] 6.740887e-13 1.386815e-12 6.291435e-12 2.908199e-11 1.357813e-01 1.622843e-03 9.308698e-12
# [22] 2.155277e-12 9.451587e-13 1.667971e-12 3.400467e-12 9.932893e-12 2.856614e-11 2.021765e-10
# [29] 1.391612e-02 5.316106e-02 3.362397e-12 4.472526e-03 1.010816e-12 1.380213e-12 1.772973e-12
# [36] 2.311095e-12 1.546197e-12 2.054166e-12 3.102703e-03 1.708751e-04 1.568014e-12 3.345721e-12
# [43] 6.057098e-03 1.357375e-03 8.386978e-13 4.098344e-13 3.229903e-13 2.314049e-13 7.057686e-14
# [50] 7.922073e-14 9.256059e-14 1.399094e-13 1.034357e-13 6.446252e-14 4.453905e-14 3.239196e-14
# [57] 2.504603e-14 2.716369e-14 2.999194e-14 3.893903e-14 3.541981e-14 2.749726e-14 2.111025e-14
#[64] 1.635950e-14
w_res2$pihat
# [1]  2.939303e-01  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  2.391859e-09
# [7]  0.000000e+00  2.525254e-11  0.000000e+00  0.000000e+00  3.027992e-01  4.030501e-01
# [13] -3.337642e-17  0.000000e+00  0.000000e+00  0.000000e+00  4.408763e-18  4.202175e-18
# [19]  0.000000e+00  0.000000e+00 -8.699196e-19  0.000000e+00  0.000000e+00  0.000000e+00
# [25]  1.866865e-15  1.573426e-15  1.971372e-15  3.017748e-15  2.479745e-15 -7.846053e-10
# [31]  0.000000e+00  0.000000e+00  6.193780e-12  0.000000e+00  0.000000e+00  0.000000e+00
# [37]  0.000000e+00  4.546087e-10 -7.375361e-11 -2.955745e-10  0.000000e+00  0.000000e+00
# [43]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  7.122128e-10
# [49]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [55]  0.000000e+00 -9.722126e-10  5.758730e-12  0.000000e+00 -5.766104e-11  1.646584e-10
# [61]  0.000000e+00  0.000000e+00 -2.584939e-26  4.722525e-10
w_res1$pihat - w_res2$pihat
# [1] -7.847781e-02  3.944225e-10  2.274752e-01  3.372971e-02  3.191436e-12 -2.391016e-09
# [7]  5.594566e-13 -2.220572e-11  9.408704e-12  8.706004e-11 -1.215537e-01 -2.805948e-01
# [13]  4.778741e-12  1.237769e-12  6.740887e-13  1.386815e-12  6.291431e-12  2.908199e-11
# [19]  1.357813e-01  1.622843e-03  9.308699e-12  2.155277e-12  9.451587e-13  1.667971e-12
# [25]  3.398600e-12  9.931320e-12  2.856416e-11  2.021735e-10  1.391612e-02  5.316106e-02
# [31]  3.362397e-12  4.472526e-03 -5.182963e-12  1.380213e-12  1.772973e-12  2.311095e-12
# [37]  1.546197e-12 -4.525545e-10  3.102703e-03  1.708754e-04  1.568014e-12  3.345721e-12
# [43]  6.057098e-03  1.357375e-03  8.386978e-13  4.098344e-13  3.229903e-13 -7.119814e-10
# [49]  7.057686e-14  7.922073e-14  9.256059e-14  1.399094e-13  1.034357e-13  6.446252e-14
# [55]  4.453905e-14  9.722450e-10 -5.733684e-12  2.716369e-14  5.769103e-11 -1.646194e-10
# [61]  3.541981e-14  2.749726e-14  2.111025e-14 -4.722361e-10
