context("simulating multilevel HMM data")

###### setting parameter values

n_t     <- 100
n       <- 10
m       <- 3
q_emiss <- 4
gamma   <- matrix(c(0.8, 0.1, 0.1,
                    0.2, 0.7, 0.1,
                    0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
                        0.1, 0.1, 0.8, 0.0,
                        0.0, 0.0, 0.1, 0.9), nrow = m, ncol = q_emiss, byrow = TRUE)

beta      <- rep(list(NULL), 2)
beta[[1]] <- matrix(c(0.5, 1.0,
                      -0.5, 0.5,
                      0.0, 1.0), byrow = TRUE, ncol = 2)
beta2      <- rep(list(NULL), 2)
beta2[[1]] <- matrix(c(0.5, 1.0,
                      -0.5, 0.5,
                      0.0, 1.0), byrow = TRUE, ncol = 2)
beta2[[2]] <- matrix(c(0.5, 1.0,
                       -0.5, 0.5,
                       0.0, 1.0), byrow = TRUE, ncol = 2)

xx_vec      <- rep(list(NULL),2)
xx_vec[[1]] <-  c(rep(0,5), rep(1,5))

xx_vec2      <- rep(list(NULL),2)
xx_vec2[[1]] <-  c(rep(0,5), rep(1,5))
xx_vec2[[2]] <-  c(rep(0,5), rep(1,5))


####################
## TESTING
###############

test_that("expected errors simulating data", {

  # wrong input of either gamma or the emission distribution
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma =
                          matrix(c(.3, .4, .2,
                                   .1, .1, .8,
                                   .3, .2, .5), byrow = TRUE, ncol = m),
                        emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1),
               "row of the transition probability matrix gamma should sum up to 1")
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                        emiss_distr = matrix(c(0.5, 0.5, 0.0, 0.0,
                                        0.1, 0.1, 0.8, 0.2,
                                        0.0, 0.0, 0.1, 0.9), nrow = m, ncol = q_emiss, byrow = TRUE),
                        var_gamma = 1, var_emiss = 1),
               "row of the emission distribution matrix should sum up to 1")
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma =
                          matrix(c(.5, .5,
                                   .5, .5), byrow = TRUE, ncol = 2),
                        emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1),
               paste("matrix gamma should be a", m, "by", m, "matrix"))
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                        emiss_distr = matrix(c(0.5, 0.5, 0.0, 0.0,
                                               0.1, 0.1, 0.8, 0.2), nrow = m-1, ncol = q_emiss, byrow = TRUE),
                        var_gamma = 1, var_emiss = 1),
               "rows of the emission distribution matrix should be equal")
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                        emiss_distr = matrix(c(0.5, 0.5, 0.0,
                                               0.1, 0.1, 0.8,
                                               0.1, 0.1, 0.8), nrow = m, ncol = q_emiss-1, byrow = TRUE),
                        var_gamma = 1, var_emiss = 1),
               "columns of the emission distribution matrix should be equal")

  # wrong input of covariates and/or beta
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec, var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified.")
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                        emiss_distr = emiss_distr, beta = beta, var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified.")
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec, beta = beta2,
                        var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified in one of the elements")
  expect_error(sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec2, beta = beta,
                        var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified in one of the elements")



})

set.seed(2432)
data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                  emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1)


data2 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                  emiss_distr = emiss_distr, beta = beta, xx_vec = xx_vec,
                  var_gamma = 1, var_emiss = 1)



