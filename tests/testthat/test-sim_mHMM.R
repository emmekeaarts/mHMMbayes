context("simulating multilevel HMM data")

###### setting parameter values

n_t     <- 100
n       <- 10
m       <- 3
q_emiss <- 4

# Categorical data
gamma   <- matrix(c(0.8, 0.1, 0.1,
                    0.2, 0.7, 0.1,
                    0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
emiss_distr <- list(matrix(c(0.45, 0.45, 0.05, 0.05,
                             0.1, 0.05, 0.8, 0.05,
                             0.1, 0.1, 0.2, 0.6), nrow = m, ncol = q_emiss[1], byrow = TRUE))

beta      <- rep(list(NULL), 2)
beta[[1]] <- matrix(c(0.5, 1.0,
                      -0.5, 0.5,
                      0.0, 1.0), byrow = TRUE, ncol = 2)
beta2      <- rep(list(NULL), 2)
beta2[[1]] <- matrix(c(0.5, 1.0,
                      -0.5, 0.5,
                      0.0, 1.0), byrow = TRUE, ncol = 2)
beta2[[2]] <- matrix(c(0.5, 1.0, 0.0,
                       -0.5, 0.5, 1.0,
                       0.0, 1.0, 1.0), byrow = TRUE, ncol = 3)
beta3     <- rep(list(NULL), 2)
beta3[[1]] <- matrix(c(0.5, 1.0, 0.0,
                      -0.5, 0.5, 1.0,
                      0.0, 1.0, 0.0), byrow = TRUE, ncol = 3)
beta4     <- rep(list(NULL), 2)
beta4[[1]] <- matrix(c(0.5, 1.0, 0.0,
                       -0.5, 0.5, 1.0), byrow = TRUE, ncol = 3)

beta5     <- rep(list(NULL), 2)
beta5[[1]] <- matrix(c(0.5, 1.0,
                       -0.5, 0.5,
                       0.0, 1.0), byrow = TRUE, ncol = 2)
beta5[[2]] <- matrix(c(0.5, 1.0, 1.0, 0.0,
                       -0.5, 0.5, 1.0, 1.0,
                       0.0, 1.0, 1.0, 0.0), byrow = TRUE, ncol = 4)

beta6     <- rep(list(NULL), 2)
beta6[[1]] <- matrix(c(0.5, 1.0,
                       -0.5, 0.5,
                       0.0, 1.0), byrow = TRUE, ncol = 2)
beta6[[2]] <- matrix(c(0.5, 1.0, 0.0,
                       -0.5, 0.5, 1.0,
                       0.0, 1.0, 1.0,
                       0.0, 1.0, 1.0), byrow = TRUE, ncol = 3)

xx_vec      <- rep(list(NULL),2)
xx_vec[[1]] <-  c(rep(0,5), rep(1,5))

xx_vec2      <- rep(list(NULL),2)
xx_vec2[[1]] <-  c(rep(0,5), rep(1,5))
xx_vec2[[2]] <-  c(rep(0,5), rep(1,5))

xx_vec3     <- rep(list(NULL),2)
xx_vec3[[1]] <-  c(rep(0,4), rep(1,5))

xx_vec4      <- rep(list(NULL),3)
xx_vec4[[1]] <-  c(rep(0,5), rep(1,5))

beta7      <- rep(list(NULL), 3)
beta7[[2]] <- matrix(c(-1,
                       1,
                       1), byrow = TRUE, ncol = 1)
beta7[[3]] <- matrix(c(-0.1,
                       0.2,
                       0.2), byrow = TRUE, ncol = 1)

beta8      <- rep(list(NULL), 3)
beta8[[1]] <- matrix(c(0.5, 1.0,
                      -0.5, 0.5,
                      0.0, 1.0), byrow = TRUE, ncol = 2)

xx_vec7      <- rep(list(NULL),3)
xx_vec7[[2]] <-  c(rep(0,5), rep(1,5))
xx_vec7[[3]] <-  c(0.1, 0.0, 1.0, 0.5, 1.0, 0.1, 0.5, 1.0, 0.0, 0.5)

# Continuous data
n_dep_cont   <- 2

emiss_distr_cont <- list(matrix(c( 5, 1,
                              10, 1,
                              15, 1), nrow = m, byrow = TRUE),
                    matrix(c(0.5, 0.2,
                             1.0, 0.5,
                             2.0, 0.3), nrow = m, byrow = TRUE))
data_cont <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_cont),
                      data_distr = 'continuous', gamma = gamma, emiss_distr = emiss_distr_cont,
                      var_gamma = .5, var_emiss = c(.5, 0.01))

data_cont2 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_cont),
                       data_distr = 'continuous', gamma = gamma, emiss_distr = emiss_distr_cont, xx_vec = xx_vec4, beta = beta8,
                       var_gamma = .5, var_emiss = c(0, 0), return_ind_par = TRUE)

# Count data
n_dep_count   <- 2

emiss_distr_count <- list(matrix(c( 5,
                                   10,
                                   15), nrow = m, byrow = TRUE),
                         matrix(c(0.5,
                                  1.0,
                                  2.0), nrow = m, byrow = TRUE))

emiss_distr_count_log <- lapply(emiss_distr_count, function(q) log(q))

data_count <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_count),
                      data_distr = 'count', gamma = gamma, emiss_distr = emiss_distr_count,
                      var_gamma = .5, var_emiss = c(.5, 0.01))

data_count2 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_count),
                       data_distr = 'count', gamma = gamma, emiss_distr = emiss_distr_count_log, xx_vec = xx_vec4, beta = beta8,
                       var_gamma = .5, var_emiss = c(0, 0), return_ind_par = TRUE, log_scale = TRUE)


####################
## TESTING
###############

test_that("expected errors simulating data", {

  # wrong input of either gamma or the emission distribution
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma =
                          matrix(c(.3, .4, .2,
                                   .1, .1, .8,
                                   .3, .2, .5), byrow = TRUE, ncol = m),
                        emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1),
               "row of the transition probability matrix gamma should sum up to 1")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = list(matrix(c(0.5, 0.5, 0.0, 0.0,
                                        0.1, 0.1, 0.8, 0.2,
                                        0.0, 0.0, 0.1, 0.9), nrow = m, ncol = q_emiss, byrow = TRUE)),
                        var_gamma = 1, var_emiss = 1),
               "row of the emission distribution matrix should sum up to 1")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma =
                          matrix(c(.5, .5,
                                   .5, .5), byrow = TRUE, ncol = 2),
                        emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1),
               paste("matrix gamma should be a", m, "by", m, "matrix"))
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = list(matrix(c(0.5, 0.5, 0.0, 0.0,
                                               0.1, 0.1, 0.8, 0.2), nrow = m-1, ncol = q_emiss, byrow = TRUE)),
                        var_gamma = 1, var_emiss = 1),
               "rows of emission distribution matrix in element")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = list(matrix(c(0.5, 0.5, 0.0,
                                               0.1, 0.1, 0.8,
                                               0.1, 0.1, 0.8), nrow = m, ncol = q_emiss-1, byrow = TRUE)),
                        var_gamma = 1, var_emiss = 1),
               "columns of the emission distribution matrix should be")

  # wrong input of covariates and/or beta
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec, var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified.")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, beta = beta, var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified.")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec, beta = beta2,
                        var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified in one of the elements")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec2, beta = beta,
                        var_gamma = 1, var_emiss = 1),
               "Either only xx_vec or only beta is specified in one of the elements")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec3, beta = beta,
                        var_gamma = 1, var_emiss = 1),
               "length of the vectors in xx_vec should be equal")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec, beta = beta3,
                        var_gamma = 1, var_emiss = 1),
               "first element of beta to predict the transiton probability matrix gamma should be")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec, beta = beta4,
                        var_gamma = 1, var_emiss = 1),
               "first element of beta to predict the transiton probability matrix gamma should be")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec2, beta = beta5,
                        var_gamma = 1, var_emiss = 1),
               "second element of beta to predict the emission")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                        emiss_distr = emiss_distr, xx_vec = xx_vec2, beta = beta6,
                        var_gamma = 1, var_emiss = 1),
               "second element of beta to predict the emission")
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_count),
                        data_distr = 'count', gamma = gamma, emiss_distr = emiss_distr_count_log, xx_vec = xx_vec4, beta = beta8,
                        var_gamma = .5, var_emiss = c(0, 0), log_scale = FALSE),
               "Covariates have been used to predict the count emission distribution")

  # wrong input of between-subject variance
  expect_error(sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_count),
                        data_distr = 'count', gamma = gamma, emiss_distr = emiss_distr_count,
                        var_gamma = .5, var_emiss = list(matrix(.5), matrix(0.01))),
               "The number of rows of the between-subject variance for the emission distribution should be")
})


test_that("output simulated data", {
  set.seed(2432)
  data1 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                    emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1)


  data2 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = 1, q_emiss = q_emiss), gamma = gamma,
                    emiss_distr = emiss_distr, beta = beta, xx_vec = xx_vec,
                    var_gamma = 1, var_emiss = 1)

  # dimensions
  expect_equal(dim(data1$obs), c(n * n_t, 2))
  expect_equal(dim(data1$states), c(n * n_t, 2))
  expect_equal(unique(data1$obs[,1]), 1:n)
  expect_equal(sort(unique(data1$obs[,2])), 1:q_emiss)
  expect_equal(unique(data1$states[,1]), 1:n)
  expect_equal(sort(unique(data1$states[,2])), 1:m)
  expect_equal(dim(data2$obs), c(n * n_t, 2))
  expect_equal(dim(data2$states), c(n * n_t, 2))
  expect_equal(unique(data2$obs[,1]), 1:n)
  expect_equal(sort(unique(data2$obs[,2])), 1:q_emiss)
  expect_equal(unique(data2$states[,1]), 1:n)
  expect_equal(sort(unique(data2$states[,2])), 1:m)

  # calculations
  expect_equal(data1$obs[10:15,1], rep(1,6))
  expect_equal(data1$states[101:105,1], rep(2,5))
  expect_equal(data2$obs[10:15,1], rep(1,6))
  expect_equal(data2$states[101:105,1], rep(2,5))
})

test_that("output simulated count data", {
  set.seed(2432)
  data1 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_count),
                    data_distr = 'count', gamma = gamma, emiss_distr = emiss_distr_count,
                    var_gamma = .5, var_emiss = c(.5, 0.01))

  data2 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep_count),
                    data_distr = 'count', gamma = gamma, emiss_distr = emiss_distr_count_log, xx_vec = xx_vec4, beta = beta8,
                    var_gamma = .5, var_emiss = c(0.01, 0.001), return_ind_par = TRUE, log_scale = TRUE)

  # dimensions
  expect_equal(dim(data1$obs), c(n * n_t, 3))
  expect_equal(dim(data1$states), c(n * n_t, 2))
  expect_equal(unique(data1$obs[,1]), 1:n)
  expect_equal(sort(unique(data1$obs[,2])), 0:26)
  expect_equal(unique(data1$states[,1]), 1:n)
  expect_equal(sort(unique(data1$states[,2])), 1:m)
  expect_equal(dim(data2$obs), c(n * n_t, 3))
  expect_equal(dim(data2$states), c(n * n_t, 2))
  expect_equal(unique(data2$obs[,1]), 1:n)
  expect_equal(sort(unique(data2$obs[,2])), 0:25)
  expect_equal(unique(data2$states[,1]), 1:n)
  expect_equal(sort(unique(data2$states[,2])), 1:m)

  # calculations
  expect_equal(data1$obs[10:15,1], rep(1,6))
  expect_equal(data1$states[101:105,1], rep(2,5))
  expect_equal(data2$obs[10:15,1], rep(1,6))
  expect_equal(data2$states[101:105,1], rep(2,5))
})


