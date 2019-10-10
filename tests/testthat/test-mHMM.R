context("function mHMM and S3 print and summary methods")


################
### CREATE OUTPUT OBJECTS TO TEST
################

## general properties tested model
n_t <- 100
n <- 10
m <- 3
J = 11
burn_in = 5

##### create test data wit 1 dependent variable
q_emiss1 <- 4
n_dep1 <- 1

gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.2, 0.6, 0.2,
                  0.1, 0.2, 0.7), ncol = m, byrow = TRUE)
emiss_distr1 <- matrix(c(0.5, 0.5, 0.0, 0.0,
                        0.1, 0.1, 0.8, 0.0,
                        0.1, 0.1, 0.2, 0.6), nrow = m, ncol = q_emiss1, byrow = TRUE)
set.seed(4231)
data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss1, gamma = gamma,
                  emiss_distr = emiss_distr1, var_gamma = .5, var_emiss = .5)

# Fit the mHMM on 1 dep variable data
set.seed(3523)
out_2st_sim <- mHMM(s_data = data1$obs,
                    gen = list(m = m, n_dep = n_dep1, q_emiss = q_emiss1),
                    start_val = list(gamma, emiss_distr1),
                    mcmc = list(J = J, burn_in = burn_in))

##### create test data wit 2 dependent variable
n_dep2 <- 2
q_emiss2 <- c(4,2)
emiss_distr2 <- matrix(c(0.7, 0.3,
                         0.9, 0.1,
                         0.8, 0.2), nrow = m, ncol = q_emiss2[2], byrow = TRUE)
set.seed(4231)
data2 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss2[2], gamma = gamma,
                  emiss_distr = emiss_distr2, var_gamma = .5, var_emiss = .5)
data3 <- list(states = data1$states, obs = cbind(data1$obs, data2$obs[,2]))
colnames(data3$obs) <- c("subj", "output_1", "output_2")

# Fit the mHMM on 2 dep variable data
set.seed(3523)
out_2st_simb <- mHMM(s_data = data3$obs,
                     gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                     start_val = list(gamma, emiss_distr1, emiss_distr2),
                     mcmc = list(J = J, burn_in = burn_in))


##### Create model with covariates
xx <- rep(list(matrix(1, ncol = 1, nrow = n)), (n_dep2 + 1))
for(i in 2:(n_dep2 + 1)){
  xx[[i]] <- cbind(xx[[i]], nonverbal_cov$std_CDI_change)
}
set.seed(3523)
out_2st_sim_cov1 <- mHMM(s_data = data3$obs, xx = xx,
                     gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                     start_val = list(gamma, emiss_distr1, emiss_distr2),
                     mcmc = list(J = J, burn_in = burn_in))

xx_gamma_only <- rep(list(matrix(1, ncol = 1, nrow = n)), (n_dep2 + 1))
xx_gamma_only[[1]] <- cbind(xx_gamma_only[[i]], nonverbal_cov$std_CDI_change)
set.seed(3523)
out_2st_sim_cov2 <- mHMM(s_data = data3$obs, xx = xx_gamma_only,
                         gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                         start_val = list(gamma, emiss_distr1, emiss_distr2),
                         mcmc = list(J = J, burn_in = burn_in))

xx_gamma_only_V2 <- rep(list(NULL),  (n_dep2 + 1))
xx_gamma_only_V2[[1]] <- cbind(rep(1,n), nonverbal_cov$std_CDI_change)
set.seed(3523)
out_2st_sim_cov3 <- mHMM(s_data = data3$obs, xx = xx_gamma_only_V2,
                         gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                         start_val = list(gamma, emiss_distr1, emiss_distr2),
                         mcmc = list(J = J, burn_in = burn_in))

xx_dep1_only <- rep(list(matrix(1, ncol = 1, nrow = n)), (n_dep2 + 1))
xx_dep1_only[[2]] <- cbind(xx_dep1_only[[i]], nonverbal_cov$std_CDI_change)
set.seed(3523)
out_2st_sim_cov4 <- mHMM(s_data = data3$obs, xx = xx_dep1_only,
                         gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                         start_val = list(gamma, emiss_distr1, emiss_distr2),
                         mcmc = list(J = J, burn_in = burn_in))

xx_wrong1 <- rep(list(NULL),  (n_dep2 + 1))
for(i in 2:(n_dep2 + 1)){
  xx_wrong1[[i]] <- matrix(nonverbal_cov$std_CDI_change, ncol = 1)
}
xx_wrong2 <- rep(list(NULL),  (n_dep2 + 1))
xx_wrong2[[1]] <- matrix(nonverbal_cov$std_CDI_change, ncol = 1)
xx_wrong3 <- list(cbind(rep(1,n), nonverbal_cov$std_CDI_change))
xx_wrong4 <- cbind(rep(1,n), nonverbal_cov$std_CDI_change)
xx_wrong5 <- rep(list(NULL),  (n_dep2 + 1))
xx_wrong5[[1]] <- cbind(rep(1,n), c(rep(2,n/2), rep(1, n/2)))
xx_wrong6 <- rep(list(NULL),  (n_dep2 + 1))
xx_wrong6[[1]] <- data.frame(rep(1,n), as.factor(c(rep(2,n/2), rep(1, n/2))))

####################
## TESTING
###############


test_that("errors mHMM input", {
  expect_error(mHMM(s_data = data3$obs,
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, list(emiss_distr1, emiss_distr2)),
                    mcmc = list(J = J, burn_in = burn_in)), "number of elements in the list start_val")
  expect_error(mHMM(s_data = data.frame(data3$obs[,1], as.factor(data3$obs[,2]), data3$obs[,3]),
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, emiss_distr1, emiss_distr2),
                    mcmc = list(J = J, burn_in = burn_in)), "factorial variables")
  expect_error(mHMM(s_data = data3$obs, xx = xx_wrong1,
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, emiss_distr1, emiss_distr2),
                    mcmc = list(J = J, burn_in = burn_in)), "first column in each element of xx has to represent the intercept")
  expect_error(mHMM(s_data = data3$obs, xx = xx_wrong2,
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, emiss_distr1, emiss_distr2),
                    mcmc = list(J = J, burn_in = burn_in)), "first column in each element of xx has to represent the intercept")
  expect_error(mHMM(s_data = data3$obs, xx = xx_wrong3,
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, emiss_distr1, emiss_distr2),
                    mcmc = list(J = J, burn_in = burn_in)), "xx should be a list, with the number of elements")
  expect_error(mHMM(s_data = data3$obs, xx = xx_wrong4,
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, emiss_distr1, emiss_distr2),
                    mcmc = list(J = J, burn_in = burn_in)), "xx should be a list, with the number of elements")
  expect_error(mHMM(s_data = data3$obs, xx = xx_wrong5,
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, emiss_distr1, emiss_distr2),
                    mcmc = list(J = J, burn_in = burn_in)), "Dichotomous covariates")
  expect_error(mHMM(s_data = data3$obs, xx = xx_wrong6,
                    gen = list(m = m, n_dep = n_dep2, q_emiss = q_emiss2),
                    start_val = list(gamma, emiss_distr1, emiss_distr2),
                    mcmc = list(J = J, burn_in = burn_in)), "Factors currently cannot be used as covariates")
})


test_that("class inherit", {
  expect_s3_class(out_2st_sim, "mHMM")
  expect_s3_class(out_2st_simb, "mHMM")
})

test_that("print mHMM", {
  # 1 dependent var
  expect_output(print(out_2st_sim), paste("Number of subjects:", n))
  expect_output(print(out_2st_sim), paste(J, "iterations"))
  expect_output(print(out_2st_sim), paste("likelihood over all subjects:", -123))
  expect_output(print(out_2st_sim), paste("AIC over all subjects:", 277))
  expect_output(print(out_2st_sim), paste("states used:", m))
  expect_output(print(out_2st_sim), paste("dependent variables used:", n_dep1))
  # 2 dependent var
  expect_output(print(out_2st_simb), paste("Number of subjects:", n))
  expect_output(print(out_2st_simb), paste(J, "iterations"))
  expect_output(print(out_2st_simb), paste("likelihood over all subjects:", -174))
  expect_output(print(out_2st_simb), paste("AIC over all subjects:", 384))
  expect_output(print(out_2st_simb), paste("states used:", m))
  expect_output(print(out_2st_simb), paste("dependent variables used:", n_dep2))
})

test_that("summary mHMM", {
  # 1 dependent var
  expect_output(summary(out_2st_sim), "From state 1      0.688      0.160      0.140")
  expect_output(summary(out_2st_sim), "From state 3      0.197      0.271      0.525")
  expect_output(summary(out_2st_sim), "State 1      0.450      0.418")
  expect_output(summary(out_2st_sim), "State 3      0.148      0.195")
  # 2 dependent var
  expect_output(summary(out_2st_simb), "From state 1      0.714      0.152      0.124")
  expect_output(summary(out_2st_simb), "From state 3      0.214      0.224      0.553")
  expect_output(summary(out_2st_simb), "output_1")
  expect_output(summary(out_2st_simb), "output_2")
})


test_that("output PD_subj", {
  # test output dimensions
  expect_equal(length(out_2st_sim$PD_subj), n)
  expect_equal(dim(out_2st_sim$PD_subj[[1]]), c(J, m*q_emiss1 + m * m + 1))
  expect_equal(dim(out_2st_sim$PD_subj[[1]]), dim(out_2st_sim$PD_subj[[n]]))
  # start values
  expect_equal(as.vector(out_2st_sim$PD_subj[[2]][1,]), c(as.vector(t(emiss_distr1)),
                                                          as.vector(t(gamma)), NA))
  expect_equal(as.vector(out_2st_simb$PD_subj[[2]][1,]), c(as.vector(t(emiss_distr1)),
                                                           as.vector(t(emiss_distr2)),as.vector(t(gamma)), NA))
  # calculations
  expect_equal(as.numeric(out_2st_sim$PD_subj[[2]][10,(m*q_emiss1 + 2)]), 0.0652)

})

test_that("output emis_cov_bar", {
  # test output dimensions
  expect_true(is.character(out_2st_sim$emiss_cov_bar))
  expect_equal(dim(out_2st_sim_cov1$emiss_cov_bar[[1]]), c(J, m * (q_emiss2[1] - 1)))
  expect_equal(dim(out_2st_sim_cov1$emiss_cov_bar[[2]]), c(J, m * (q_emiss2[2] - 1)))
  expect_true(is.character(out_2st_sim_cov2$emiss_cov_bar))
  expect_true(is.character(out_2st_sim_cov3$emiss_cov_bar))
  expect_true(is.character(out_2st_sim_cov3$emiss_cov_bar))
  expect_equal(dim(out_2st_sim_cov4$emiss_cov_bar[[1]]), c(J, m * (q_emiss2[1] - 1)))
  expect_true(is.character(out_2st_sim_cov4$emiss_cov_bar[[2]]))
  # test calcultations

})

test_that("output emis_int_bar", {
  # test output dimensions
})

test_that("output emiss_int_subj", {
  # test output dimensions
})

test_that("output emiss_naccept", {
  # test output dimensions
})

test_that("output emiss_prob_bar", {
  # test output dimensions

  # start values
  expect_equal(as.vector(out_2st_sim$emiss_prob_bar[[1]][1,]), as.vector(t(emiss_distr1)))
  expect_equal(as.vector(out_2st_simb$emiss_prob_bar[[1]][1,]), as.vector(t(emiss_distr1)))
  expect_equal(as.vector(out_2st_simb$emiss_prob_bar[[2]][1,]), as.vector(t(emiss_distr2)))

})

test_that("output gamma_cov_bar", {
  # test output dimensions
})

test_that("output gamma_int_bar", {
  # test output dimensions
})

test_that("output gamma_int_subj", {
  # test output dimensions
})

test_that("output gamma_naccept", {
  # test output dimensions
})

test_that("output gamma_prob_bar", {
  # test output dimensions

  # start values
  expect_equal(as.vector(out_2st_sim$gamma_prob_bar[1,]), as.vector(t(gamma)))
  expect_equal(as.vector(out_2st_simb$gamma_prob_bar[1,]), as.vector(t(gamma)))
  # calculations
  expect_equal(as.numeric(out_2st_sim$gamma_prob_bar[11,1]), 0.6769550)
  expect_equal(as.numeric(out_2st_sim$gamma_prob_bar[10,8]), 0.3251220)
})

test_that("output input", {
  # test output dimensions
})








