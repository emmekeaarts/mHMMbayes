context("obtaining state sequence using Viterbi")

## general properties tested model
n_t <- 100
n <- 10
m <- 3
J = 11
burn_in = 5

##### create test data wit 1 dependent variable
n_dep2 <- 2
q_emiss2 <- c(4,2)

gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.2, 0.6, 0.2,
                  0.1, 0.2, 0.7), ncol = m, byrow = TRUE)
emiss_distr1 <- matrix(c(0.5, 0.5, 0.0, 0.0,
                         0.1, 0.1, 0.8, 0.0,
                         0.1, 0.1, 0.2, 0.6), nrow = m, ncol = q_emiss2[1], byrow = TRUE)
emiss_distr2 <- matrix(c(0.7, 0.3,
                         0.9, 0.1,
                         0.8, 0.2), nrow = m, ncol = q_emiss2[2], byrow = TRUE)

# (small workaround as cannot simulate multivariate data yet)
set.seed(4231)
data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss2[1], gamma = gamma,
                  emiss_distr = emiss_distr1, var_gamma = .5, var_emiss = .5)

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


####################
## TESTING
###############

test_that("expected errors viterbi", {
  expect_error(vit_mHMM(out_2st_simb, s_data = data3$obs, burn_in = 10), "burn in period should be at least 2 points smaller")
  ab <- c(2,3,4)
  expect_error(vit_mHMM(ab, s_data = data3$obs), "should be from the class mHMM")
  expect_error(vit_mHMM(out_2st_simb, s_data = data3$obs[1:200,]), "number of subjects in the datasets")
})

test_that("output viterbi", {
  states1 <- vit_mHMM(out_2st_simb, s_data = data3$obs)
  expect_equal(dim(states1), c(n_t, n))
  expect_equal(sort(unique(states1[,1])), c(1:m))
  expect_equal(sort(unique(states1[,n])), c(1:m))
  expect_equal(sum(states1[,2]), 177)
  expect_equal(sum(states1[2,]), 18)
})
