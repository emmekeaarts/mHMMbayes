context("obtain gamma and emiss")

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
                     mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)


####################
## TESTING
###############

test_that("expected errors obtain gamma and emiss", {
  expect_error(obtain_gamma(out_2st_simb, level = 'a'), " should be set to either group or subject")
  ab <- c(2,3,4)
  expect_error(obtain_gamma(ab), "should be from the class mHMM")
  expect_error(obtain_gamma(out_2st_simb, burn_in = 10), "burn in period should be at least 2 points smaller")
  expect_error(obtain_emiss(out_2st_simb, level = 'a'), " should be set to either group or subject")
  expect_error(obtain_emiss(ab), "should be from the class mHMM")
  expect_error(obtain_emiss(out_2st_simb, burn_in = 10), "burn in period should be at least 2 points smaller")
})

test_that("output obtain gamma", {
  gamma1_g    <- obtain_gamma(out_2st_simb, level = 'group')
  gamma1_subj <- obtain_gamma(out_2st_simb, level = 'subject')
  # test dimensions
  expect_equal(dim(gamma1_g), c(m,m))
  expect_equal(length(gamma1_subj), n)
  expect_equal(dim(gamma1_subj[[1]]), c(m,m))
  expect_equal(dim(gamma1_subj[[1]]), dim(gamma1_subj[[n]]))
  # calculations
  expect_equal(as.vector(gamma1_g[2,]), c(0.304, 0.488, 0.197))
  expect_equal(as.vector(gamma1_subj[[1]][3,]), c(0.16, 0.22, 0.62))
  expect_equal(as.vector(gamma1_subj[[n]][1,]), c(0.611, 0.277, 0.112))
})

test_that("output obtain emiss", {
  emiss1_g    <- obtain_emiss(out_2st_simb, level = 'group')
  emiss1_subj <- obtain_emiss(out_2st_simb, level = 'subject')
  # test dimensions
  expect_equal(dim(emiss1_g[[1]]), c(m,q_emiss2[1]))
  expect_equal(dim(emiss1_g[[2]]), c(m,q_emiss2[2]))
  expect_equal(length(emiss1_subj[[1]]), n)
  expect_equal(length(emiss1_subj[[2]]), n)
  expect_equal(dim(emiss1_subj[[1]][[1]]), c(m,q_emiss2[1]))
  expect_equal(dim(emiss1_subj[[2]][[1]]), dim(emiss1_subj[[2]][[n]]))
  # calculations
  expect_equal(as.vector(emiss1_g[[1]][2,]), c(0.125, 0.128, 0.654, 0.084))
  expect_equal(as.vector(emiss1_g[[2]][2,]), c( 0.807, 0.193))
  expect_equal(as.vector(emiss1_subj[[1]][[1]][3,]), c(0.146, 0.098, 0.220, 0.537))
  expect_equal(as.vector(emiss1_subj[[2]][[n]][1,]), c(0.774, 0.226))
})
