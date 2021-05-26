context("plot mHMM and mHMM.gamma")

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
                  emiss_distr = list(emiss_distr1=emiss_distr1), var_gamma = .5, var_emiss = .5)

set.seed(4231)
data2 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss2[2], gamma = gamma,
                  emiss_distr = list(emiss_distr2=emiss_distr2), var_gamma = .5, var_emiss = .5)
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

test_that("plotting functions don't throw unexpected errors", {
  gamma1_g <- obtain_gamma(out_2st_simb)
  gamma1_subj <- obtain_gamma(out_2st_simb, level = "subject")
  emmis1_g <- obtain_emiss(out_2st_simb)

  # plot.mHMM
  expect_silent(plot(out_2st_simb))
  expect_silent(plot(out_2st_simb, component = "gamma"))
  expect_silent(plot(out_2st_simb, component = "emiss"))
  expect_silent(plot(out_2st_simb, component = "emiss", dep = 1))
  expect_silent(plot(out_2st_simb, component = "emiss", dep = 2, burn_in = 4, cat_lab = c("active", "not active")))
  expect_silent(plot(out_2st_simb, component = "gamma", col = c("goldenrod", "steelblue", "indianred"),
                     cat_lab = c("active", "not active")))
  expect_silent(plot(out_2st_simb, component = "emiss", bty = "l"))
  expect_error(plot(out_2st_simb, component = "a"), "component should be a string")
  ab <- c(2,3,4)
  expect_error(plot.mHMM(ab), "should be from the class mHMM")
  expect_error(plot(out_2st_simb, burn_in = 10), "burn in period should be at least 2 points smaller")

  # plot.mHMM_gamma
  expect_silent(plot(gamma1_g))
  expect_silent(plot(gamma1_subj, subj_nr = 3))
  expect_silent(plot(gamma1_subj, subj_nr = 3, col = rep( c("goldenrod", "steelblue", "indianred"), each = 3)))
  expect_warning(plot(gamma1_g, subj_nr = 2), "only be specified when plotting the subject level")
  expect_error(plot(gamma1_subj), "specified with the input variable -subj_nr-")
  expect_error(plot.mHMM_gamma(out_2st_simb), "should be from the class mHMM_gamma")
  expect_error(plot.mHMM_gamma(emmis1_g), "should be from the class mHMM_gamma")
})




