context("plot mHMM, mHMM.gamma and mHMM.emiss")

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
emiss_distr <- list(matrix(c( 5, 1,
                              10, 1,
                              15, 1), nrow = m, byrow = TRUE),
                    matrix(c(0.5, 0.1,
                             1.0, 0.2,
                             2.0, 0.1), nrow = m, byrow = TRUE))

# (small workaround as cannot simulate multivariate data yet)
set.seed(4231)
data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss2[1], gamma = gamma,
                  emiss_distr = list(emiss_distr1=emiss_distr1), var_gamma = .5, var_emiss = .5)

set.seed(4231)
data2 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss2[2], gamma = gamma,
                  emiss_distr = list(emiss_distr2=emiss_distr2), var_gamma = .5, var_emiss = .5)

data3 <- list(states = data1$states, obs = cbind(data1$obs, data2$obs[,2]))
colnames(data3$obs) <- c("subj", "output_1", "output_2")

set.seed(4231)
data_cont <- sim_mHMM(n_t = n_t, n = n, m = m, n_dep = n_dep2, data_distr = c('continuous','continuous'),
                      gamma = gamma, emiss_distr = emiss_distr, var_gamma = .1, var_emiss = c(.5, 0.01))
# Specify hyper-prior for the continuous emission distribution
hyp_pr <- list(
                emiss_mu0 = list(matrix(c(3,7,17), nrow = 1), matrix(c(0.7, 0.8, 1.8), nrow = 1)),
                emiss_K0  = list(1, 1),
                emiss_nu  = list(1, 1),
                emiss_V   = list(rep(2, m), rep(1, m)),
                emiss_a0  = list(rep(1, m), rep(1, m)),
                emiss_b0  = list(rep(1, m), rep(1, m)))

# Fit the mHMM_cont on 2 dep variable data:
out_3st_cont_sim <- mHMM_cont(s_data = data_cont$obs,
                              gen = list(m = m, n_dep = n_dep2),
                              start_val = c(list(gamma), emiss_distr),
                              emiss_hyp_prior = hyp_pr,
                              mcmc = list(J = J, burn_in =burn_in, show_progress = FALSE))
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
  emiss1_g_cont <- obtain_emiss(out_3st_cont_sim)
  emiss1_subj_cont<-obtain_emiss(out_3st_cont_sim, level = "subject")
  emiss1_g<- obtain_emiss(out_2st_simb)
  emiss1_subj<-obtain_emiss(out_2st_simb, level = "subject")


  # plot.mHMM
  expect_silent(plot(out_2st_simb))
  expect_silent(plot(out_2st_simb, component = "gamma"))
  expect_silent(plot(out_2st_simb, component = "emiss"))
  expect_silent(plot(out_2st_simb, component = "emiss", dep = 1))
  expect_silent(plot(out_2st_simb, component = "emiss", dep = 2, burn_in = 4, cat_lab = list(c("one", "two","three","four"),c("one","two"))))
  expect_silent(plot(out_2st_simb, component = "gamma", col = c("goldenrod", "steelblue", "indianred"),
                     cat_lab = c("active", "not active")))
  expect_silent(plot(out_2st_simb, component = "emiss", bty = "l"))
  expect_error(plot(out_2st_simb, component = "emiss", dep = 2, burn_in = 4, cat_lab = list(c("one", "two","three"))),"cat_lab should be a list with n_dep")
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
  expect_error(plot.mHMM_gamma(emiss1_g), "should be from the class mHMM_gamma")

  #plot.mHMM_emiss
  expect_silent(plot(emiss1_g))
  expect_silent(plot(emiss1_subj, subj_nr = 3))
  expect_silent(plot.mHMM_emiss(emiss1_subj_cont, subj_nr = 3))
  expect_silent(plot.mHMM_emiss(emiss1_g_cont))
  expect_silent(plot(emiss1_subj, subj_nr = 3,col=c("orange","steelblue","indianred","goldenrod")))
  expect_silent(plot.mHMM_emiss(emiss1_g_cont,col=c("red","blue")))


  expect_warning(plot(emiss1_g, subj_nr = 2), "only be specified when plotting the subject level")
  expect_warning(plot.mHMM_emiss(emiss1_g_cont, subj_nr = 2), "only be specified when plotting the subject level")

  expect_error(plot(emiss1_subj), "specified with the input variable -subj_nr-")
  expect_error(plot.mHMM_emiss(emiss1_subj_cont), "specified with the input variable -subj_nr-")


  expect_error(plot.mHMM_emiss(out_2st_simb), "should be from the mHMM_emiss")
  expect_error(plot.mHMM_emiss(gamma1_g), "should be from the mHMM_emiss")

})




