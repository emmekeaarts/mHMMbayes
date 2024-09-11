################
### CREATE OUTPUT OBJECTS TO TEST COUNT DATA
################

# Count data
set.seed(1402)

## general properties tested model
n_t <- 100
n <- 10
m <- 3
J <- 11
burn_in <- 5
n_dep <- 2

gamma   <- matrix(c(0.8, 0.1, 0.1,
                    0.2, 0.7, 0.1,
                    0.2, 0.2, 0.6), ncol = m, byrow = TRUE)

emiss_distr <- list(matrix(c(30, 70, 170), nrow = m),
                    matrix(c(7, 8, 18), nrow = m))

emiss_distr_log <- lapply(emiss_distr, function(q) log(q))
emiss_V         <- list(rep(16, m), rep(4, m))
emiss_V_log     <- var_to_logvar(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu = emiss_distr,
                                 var_emiss = emiss_V,
                                 byrow = FALSE)
emiss_V_log <- lapply(emiss_V_log, function(q) matrix(q, nrow = m))

# Simulate count data:
data_count1 <- sim_mHMM(n_t = n_t,
                        n = n,
                        data_distr = "count",
                        gen = list(m = m, n_dep = n_dep),
                        gamma = gamma,
                        emiss_distr = emiss_distr_log,
                        var_gamma = 0.1,
                        var_emiss = emiss_V_log,
                        return_ind_par = TRUE,
                        log_scale = TRUE)

# correct specification
emiss_mu0_log <- list(matrix(log(c(30, 70, 170)), ncol = m, byrow = TRUE),
                      matrix(log(c(7, 8, 18)), ncol = m, byrow = TRUE))
emiss_K0  <- list(rep(1, 1), rep(1, 1))

emiss_mu0_log_cov <- list(matrix(c(log(c(30, 70, 170)),
                                   1, -1, 0), ncol = m, byrow = TRUE),
                          matrix(c(log(c(7, 8, 18)),
                                   2, 0, -2), ncol = m, byrow = TRUE))
emiss_K0_cov  <- list(rep(1, 2), rep(1, 2))
emiss_V   <- list(rep(16, m), rep(4, m))
emiss_nu  <- list(0.1, 0.1)

emiss_V_log  <- var_to_logvar(gen = list(m = m, n_dep = n_dep),
                              emiss_mu = emiss_distr,
                              var_emiss = emiss_V,
                              byrow = FALSE)

manual_prior_emiss_log <- prior_emiss_count(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = emiss_mu0_log,
  emiss_K0 = emiss_K0,
  emiss_V =  emiss_V_log,
  emiss_nu = emiss_nu,
  log_scale = TRUE)

out_count <- mHMM(s_data = data_count1$obs,
                  gen = list(m = m, n_dep = n_dep),
                  start_val = c(list(gamma), emiss_distr), data_distr = "count", emiss_hyp_prior = manual_prior_emiss_log,
                  mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)


####################
## TESTING
###############

test_that("plotting functions don't throw unexpected errors", {
  gamma1_g <- obtain_gamma(out_count)
  gamma1_subj <- obtain_gamma(out_count, level = "subject")
  emmis1_g <- obtain_emiss(out_count)
  emmis1_subj <- obtain_emiss(out_count, level = "subject")

  # plot.mHMM
  expect_silent(plot(out_count))
  expect_silent(plot(out_count, component = "gamma"))
  expect_silent(plot(out_count, component = "emiss"))
  expect_silent(plot(out_count, component = "emiss", dep = 1))
  expect_silent(plot(out_count, component = "emiss", dep = 2, burn_in = 4))
  expect_silent(plot(out_count, component = "gamma", col = c("goldenrod", "steelblue", "indianred"),
                     cat_lab = c("active", "not active")))
  expect_silent(plot(out_count, component = "emiss", bty = "l"))
  expect_error(plot(out_count, component = "a"), "component should be a string")
  ab <- c(2,3,4)
  expect_error(plot.mHMM(ab), "should be from the class mHMM")
  expect_error(plot(out_count, burn_in = 10), "burn in period should be at least 2 points smaller")

  # plot.mHMM_gamma
  expect_silent(plot(gamma1_g))
  expect_silent(plot(gamma1_subj, subj_nr = 3))
  expect_silent(plot(gamma1_subj, subj_nr = 3, col = rep( c("goldenrod", "steelblue", "indianred"), each = 3)))
  expect_warning(plot(gamma1_g, subj_nr = 2), "only be specified when plotting the subject level")
  expect_error(plot(gamma1_subj), "specified with the input variable -subj_nr-")
  expect_error(plot.mHMM_gamma(out_count), "should be from the class mHMM_gamma")
  expect_error(plot.mHMM_gamma(emmis1_g), "should be from the class mHMM_gamma")
})
