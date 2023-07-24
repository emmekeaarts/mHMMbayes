context("plot mHMM and mHMM.gamma")

### ==============================
### CREATE OUTPUT OBJECTS TO TEST
### ==============================

## general properties tested model
n_t <- 100
n <- 10
m <- 3
J = 11
burn_in = 5
n_dep <- 2
q_emiss <- c(4,2)

gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.2, 0.6, 0.2,
                  0.1, 0.2, 0.7), ncol = m, byrow = TRUE)
emiss_distr <- list(matrix(c(0.5, 0.5, 0.0, 0.0,
                             0.1, 0.1, 0.8, 0.0,
                             0.1, 0.1, 0.2, 0.6), nrow = m, ncol = q_emiss[1], byrow = TRUE),
                    matrix(c(0.7, 0.3,
                             0.9, 0.1,
                             0.8, 0.2), nrow = m, ncol = q_emiss[2], byrow = TRUE)
)

set.seed(4231)
data_sim1 <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss), gamma = gamma,
                     emiss_distr = emiss_distr, var_gamma = .5, var_emiss = c(.5, 0.5))
colnames(data_sim1$obs) <- c("subj", "output_1", "output_2")


# Fit the mHMM on 2 dep variable data
set.seed(3523)
out_2st_simb <- mHMM(s_data = data_sim1$obs,
                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                     start_val = c(list(gamma), emiss_distr),
                     mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)
n_subj <- out_2st_simb$input$n_subj
n_dep <- out_2st_simb$input$n_dep

xx_cont <- c(list(matrix(c(rep(1, n_subj),rnorm(n_subj)),
                         ncol = 2, nrow = n_subj)),
             rep(list(matrix(c(rep(1, n_subj),rnorm(n_subj)),
                             ncol = 2, nrow = n_subj)),n_dep))

xx_dich <- c(list(matrix(c(rep(1, n_subj), sample(c(0,1), n_subj, replace = T)),
                         ncol = 2, nrow = n_subj)),
             rep(list(matrix(c(rep(1, n_subj), sample(c(0,1), n_subj, replace = T)),
                             ncol = 2, nrow = n_subj)),n_dep))

# Fit the mHMM on 2 dep variable data with continuous covariates
set.seed(3523)
out_2st_simb_ccov <- mHMM(s_data = data_sim1$obs, xx = xx_cont,
                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                     start_val = c(list(gamma), emiss_distr),
                     mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)
# Fit the mHMM on 2 dep variable data with continuous covariates
set.seed(3523)
out_2st_simb_dcov <- mHMM(s_data = data_sim1$obs, xx = xx_dich,
                          gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                          start_val = c(list(gamma), emiss_distr),
                          mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)

### ==============================
## TESTING
### ==============================

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
  expect_silent(plot(out_2st_simb, component = "gamma", col = c("goldenrod", "steelblue", "indianred"), cat_lab = c("active", "not active")))
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

  # plot_pred
  expect_error(plot_pred(out_2st_simb), "covariate is required for predicting")
  expect_no_error(plot_pred(out_2st_simb_ccov, component = "gamma"))
  expect_no_error(plot_pred(out_2st_simb_dcov, component = "gamma"))
  expect_no_error(plot_pred(out_2st_simb_ccov, component = "emiss"))
  expect_no_error(plot_pred(out_2st_simb_dcov, component = "emiss"))
  expect_no_error(plot_pred(out_2st_simb_ccov, component = "emiss", dep = 1))
  expect_no_error(plot_pred(out_2st_simb_dcov, component = "emiss", dep = 1))
  expect_no_error(plot_pred(out_2st_simb_ccov, component = "emiss", dep = 2, burn_in = 4, cat_lab = c("active", "not active")))
  expect_no_error(plot_pred(out_2st_simb_ccov, component = "gamma", col = c("goldenrod", "steelblue", "indianred"), cat_lab = c("active", "not active")))
  expect_error(plot_pred(out_2st_simb_ccov, component = "a"), "component should be a string")
  expect_error(plot_pred(gamma1_subj), "should be from the class mHMM")

  # traceplot
  L <- list(out_2st_simb, out_2st_simb_ccov, out_2st_simb_dcov)
  expect_silent(traceplot(L))
  expect_silent(traceplot(L, component = "gamma"))
  expect_silent(traceplot(L, component = "emiss"))
  expect_silent(traceplot(L, component = "emiss", dep = 2))
  expect_silent(traceplot(L, component = "emiss", dep = 1))
  expect_silent(traceplot(L, component = "emiss", dep = 2, burn_in = 4, cat_lab = c("active", "not active")))
  expect_silent(traceplot(L, component = "gamma", col = c("goldenrod", "steelblue", "indianred"), cat_lab = c("active", "not active")))
  expect_error(traceplot(L, component = "a"), "component should be a string")
  expect_error(traceplot(gamma1_subj), "should be from the class mHMM")
})




