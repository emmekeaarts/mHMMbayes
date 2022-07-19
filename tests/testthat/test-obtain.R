context("obtain gamma and emiss")

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
data_sim <- sim_mHMM(n_t = n_t, n = n, m = m, n_dep = n_dep, q_emiss = q_emiss, gamma = gamma,
                     emiss_distr = emiss_distr, var_gamma = .5, var_emiss = c(.5, 0.5))
colnames(data_sim$obs) <- c("subj", "output_1", "output_2")

# Fit the mHMM on 2 dep variable data
set.seed(3523)
out_2st_simb <- mHMM(s_data = data_sim$obs,
                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                     start_val = c(list(gamma), emiss_distr),
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
  expect_equal(as.vector(gamma1_g[2,]), c(0.231,  0.553, 0.216))
  expect_equal(as.vector(gamma1_subj[[1]][3,]), c(0.235, 0.250, 0.511))
  expect_equal(as.vector(gamma1_subj[[n]][1,]), c(0.605, 0.212, 0.163))
})

test_that("output obtain emiss", {
  emiss1_g    <- obtain_emiss(out_2st_simb, level = 'group')
  emiss1_subj <- obtain_emiss(out_2st_simb, level = 'subject')
  # test dimensions
  expect_equal(dim(emiss1_g[[1]]), c(m,q_emiss[1]))
  expect_equal(dim(emiss1_g[[2]]), c(m,q_emiss[2]))
  expect_equal(length(emiss1_subj[[1]]), n)
  expect_equal(length(emiss1_subj[[2]]), n)
  expect_equal(dim(emiss1_subj[[1]][[1]]), c(m,q_emiss[1]))
  expect_equal(dim(emiss1_subj[[2]][[1]]), dim(emiss1_subj[[2]][[n]]))
  # calculations
  expect_equal(as.vector(emiss1_g[[1]][2,]), c(0.109, 0.139, 0.668, 0.097))
  expect_equal(as.vector(emiss1_g[[2]][2,]), c( 0.824, 0.176))
  expect_equal(as.vector(emiss1_subj[[1]][[1]][3,]), c(0.036, 0.072, 0.223, 0.660))
  expect_equal(as.vector(emiss1_subj[[2]][[n]][1,]), c(0.484, 0.516))
})
