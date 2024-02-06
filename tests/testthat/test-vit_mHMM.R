context("obtaining state sequence using Viterbi")

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
emiss_distr <- list(matrix(c(0.45, 0.45, 0.05, 0.05,
                             0.1, 0.05, 0.8, 0.05,
                             0.1, 0.1, 0.2, 0.6), nrow = m, ncol = q_emiss[1], byrow = TRUE),
                    matrix(c(0.7, 0.3,
                             0.9, 0.1,
                             0.8, 0.2), nrow = m, ncol = q_emiss[2], byrow = TRUE)
)

set.seed(4231)
data_sim <- sim_mHMM(n_t = n_t, n = n, gen = list(m= m, n_dep = n_dep, q_emiss = q_emiss), gamma = gamma,
                     emiss_distr = emiss_distr, var_gamma = .5, var_emiss = c(.5, 0.5))
colnames(data_sim$obs) <- c("subj", "output_1", "output_2")

# Fit the mHMM on 2 dep variable data
set.seed(3523)
out_2st_simb <- mHMM(s_data = data_sim$obs,
                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                     start_val = c(list(gamma), emiss_distr),
                     mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)


# Count data
set.seed(0602)

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

# Simulate count data:
data_count <- sim_mHMM(n_t = n_t,
                     n = n,
                     data_distr = "count",
                     gen = list(m = m, n_dep = n_dep),
                     gamma = gamma,
                     emiss_distr = emiss_distr,
                     var_gamma = 0.1,
                     var_emiss = c(5,2),
                     return_ind_par = TRUE)

# correct specification
emiss_mu0 <- list(matrix(c(30, 70, 170), nrow = 1),
                  matrix(c(7, 8, 18), nrow = 1))
emiss_K0  <- list(1, 1)
emiss_V   <- list(rep(16, m), rep(4, m))
emiss_nu  <- list(0.1, 0.1)

manual_prior_emiss1 <- prior_emiss_count(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = emiss_mu0,
  emiss_K0 = emiss_K0,
  emiss_V =  emiss_V,
  emiss_nu = emiss_nu)

out_count <- mHMM(s_data = data_count$obs,
                  gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                  start_val = c(list(gamma), emiss_distr),data_distr = "count",emiss_hyp_prior = manual_prior_emiss1,
                  mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)


####################
## TESTING
###############

test_that("expected errors viterbi", {
  expect_error(vit_mHMM(out_2st_simb, s_data = data_sim$obs, burn_in = 10), "burn in period should be at least 2 points smaller")
  ab <- c(2,3,4)
  expect_error(vit_mHMM(ab, s_data = data_sim$obs), "should be from the class mHMM")
  expect_error(vit_mHMM(out_2st_simb, s_data = data_sim$obs[1:200,]), "number of subjects in the datasets")
})

test_that("output viterbi", {
  states1 <- vit_mHMM(out_2st_simb, s_data = data_sim$obs)
  expect_equal(dim(states1), c(n_t * n, 2))
  expect_equal(sort(unique(states1[,2])), c(1:m))
  expect_equal(sum(states1[,2]), 1890)
})

test_that("output count viterbi", {
  states1 <- vit_mHMM(out_count, s_data = data_count$obs)
  expect_equal(dim(states1), c(n_t * n, 2))
  expect_equal(sort(unique(states1[,2])), c(1:m))
  expect_equal(sum(states1[,2]), 1710)
})

