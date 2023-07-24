context("predict transition/emission probabilities")

## ===================================
## General properties tested model
## ===================================
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
data_sim <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss), gamma = gamma,
                     emiss_distr = emiss_distr, var_gamma = .5, var_emiss = c(.5, 0.5))
colnames(data_sim$obs) <- c("subj", "output_1", "output_2")

## ===================================
## Fit the mHMM on 2 dep variable data
## ===================================

set.seed(3523)
out_2st_simb <- mHMM(s_data = data_sim$obs,
                     gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                     start_val = c(list(gamma), emiss_distr),
                     mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)

## ===================================
## Create model with covariates
## ===================================
xx <- rep(list(matrix(1, ncol = 1, nrow = n)), (n_dep + 1))

for(i in 2:(n_dep + 1)){
  xx[[i]] <- cbind(xx[[i]], nonverbal_cov$std_CDI_change)
}
set.seed(3523)
out_2st_sim_cov1 <- mHMM(s_data = data_sim$obs, xx = xx,
                         gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                         start_val = c(list(gamma), emiss_distr),
                         mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)

xx_gamma_only <- rep(list(matrix(1, ncol = 1, nrow = n)), (n_dep + 1))
xx_gamma_only[[1]] <- cbind(xx_gamma_only[[i]], nonverbal_cov$std_CDI_change)
set.seed(3523)
out_2st_sim_cov2 <- mHMM(s_data = data_sim$obs, xx = xx_gamma_only,
                         gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                         start_val = c(list(gamma), emiss_distr),
                         mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)


## ===================================
## Testing
## ===================================

test_that("expected errors", {
  expect_error(pred_prob(out_2st_simb), "A covariate is required")
  expect_error(pred_prob(gamma), "should be from the class mHMM")
  expect_error(pred_prob(out_2st_simb, component = 'a'), " should be a string")
  ab <- c(2,3,4)
  expect_error(pred_prob(ab), "should be from the class mHMM")
  expect_error(pred_prob(out_2st_sim_cov1, component = "gamma"), "Please include a covariate")
  expect_error(pred_prob(out_2st_sim_cov2, component = "emiss"), "Please include a covariate")
})

test_that("output of pred_prob", {
  pre_gamma    <- pred_prob(out_2st_sim_cov2, component = "gamma")
  pre_emiss    <- pred_prob(out_2st_sim_cov1, component = "emiss")
  pre_emiss_list    <- pred_prob(out_2st_sim_cov1, component = "emiss", print.df = FALSE)

  # test dimensions
  expect_equal(dim(pre_gamma), c(m*n, m+2))
  expect_equal(dim(pre_emiss), c(m*n, q_emiss[1]+2))
  expect_equal(length(pre_emiss_list), n)
  expect_equal(dim(pre_emiss_list[[1]]), c(m, q_emiss[1]+2))
  # calculations
  expect_equal(unname(unlist(pre_gamma[2,-c(1,2)])), c(0.708, 0.147, 0.145))
  expect_equal(unname(unlist(pre_emiss[3,-c(1,2)])), c(0.171, 0.177, 0.159, 0.493))
  expect_equal(unname(unlist(pre_emiss_list[[1]][1,-c(1,2)])), c(0.478, 0.442, 0.028, 0.053))
})
