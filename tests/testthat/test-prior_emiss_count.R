context("prior emiss count")


################
### CREATE OUTPUT OBJECTS TO TEST
################

## general properties tested model
n_t <- 100
n <- 10
m <- 3
J <- 11
burn_in <- 5
n_dep <- 2
n_xx_emiss <- c(1,1)

gamma   <- matrix(c(0.8, 0.1, 0.1,
                    0.2, 0.7, 0.1,
                    0.2, 0.2, 0.6), ncol = m, byrow = TRUE)

emiss_distr <- list(matrix(c(30, 70, 170), nrow = m),
                    matrix(c(7, 8, 18), nrow = m))

# Simulate count data:
data_sim <- sim_mHMM(n_t = n_t,
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

emiss_mu0_log <- lapply(emiss_mu0, function(q) log(q))
emiss_V_log <- obtain_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0,
                             var_emiss = emiss_V,
                             byrow = TRUE)

emiss_mu0_log_cov <- lapply(emiss_mu0, function(q) log(rbind(q,rep(1,m))))
emiss_K0_cov  <- list(rep(1,2), rep(1,2))



manual_prior_emiss1 <- prior_emiss_count(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = emiss_mu0,
  emiss_K0 = emiss_K0,
  emiss_V =  emiss_V,
  emiss_nu = emiss_nu)

manual_prior_emiss2 <- prior_emiss_count(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = emiss_mu0_log,
  emiss_K0 = emiss_K0,
  emiss_V =  emiss_V_log,
  emiss_nu = emiss_nu,
  log_scale = TRUE)

manual_prior_emiss3 <- prior_emiss_count(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = emiss_mu0_log_cov,
  emiss_K0 = emiss_K0_cov,
  emiss_V =  emiss_V_log,
  emiss_nu = emiss_nu,
  log_scale = TRUE,
  n_xx_emiss = c(1,1))


# wrong specifications
m_wrong <- 2
n_dep_wrong <- 3

emiss_mu0_wrong1 <- emiss_mu0[[1]]
emiss_mu0_wrong2 <- c(emiss_mu0,emiss_mu0[2])
emiss_mu0_wrong3 <- lapply(emiss_mu0, function(q) list(q))
emiss_mu0_wrong4 <- lapply(emiss_mu0, function(q) matrix(q[-1], nrow = 1))

emiss_K0_wrong1 <- rep(list(c(1)), n_dep + 1)
emiss_K0_wrong2 <- rep((c(1)), n_dep)
emiss_nu_wrong1 <- list(c(5), c(4), c(5))
emiss_nu_wrong2 <- (c(5, 4, 5))
emiss_V_wrong1 <- list(diag(5, m+1),
                       diag(4, m-1))

manual_prior_emiss_wrongm <- manual_prior_emiss1
manual_prior_emiss_wrongm$gen$m <- 2
manual_prior_emiss_wrongn_dep <- manual_prior_emiss1
manual_prior_emiss_wrongn_dep$gen$n_dep <- 3


####################
## TESTING
###############

test_that("errors prior_emiss_count input", {
  expect_error(prior_emiss_count(gen = list(m = m_wrong, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "The matrix relating to dependent variable")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep_wrong),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_mu0 should be a list containing")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0_wrong1, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_mu0 should be a list containing")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0_wrong2, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_mu0 should be a list containing")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0_wrong3, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_mu0 should be a list containing")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0_wrong4, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "The matrix relating to dependent variable")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0_wrong1,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_K0 should be a list containing n_dep, here 2")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0_wrong2,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_K0 should be a list containing n_dep, here 2")
  expect_warning(expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                                emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                                emiss_nu = emiss_nu_wrong1, emiss_V = emiss_V),
                              "emiss_nu should be a list containing n_dep"))
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu_wrong2, emiss_V = emiss_V),
               "emiss_nu should be a list containing n_dep")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V_wrong1),
               "emiss_V should be a list containing n_dep, here 2")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V, n_xx_emiss = 2, log_scale = TRUE),
               "n_xx_emiss should be a numeric vector with length n_dep")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V, n_xx_emiss = 2, log_scale = FALSE),
               "According to n_xx_emiss covariates have been used to predict")
  expect_error(prior_emiss_count(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                 emiss_nu = emiss_nu, emiss_V = emiss_V, n_xx_emiss = c(2,2), log_scale = TRUE),
               "According to the input paramter n_xx_emiss")
})

test_that("outputs with log_scale TRUE or FALSE",{
  expect_equal(manual_prior_emiss1,manual_prior_emiss2)
})

test_that("output dim prior_emiss_count", {
  expect_equal(length(manual_prior_emiss1), 6)
  expect_equal(as.numeric(sapply(manual_prior_emiss1, length)), c(2, n_dep, n_dep, n_dep, n_dep, 0))
  expect_equal(sapply(manual_prior_emiss1$emiss_mu0, length), rep(m, n_dep))
  expect_equal(sapply(manual_prior_emiss1$emiss_K0, length), rep(1, n_dep))
  expect_equal(sapply(manual_prior_emiss1$emiss_nu, length), rep(1, n_dep))
  expect_equal(sapply(manual_prior_emiss1$emiss_V, length), rep(m, n_dep))
  expect_equal(manual_prior_emiss3$n_xx_emiss, n_xx_emiss)
  expect_equal(sapply(manual_prior_emiss3$emiss_mu0, dim)[1,], sapply(n_xx_emiss, function(q) 1+q))
  expect_equal(sapply(manual_prior_emiss3$emiss_mu0, dim)[2,], rep(m, n_dep))
  expect_equal(sapply(manual_prior_emiss3$emiss_K0, dim)[,1], rep(1+n_xx_emiss[1], n_dep))
  expect_equal(sapply(manual_prior_emiss3$emiss_K0, dim)[,2], rep(1+n_xx_emiss[2], n_dep))
  expect_equal(sapply(manual_prior_emiss3$emiss_nu, length), rep(1, n_dep))
  expect_equal(sapply(manual_prior_emiss3$emiss_V, length), rep(m, n_dep))
})

test_that("using prior_emiss object in mHMM", {
  expect_failure(expect_error(mHMM(s_data = data_sim$obs,
                                   gen = list(m = m, n_dep = n_dep),
                                   start_val = c(list(gamma), emiss_distr),
                                   emiss_hyp_prior = manual_prior_emiss1,
                                   mcmc = list(J = 11, burn_in = 5),
                                   show_progress = FALSE, data_distr = "count")))
  expect_error(mHMM(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep),
                    start_val = c(list(gamma), emiss_distr),
                    emiss_hyp_prior = manual_prior_emiss_wrongm,
                    mcmc = list(J = 11, burn_in = 5),
                    show_progress = FALSE, data_distr = "count"),
               "number of states specified in m")
  expect_error(mHMM(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep),
                    start_val = c(list(gamma), emiss_distr),
                    emiss_hyp_prior = manual_prior_emiss_wrongn_dep,
                    mcmc = list(J = 11, burn_in = 5),
                    show_progress = FALSE, data_distr = "count"),
               "number of dependent variables")

})
