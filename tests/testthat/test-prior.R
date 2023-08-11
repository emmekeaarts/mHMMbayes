context("prior and RW sampler functions")


################
### CREATE OUTPUT OBJECTS TO TEST
################

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
data_sim <- sim_mHMM(n_t = n_t, n = n, gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss), gamma = gamma,
                     emiss_distr = emiss_distr, var_gamma = .5, var_emiss = c(.5, 0.5))
colnames(data_sim$obs) <- c("subj", "output_1", "output_2")

# hypothesized mean emission probabilities
prior_prob_emiss_cat <- list(matrix(c(0.40, 0.40, 0.10, 0.10,
                                      0.10, 0.10, 0.70, 0.10,
                                      0.10, 0.10, 0.10, 0.70), byrow = TRUE,
                                    nrow = m, ncol = q_emiss[1]),
                             matrix(c(0.60, 0.40,
                                      0.80, 0.20,
                                      0.70, 0.30), byrow = TRUE, nrow = m,
                                    ncol = q_emiss[2]))

# correct specifications
prior_int_emiss <- sapply(prior_prob_emiss_cat, prob_to_int)
emiss_mu0 <- rep(list(vector(mode = "list", length = m)), n_dep)
for(k in 1:n_dep){
  for(i in 1:m){
    emiss_mu0[[k]][[i]] <- matrix(prior_int_emiss[[k]][i,], nrow = 1)
  }
}

emiss_K0 <- rep(list(c(1)), n_dep)
emiss_nu <- list(c(5), c(4))
emiss_V <- list(diag(5, q_emiss[1] - 1),
                diag(4, q_emiss[2] - 1))

manual_prior_emiss <- prior_emiss_cat(gen = list(m = m,
                                                 n_dep = n_dep,
                                                 q_emiss = q_emiss),
                                      emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                      emiss_nu = emiss_nu, emiss_V = emiss_V)

# wrong specifications
m_wrong <- 2
n_dep_wrong <- 3
q_emiss_wrong <- c(2,3,2)

emiss_mu0_wrong1 <- emiss_mu0[[1]]
emiss_mu0_wrong2 <- rep(list(vector(mode = "list", length = m)), n_dep)
for(k in 1:n_dep){
  for(i in 1:m){
    emiss_mu0_wrong2[[k]][[i]] <- matrix(prior_int_emiss[[1]][i,], nrow = 1)
  }
}
emiss_mu0_wrong3 <- rep(list(vector(mode = "list", length = m)), n_dep)
for(k in 1:n_dep){
  for(i in 1:m){
    emiss_mu0_wrong3[[k]][[i]] <- matrix(prior_int_emiss[[2]][i,], nrow = 1)
  }
}
emiss_mu0_wrong4 <- prior_int_emiss

emiss_K0_wrong1 <- rep(list(c(1)), n_dep + 1)
emiss_K0_wrong2 <- rep((c(1)), n_dep)
emiss_nu_wrong1 <- list(c(5), c(4), c(5))
emiss_nu_wrong2 <- (c(5, 4, 5))
emiss_V_wrong1 <- list(diag(5, q_emiss[1]),
                diag(4, q_emiss[2]))

manual_prior_emiss_wrongm <- manual_prior_emiss
manual_prior_emiss_wrongm$gen$m <- 2
manual_prior_emiss_wrongn_dep <- manual_prior_emiss
manual_prior_emiss_wrongn_dep$gen$n_dep <- 3
manual_prior_emiss_wrongq_emiss1 <- manual_prior_emiss
manual_prior_emiss_wrongq_emiss1$gen$q_emiss <- c(3,3)
manual_prior_emiss_wrongq_emiss2 <- manual_prior_emiss
manual_prior_emiss_wrongq_emiss2$gen$q_emiss <- c(4,2,3)

# representing a prior belief that switching to state 3 does not occur often and
# state 3 has a relative short duration
prior_prob_gamma <- matrix(c(0.70, 0.25, 0.05,
                             0.25, 0.70, 0.05,
                             0.30, 0.30, 0.40), nrow = m, ncol = m, byrow = TRUE)

# using the function prob_to_int to obtain intercept values for the above specified
# transition probability matrix gamma
prior_int_gamma <- prob_to_int(prior_prob_gamma)
gamma_mu0 <- list(matrix(prior_int_gamma[1,], nrow = 1, ncol = m-1),
                  matrix(prior_int_gamma[2,], nrow = 1, ncol = m-1),
                  matrix(prior_int_gamma[3,], nrow = 1, ncol = m-1))

gamma_K0 <- 1
gamma_nu <- 5
gamma_V <- diag(5, m - 1)

manual_prior_gamma <- prior_gamma(m = m, gamma_mu0 = gamma_mu0,
                                  gamma_K0 = gamma_K0, gamma_nu = gamma_nu,
                                  gamma_V = gamma_V)

gamma_mu0_wrong1 <- prior_int_gamma
gamma_mu0_wrong2 <- list(prior_int_gamma)
gamma_mu0_wrong3 <- gamma_mu0
gamma_mu0_wrong3[[1]] <- c(1)
gamma_mu0_wrong4 <- gamma_mu0
gamma_mu0_wrong4[[1]] <- matrix(1, ncol = 1, nrow = 1)

gamma_K0_wrong1 <- c(1,1)
gamma_K0_wrong2 <- list(c(1))
gamma_nu_wrong1 <- c(1,1)
gamma_nu_wrong2 <- list(c(1))
gamma_V_wrong1 <- diag(5, m)
gamma_V_wrong2 <- c(5)

manual_prior_gamma_wrongm <- manual_prior_gamma
manual_prior_gamma_wrongm$m <- 2

####################
## TESTING
###############

test_that("errors prior_emiss input", {
  expect_error(prior_emiss_cat(gen = list(m = m_wrong, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_mu0 should consist of m , here 2")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep_wrong, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "should equal the number of dependent variables specified in n_dep")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss_wrong),
                               emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "should equal the number of dependent variables specified in n_dep")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0_wrong1, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_mu0 should be a list of lists containing 2 lists")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0_wrong2, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "For the list relating to dependent variable 2 of the input argument emiss_mu0")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0_wrong3, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "For the list relating to dependent variable 1 of the input argument emiss_mu0")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0_wrong4, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_mu0 should be a list of lists containing 2 lists")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0_wrong1,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_K0 should be a list containing n_dep, here 2")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0_wrong2,
                               emiss_nu = emiss_nu, emiss_V = emiss_V),
               "emiss_K0 should be a list containing n_dep, here 2")
  expect_warning(expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu_wrong1, emiss_V = emiss_V),
               "emiss_nu should be a list containing n_dep"))
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                                              emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                                              emiss_nu = emiss_nu_wrong2, emiss_V = emiss_V),
                              "emiss_nu should be a list containing n_dep")
  expect_error(prior_emiss_cat(gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                               emiss_mu0 = emiss_mu0, emiss_K0 = emiss_K0,
                               emiss_nu = emiss_nu, emiss_V = emiss_V_wrong1),
               "emiss_V should be a list containing n_dep, here 2")
})

test_that("output dim prior_emiss", {
  expect_equal(length(manual_prior_emiss), 6)
  expect_equal(as.numeric(sapply(manual_prior_emiss, length)), c(3, n_dep, n_dep, n_dep, n_dep, 0))
  expect_equal(sapply(manual_prior_emiss$emiss_mu0, length), rep(m, n_dep))
  expect_equal(sapply(manual_prior_emiss$emiss_mu0[[1]], dim), matrix(c(1, q_emiss[1] - 1), ncol = m, nrow = 2))
  expect_equal(sapply(manual_prior_emiss$emiss_mu0[[2]], dim), matrix(c(1,  q_emiss[2] - 1), ncol = m, nrow = 2))
  expect_equal(sapply(manual_prior_emiss$emiss_K0, length), rep(1, n_dep))
  expect_equal(sapply(manual_prior_emiss$emiss_nu, length), rep(1, n_dep))
  expect_equal(sapply(manual_prior_emiss$emiss_V, dim), matrix(c(rep(q_emiss[1] - 1, 2),
                                                                 rep(q_emiss[2] - 1, 2)), ncol = n_dep, nrow = 2))
})

test_that("using prior_emiss object in mHMM", {
  expect_failure(expect_error(mHMM(s_data = data_sim$obs,
                                   gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                                   start_val = c(list(gamma), emiss_distr),
                                   emiss_hyp_prior = manual_prior_emiss,
                                   mcmc = list(J = 11, burn_in = 5))))
  expect_error(mHMM(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    start_val = c(list(gamma), emiss_distr),
                    emiss_hyp_prior = manual_prior_emiss_wrongm,
                    mcmc = list(J = 11, burn_in = 5)),
               "number of states specified in m")
  expect_error(mHMM(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    start_val = c(list(gamma), emiss_distr),
                    emiss_hyp_prior = manual_prior_emiss_wrongn_dep,
                    mcmc = list(J = 11, burn_in = 5)),
               "number of dependent variables")
  expect_error(mHMM(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    start_val = c(list(gamma), emiss_distr),
                    emiss_hyp_prior = manual_prior_emiss_wrongq_emiss1,
                    mcmc = list(J = 11, burn_in = 5)),
               "number of number of observed categories")
  expect_warning(expect_error(mHMM(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    start_val = c(list(gamma), emiss_distr),
                    emiss_hyp_prior = manual_prior_emiss_wrongq_emiss2,
                    mcmc = list(J = 11, burn_in = 5)),
               "number of number of observed categories"))

})

test_that("errors prior_gamma input", {
  expect_error(prior_gamma(m = m_wrong,
                               gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0,
                               gamma_nu = gamma_nu, gamma_V = gamma_V),
               "gamma_mu0 should be a list containing 2 matrices")
  expect_error(prior_gamma(m = m,
                               gamma_mu0 = gamma_mu0_wrong1, gamma_K0 = gamma_K0,
                               gamma_nu = gamma_nu, gamma_V = gamma_V),
               "gamma_mu0 should be a list containing 3 matrices")
  expect_error(prior_gamma(m = m,
                               gamma_mu0 = gamma_mu0_wrong2, gamma_K0 = gamma_K0,
                               gamma_nu = gamma_nu, gamma_V = gamma_V),
               "gamma_mu0 should be a list containing 3 matrices")
  expect_error(prior_gamma(m = m,
                               gamma_mu0 = gamma_mu0_wrong3, gamma_K0 = gamma_K0,
                               gamma_nu = gamma_nu, gamma_V = gamma_V),
               "gamma_mu0 should be a list containing 3 matrices")
  expect_error(prior_gamma(m = m,
                               gamma_mu0 = gamma_mu0_wrong4, gamma_K0 = gamma_K0,
                               gamma_nu = gamma_nu, gamma_V = gamma_V),
               "Each matrix in the list gamma_mu0")
  expect_error(prior_gamma(m = m,
                               gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0_wrong1,
                               gamma_nu = gamma_nu, gamma_V = gamma_V),
               "gamma_K0 should be a numeric vector with length 1")
  expect_error(prior_gamma(m = m,
                               gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0_wrong2,
                               gamma_nu = gamma_nu, gamma_V = gamma_V),
               "gamma_K0 should be a numeric vector with length 1")
  expect_error(prior_gamma(m = m,
                          gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0,
                          gamma_nu = gamma_nu_wrong1, gamma_V = gamma_V),
              "gamma_nu should be a numeric vector with length 1")
  expect_error(prior_gamma(m = m,
                          gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0,
                          gamma_nu = gamma_nu_wrong2, gamma_V = gamma_V),
              "gamma_nu should be a numeric vector with length 1")
  expect_error(prior_gamma(m = m,
                               gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0,
                               gamma_nu = gamma_nu, gamma_V = gamma_V_wrong1),
               "gamma_V should be an 2 by 2 matrix")
  expect_error(prior_gamma(m = m,
                           gamma_mu0 = gamma_mu0, gamma_K0 = gamma_K0,
                           gamma_nu = gamma_nu, gamma_V = gamma_V_wrong2),
               "gamma_V should be an 2 by 2 matrix")
})

test_that("output dim prior_gamma", {
  expect_equal(length(manual_prior_gamma), 6)
  expect_equal(as.numeric(sapply(manual_prior_gamma, length)), c(1, 0, m, 1, 1, (m-1) *( m-1)))
  expect_equal(sapply(manual_prior_gamma$gamma_mu0, dim), matrix(c(1, m-1), nrow = 2, ncol = m))
  expect_equal(dim(manual_prior_gamma$gamma_K0), c(1,1))
  expect_equal(length(manual_prior_gamma$gamma_nu), c(1))
  expect_equal(dim(manual_prior_gamma$gamma_V), c(m-1, m-1))
})

test_that("using prior_gamma object in mHMM", {
  expect_failure(expect_error(mHMM(s_data = data_sim$obs,
                                   gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                                   start_val = c(list(gamma), emiss_distr),
                                   gamma_hyp_prior = manual_prior_gamma,
                                   mcmc = list(J = 11, burn_in = 5))))
  expect_error(mHMM(s_data = data_sim$obs,
                    gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
                    start_val = c(list(gamma), emiss_distr),
                    gamma_hyp_prior = manual_prior_gamma_wrongm,
                    mcmc = list(J = 11, burn_in = 5)),
               "number of states specified in m")
})
