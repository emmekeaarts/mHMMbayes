################
### CREATE OUTPUT OBJECTS TO TEST COUNT DATA
################

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

emiss_distr_log <- lapply(emiss_distr, function(q) log(q))
emiss_V         <- list(rep(16, m), rep(4, m))
emiss_V_log     <- var_to_logvar(gen = list(m = m, n_dep = n_dep),
                                 emiss_mu = emiss_distr,
                                 var_emiss = emiss_V,
                                 byrow = FALSE)
emiss_V_log <- lapply(emiss_V_log, function(q) matrix(q, nrow = m))

# Define list of vectors of covariate values
xx_vec <- c(list(NULL),rep(list(rnorm(n, mean = 0, sd = 0.1)),2))

# Define object beta with regression coefficients for the three dependent variables
beta      <- rep(list(NULL), n_dep+1)
beta[[2]] <- matrix(c(1,-1,0), byrow = TRUE, ncol = 1)
beta[[3]] <- matrix(c(2,0,-2), byrow = TRUE, ncol = 1)

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

data_count2 <- sim_mHMM(n_t = n_t,
                        n = n,
                        data_distr = "count",
                        gen = list(m = m, n_dep = n_dep),
                        gamma = gamma,
                        xx_vec = xx_vec,
                        beta = beta,
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

manual_prior_emiss_log_cov <- prior_emiss_count(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = emiss_mu0_log_cov,
  emiss_K0 = emiss_K0_cov,
  emiss_V =  emiss_V_log,
  emiss_nu = emiss_nu,
  n_xx_emiss = c(1,1),
  log_scale = TRUE)

# Matrix with covariates
xx <- rep(list(matrix(1, ncol = 1, nrow = n)), (n_dep + 1))
for(i in 2:(n_dep + 1)){
  xx[[i]] <- cbind(xx[[i]], xx_vec[[i]])
}

out_count <- mHMM(s_data = data_count2$obs,
                  gen = list(m = m, n_dep = n_dep),
                  start_val = c(list(gamma), emiss_distr),data_distr = "count",emiss_hyp_prior = manual_prior_emiss_log,
                  mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)

out_count_cov <- mHMM(s_data = data_count2$obs,
                      gen = list(m = m, n_dep = n_dep), xx = xx,
                      start_val = c(list(gamma), emiss_distr),data_distr = "count",emiss_hyp_prior = manual_prior_emiss_log_cov,
                      mcmc = list(J = J, burn_in = burn_in), show_progress = FALSE)


####################
## TESTING
###############

test_that("class inherit", {
  expect_s3_class(out_count, c("mHMM","count"))
  expect_s3_class(out_count_cov, c("mHMM","count"))
})

test_that("print mHMM", {
  # 2 dependent var, no covariates
  expect_output(print(out_count), paste("Number of subjects:", n))
  expect_output(print(out_count), paste(J, "iterations"))
  expect_output(print(out_count), paste("likelihood over all subjects:", -701.4061))
  expect_output(print(out_count), paste("AIC over all subjects:", 1426.812))
  expect_output(print(out_count), paste("states used:", m))
  expect_output(print(out_count), paste("dependent variables used:", n_dep))

  # 2 dependent var, 1 covariate each
  expect_output(print(out_count_cov), paste("Number of subjects:", n))
  expect_output(print(out_count_cov), paste(J, "iterations"))
  expect_output(print(out_count_cov), paste("likelihood over all subjects:", -680.2284))
  expect_output(print(out_count_cov), paste("AIC over all subjects:", 1384.457))
  expect_output(print(out_count_cov), paste("states used:", m))
  expect_output(print(out_count_cov), paste("dependent variables used:", n_dep))
})

test_that("summary mHMM", {
  # 2 dependent var, no covariates
  expect_output(summary(out_count), "From state 1      0.784      0.114      0.102")
  expect_output(summary(out_count), "From state 3      0.200      0.215      0.586")
  expect_output(summary(out_count), "`observation 1`")
  expect_output(summary(out_count), "`observation 2`")

  # 2 dependent var, 1 covariate each
  expect_output(summary(out_count_cov), "From state 1      0.777      0.100      0.124")
  expect_output(summary(out_count_cov), "From state 3      0.207      0.237      0.556")
  expect_output(summary(out_count_cov), "`observation 1`")
  expect_output(summary(out_count_cov), "`observation 2`")

})

