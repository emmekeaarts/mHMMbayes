head(nonverbal)
m <-2
n_dep <- 4
q_emiss <- c(3, 2, 3, 2)

# specifying starting values
start_TM <- diag(.8, m)
start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM <- list(matrix(c(0.05, 0.90, 0.05, 0.90, 0.05, 0.05), byrow = TRUE,
                        nrow = m, ncol = q_emiss[1]), # vocalizing patient
                 matrix(c(0.1, 0.9, 0.1, 0.9), byrow = TRUE, nrow = m,
                        ncol = q_emiss[2]), # looking patient
                 matrix(c(0.90, 0.05, 0.05, 0.05, 0.90, 0.05), byrow = TRUE,
                        nrow = m, ncol = q_emiss[3]), # vocalizing therapist
                 matrix(c(0.1, 0.9, 0.1, 0.9), byrow = TRUE, nrow = m,
                        ncol = q_emiss[4])) # looking therapist

# Run a model without covariate(s):
set.seed(14532)
out1 <- mHMM(s_data = nonverbal,
             gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
             start_val = c(list(start_TM), start_EM),
             mcmc = list(J = 1000, burn_in = 200))

# for this dataset, we use:
# PD_subj
# gamma_prob_bar
# emiss_prob_bar

# check, still delete not used slots

save(out1, file = "nonv_2st_1000it.rda")


## different starting values
m <-2
n_dep <- 4
q_emiss <- c(3, 2, 3, 2)

# specifying different starting values
start_TM <- diag(.8, m)
start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM_b <- list(matrix(c(0.2, 0.6, 0.2,
                            0.6, 0.2, 0.2), byrow = TRUE,
                          nrow = m, ncol = q_emiss[1]), # vocalizing patient
                   matrix(c(0.4, 0.6,
                            0.4, 0.6), byrow = TRUE, nrow = m,
                          ncol = q_emiss[2]), # looking patient
                   matrix(c(0.6, 0.2, 0.2,
                            0.2, 0.6, 0.2), byrow = TRUE,
                          nrow = m, ncol = q_emiss[3]), # vocalizing therapist
                   matrix(c(0.4, 0.6,
                            0.4, 0.6), byrow = TRUE, nrow = m,
                          ncol = q_emiss[4])) # looking therapist

set.seed(9843)
out1b <- mHMM(s_data = nonverbal,
              gen = list(m = m, n_dep = n_dep,q_emiss = q_emiss),
              start_val = c(list(start_TM), start_EM_b),
              mcmc = list(J = 1000, burn_in = 200))

out1b$PD_subj <- NULL
out1b$emiss_cov_bar <-NULL
out1b$emiss_int_bar <-NULL
out1b$emiss_int_subj <-NULL
out1b$emiss_naccept <-NULL
out1b$gamma_cov_bar <-NULL
out1b$gamma_int_bar <-NULL
out1b$gamma_int_subj <-NULL
out1b$gamma_naccept <-NULL

save(out1b, file = "nonv_2stb_1000it.rda")



#### 3 states
m <-3

# specifying starting values
start_TM <- diag(.8, m)
start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM <- list(matrix(c(0.05, 0.90, 0.05,
                          0.15, 0.75, 0.10,
                          0.90, 0.05, 0.05), byrow = TRUE,
                        nrow = m, ncol = q_emiss[1]), # vocalizing patient
                 matrix(c(0.1, 0.9,
                          0.1, 0.9,
                          0.1, 0.9), byrow = TRUE, nrow = m,
                        ncol = q_emiss[2]), # looking patient
                 matrix(c(0.90, 0.05, 0.05,
                          0.15, 0.75, 0.10,
                          0.05, 0.90, 0.05), byrow = TRUE,
                        nrow = m, ncol = q_emiss[3]), # vocalizing therapist
                 matrix(c(0.1, 0.9,
                          0.1, 0.9,
                          0.1, 0.9), byrow = TRUE, nrow = m,
                        ncol = q_emiss[4])) # looking therapist

# Run a model without covariate(s):
set.seed(24141)
out2 <- mHMM(s_data = nonverbal,
             gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
             start_val = c(list(start_TM), start_EM),
             mcmc = list(J = 1000, burn_in = 200))

# for this data, we only use the AIC, no real data, do not include in package!

save(out2, file = "nonv_3st_1000it.rda")


### 4 states
m <-4

# specifying starting values
start_TM <- diag(.8, m)
start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM <- list(matrix(c(0.05, 0.90, 0.05,
                          0.15, 0.75, 0.10,
                          0.90, 0.05, 0.05,
                          0.90, 0.05, 0.05), byrow = TRUE,
                        nrow = m, ncol = q_emiss[1]), # vocalizing patient
                 matrix(c(0.1, 0.9,
                          0.1, 0.9,
                          0.1, 0.9,
                          0.1, 0.9), byrow = TRUE, nrow = m,
                        ncol = q_emiss[2]), # looking patient
                 matrix(c(0.90, 0.05, 0.05,
                          0.15, 0.75, 0.10,
                          0.05, 0.05, 0.90,
                          0.05, 0.05, 0.90), byrow = TRUE,
                        nrow = m, ncol = q_emiss[3]), # vocalizing therapist
                 matrix(c(0.1, 0.9,
                          0.1, 0.9,
                          0.1, 0.9,
                          0.8, 0.2), byrow = TRUE, nrow = m,
                        ncol = q_emiss[4])) # looking therapist

# Run a model without covariate(s):
set.seed(65342)
out3 <- mHMM(s_data = nonverbal,
             gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),
             start_val = c(list(start_TM), start_EM),
             mcmc = list(J = 1000, burn_in = 200))

# for this dataset, we only use:
# gamma_int_bar
# emiss_int_bar

out3$emiss_cov_bar <-NULL
out3$emiss_prob_bar <-NULL
out3$emiss_int_subj <-NULL
out3$emiss_naccept <-NULL
out3$gamma_cov_bar <-NULL
out3$gamma_prob_bar <-NULL
out3$gamma_int_subj <-NULL
out3$gamma_naccept <-NULL

save(out3, file = "nonv_4st_1000it.rda")

