context("var to logvar")


################
### CREATE OUTPUT OBJECTS TO TEST
################

# Define general parameters:
m       <- 3        # Number of hidden states
n_dep   <- 2        # Number of dependent variables

# correct specification
emiss_distr <- list(matrix(c(30, 70, 170), nrow = m),
                    matrix(c(7, 8, 18), nrow = m))
emiss_mu0 <- list(matrix(c(30, 70, 170), nrow = 1),
                  matrix(c(7, 8, 18), nrow = 1))
emiss_V   <- list(rep(16, m), rep(4, m))

logvar1 <- var_to_logvar(gen = list(m = m, n_dep = n_dep),
                         emiss_mu = emiss_mu0,
                         var_emiss = emiss_V,
                         byrow = TRUE)

logvar2 <- var_to_logvar(gen = list(m = m, n_dep = n_dep),
                         emiss_mu = emiss_distr,
                         var_emiss = emiss_V,
                         byrow = FALSE)

# wrong specifications
emiss_mu0_wrong1 <- emiss_mu0
emiss_mu0_wrong1[[1]][1] <- -30

emiss_mu0_wrong2 <- emiss_mu0
emiss_mu0_wrong2[[1]][3] <- 0

# wrong specifications
m_wrong <- 2
n_dep_wrong <- 3

emiss_mu0_wrong3 <- emiss_mu0[[1]]
emiss_mu0_wrong4 <- c(emiss_mu0,emiss_mu0[2])
emiss_mu0_wrong5 <- lapply(emiss_mu0, function(q) list(q))
emiss_mu0_wrong6 <- lapply(emiss_mu0, function(q) matrix(q[-1], nrow = 1))

emiss_V_wrong1 <- list(diag(1, m+1),
                       diag(2, m-1))
emiss_V_wrong2 <- list(1,2)
emiss_V_wrong3 <- c(1,2)
emiss_V_wrong4 <- emiss_V
emiss_V_wrong4[[1]][2] <- -1

####################
## TESTING
###############

test_that("errors prior_emiss_count input", {

  # wrong specification of means
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0_wrong1,
                             var_emiss = emiss_V,
                             byrow = TRUE),
               "All values in emiss_mu should be real positive numbers")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0_wrong2,
                             var_emiss = emiss_V,
                             byrow = TRUE),
               "All values in emiss_mu should be real positive numbers")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0_wrong3,
                             var_emiss = emiss_V,
                             byrow = TRUE),
               "emiss_mu should be a list containing")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0_wrong4,
                             var_emiss = emiss_V,
                             byrow = TRUE),
               "emiss_mu should be a list containing")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0_wrong5,
                             var_emiss = emiss_V,
                             byrow = TRUE),
               "emiss_mu should be a list containing")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0_wrong6,
                             var_emiss = emiss_V,
                             byrow = TRUE),
               "The number of columns of the emission distribution matrix in element")

  # wrong specification of variances
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0,
                             var_emiss = emiss_V_wrong1,
                             byrow = TRUE),
               "var_emiss should be a list containing n_dep")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0,
                             var_emiss = emiss_V_wrong2,
                             byrow = TRUE),
               "var_emiss should be a list containing n_dep")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0,
                             var_emiss = emiss_V_wrong3,
                             byrow = TRUE),
               "var_emiss should be a list containing n_dep")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0,
                             var_emiss = emiss_V_wrong4,
                             byrow = TRUE),
               "All values in var_emiss should be real non-negative numbers")

  # wrong specification of byrow (related to mean specification)
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_distr,
                             var_emiss = emiss_V,
                             byrow = TRUE),
               "The number of columns of the emission distribution matrix in element")
  expect_error(var_to_logvar(gen = list(m = m, n_dep = n_dep),
                             emiss_mu = emiss_mu0,
                             var_emiss = emiss_V,
                             byrow = FALSE),
               "The number of rows of the emission distribution matrix in element")

})

test_that("emiss_mu by rows or columns",{
  expect_equal(logvar1,logvar2)
})

test_that("output dim var_to_logvar", {
  expect_equal(length(logvar1), n_dep)
  expect_equal(sapply(logvar1, length), rep(m, n_dep))
})
