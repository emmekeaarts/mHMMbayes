#' Simulate data using a multilevel hidden markov model
#'
#' This function simulates data for multiple subjects, for which the data have
#' categorical observations that follow a hidden Markov model (HMM) with an
#' multilevel structure. The multilevel structure means that each subject can
#' have it's own set of parameters, and that the parameters at the subject
#' level (level 1) are tied together by a population distribution at level 2 for
#' each of the corresponding parameters. The shape of the population distribution
#' for each of the parameters is a normal (i.e., Gaussian) distribution. In
#' simulating the data, having a multilevel structure means that the parameters
#' for each subject are sampled from the population level distribution of the
#' corresponding parameter. The user specifies the population distribution for
#' each parameter: the average population transition probability matrix and its
#' variance, and the average population emission distribution and its variance.
#' For now, the variance is assumed fixed for all components of the transition
#' probability matrix and for all components of the emissiondistribution, and
#' the simulated data can only consist of one dependent variable.
#'
#' @param T The length of the observed sequence to be simulated for each
#'   subject. To only simulate subject specific transition probability matrices
#'   gamma and emission distributions (and no data), set \code{T} to 0.
#' @param n The number of subjects for which data is simulated.
#' @param m The number of hidden states in the HMM used for
#'   simulating data.
#' @param pr The number of categories of the simulated observations.
#' @param gamma A matrix with \code{m} rows and \code{m} columns containing the
#'   average population transition probability matrix used for simulating the
#'   data. That is, the probability to switch from hidden state \emph{i} (row
#'   \emph{i}) to hidden state \emph{j} (column  \emph{j}).
#' @param emiss_distr A matrix with \code{m} rows and \code{pr} columns
#'   containing the average population emission distribution of the
#'   (categorical) observations given the hidden states. That is, the
#'   probability of observing category \emph{k} (column \emph{k}) in state
#'   \emph{i} (row \emph{i}).
#' @param var_gamma An integer denoting the variance between subjects in the
#'   transition probability matrix. Note that this value corresponds to the
#'   variance of the parameters of the multinomial distribution (i.e., the
#'   intercepts of the regression equation of the multinomial distribution used
#'   to sample the transition probability matrix), see details below. In
#'   addition, only one variance value can be specified for the complete
#'   transition probability matrix, hence the variance is assumed fixed across
#'   all components. The default equals 1, which corresponds to quite some
#'   variation between subjects. A less extreme value would be 0.5. If one wants
#'   to simulate data from exactly the same HMM for all subjects,
#'   \code{var_gamma} should be set to 0.
#' @param var_emiss An integer denoting the variance between subjects in the
#'   emission distribution. Note that this value corresponds to the variance of
#'   the parameters of the multinomial distribution (i.e., the intercepts of the
#'   regression equation of the multinomial distribution used to sample the
#'   components of the emission distribution), see details below.  In addition,
#'   only one variance value can be specified for the complete emission
#'   distribution, hence the variance is assumed fixed across all components.
#'   The default equals 1, which corresponds to quite some variation between
#'   subjects. A less extreme value would be 0.5. If one wants to simulate data
#'   from exactly the same HMM for all subjects, \code{var_emiss} should be set
#'   to 0.
#' @param return_ind_par A logical scalar. Should the subject specific
#'   transition probability matrix gamma and emission probability matrix be
#'   returned by the function (\code{return_ind_par = TRUE}) or not
#'   (\code{return_ind_par = FALSE}). The default equals \code{return_ind_par =
#'   FALSE}.
#'
#' Note: the subject specific) initial state distributions (i.e., the
#' probability of each of the states at the first time point) needed to simulate
#' the data are obtained from the stationary distributions of the subject
#' specific transition probability matrices gamma.

# FUTURE:
# - add saving the subject specific gamma and emmis prob as well.
# - multivariate, so multiple dependent variables
# - covariate to predict gamma / emmis prob
# - possibility to vary variance per state / component

#' @examples
#' # simulating data for 10 subjects with each 100 observations
#' T <- 100
#' n <- 10
#' m <- 3
#' pr <- 4
#' gamma <- matrix(c(0.8, 0.1, 0.1,
#'                   0.2, 0.7, 0.1,
#'                   0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#' emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
#'                         0.1, 0.1, 0.8, 0.0,
#'                         0.0, 0.0, 0.1, 0.9), nrow = m, ncol = pr, byrow = TRUE)
#' data1 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
#'                   var_gamma = 1, var_emiss = 1)
#'
#'
#' # simulating subject specific transition probability matrices and emission distributions only
#' T <- 0
#' n <- 5
#' m <- 3
#' pr <- 4
#' gamma <- matrix(c(0.8, 0.1, 0.1,
#'                   0.2, 0.7, 0.1,
#'                   0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
#' emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
#'                         0.1, 0.1, 0.8, 0.0,
#'                         0.0, 0.0, 0.1, 0.9), nrow = m, ncol = pr, byrow = TRUE)
#' data2 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
#'                   var_gamma = 1, var_emiss = 1)
#' data2
#'
#' data3 <- sim_mHMM(T = T, n = n, m = m, pr = pr, gamma = gamma, emiss_distr = emiss_distr,
#'                   var_gamma = .5, var_emiss = .5)
#' data3
#'
#'
#'
#' @export

sim_mHMM <- function(T, n, m, pr, gamma, emiss_distr, var_gamma = 1, var_emiss = 1, return_ind_par = FALSE){
  if (dim(gamma)[1] != m){
    stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
  }
  if (dim(gamma)[2] != m){
    stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
  }
  if (dim(emiss_distr)[1] != m){
    stop(paste("The number of rows of the emission distribution matrix should be equal to the number of states, which is", m, "."))
  }
  if (dim(emiss_distr)[2] != pr){
    stop(paste("The number of columns of the emission distribution matrix should be equal to the number of observable categories, which is", pr, "."))
  }
  states <- matrix(ncol = 2, nrow = T*n)
  states[,1] <- rep(1:n, each = T)
  obs <- matrix(ncol = 2, nrow = T*n)
  obs[,1] <- rep(1:n, each = T)
  sub_gamma <- rep(list(NULL), n)
  sub_emiss <- rep(list(NULL), n)
  mnl_gamma <- prob_to_int(gamma)
  mnl_emiss <- prob_to_int(emiss_distr)
  for(j in 1:n){
    sub_gamma[[j]] <- int_to_prob(mnl_gamma +
                                    rnorm(n = m * (m-1), mean = 0, sd = sqrt(var_gamma)))
    sub_emiss[[j]] <- int_to_prob(mnl_emiss +
                                    rnorm(n = m * (pr-1), mean = 0, sd = sqrt(var_emiss)))
    if(T != 0){
      init <- solve(t(diag(m) - sub_gamma[[j]] + 1), rep(1, m))
      states[((j-1) * T + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
      obs[((j-1) * T + 1), 2] <- sample(x = 1:pr, size = 1, prob = sub_emiss[[j]][states[((j-1) * T + 1), 2],])
      for(i in 2:T){
        states[((j-1) * T + i), 2] <- sample(x = 1:m, size = 1, prob = sub_gamma[[j]][states[((j-1) * T + i - 1), 2],])
        obs[((j-1) * T + i), 2] <- sample(x = 1:pr, size = 1, prob = sub_emiss[[j]][states[((j-1) * T + i), 2],])
      }
    }
  }
  if (return_ind_par == FALSE & T != 0){
    return(list(states = states, obs = obs))
  } else if (return_ind_par == TRUE & T != 0){
    return(list(states = states, obs = obs, subject_gamma = sub_gamma, subject_emmis = sub_emiss))
  } else if (T == 0){
    return(list(subject_gamma = sub_gamma, subject_emmis = sub_emiss))
  }
}
